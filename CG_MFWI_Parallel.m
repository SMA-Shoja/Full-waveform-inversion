% Time-domain multi-scale Full Waveform Inversion with conjugate gradient
% algorithm
% AUTHOR: AYDIN SHOJA, University of Tehran 2018

clc
clear
close all
tic
delete(gcp('nocreate'))
parpool('local',16)
%% loading model
load('marmhard.dat');
v_true = reshape(marmhard,122,384);
v_true = v_true(1:2:end,1:2:end);    % reduce model dimensions
[nz,nx] = size(v_true);
h = 20;                              % model grid length

%% Temporal parameters
t = 3;
dt = .0004;
nt = round(t/dt);
k = 1:nt;
T = k.*dt;
% fd =30;
f = (1/dt)*(0:nt/2)/nt;

cfl = ((max(v_true(:))*dt)/h)*sqrt(2);  % stability factor
if cfl >= 1
    error('Unstable')
end
%% sources locations
ns_up = 3:5:nx-3;               % upper boundary
ms_up = 2*ones(1,length(ns_up));
ms_right = 3:4:nz-3;            % right side boundary
ns_right = (nx-2)*ones(1,length(ms_right));
ns_down = 10:10:nx-3;           % lower boundary
ms_down = nz-5*ones(1,length(ns_down));
ms_left = 3:7:nz-3;             % left side boundary
ns_left = 2*ones(1,length(ms_left));
ms = [ms_up];%ms_up ms_right ms_left] ;   % uncomment any boundary you desire
ns = [ns_up];%ns_up ns_right ns_left];
ks = ((ns-1)*nz)+ms; % sources locations index in model

%% recievers locations
G = 3; % recivers interval
nr_up = G:1:nx-3;
mr_up = 2*ones(1,length(nr_up));
mr_right = G:1:nz-3;
nr_right = (nx-2)*ones(1,length(mr_right));
nr_down = nx-3:-1:G;
mr_down = nz-3*ones(1,length(nr_down));
mr_left = G:1:nz-3;
nr_left = 2*ones(1,length(mr_left));
mr = [mr_up];% mr_down mr_right];
nr = [nr_up];% nr_down nr_right];
kr = ((nr-1)*nz)+mr; % recivers locations index in model


%% initial model
v = velsmooth(v_true,20,20,1000);  % smoothing true velocity to creat an initial model
imagesc(v)
v_initial = v;
%% Wavelete definition
fd = 20;
fd_old = fd;
Source = 100*( 1-2.*pi.^2.*fd.^2.*(k.*dt-(0.5)).^2 ).*exp( -pi.^2.*fd.^2.*(k.*dt-(0.5)).^2 );
Source= Source';
%% data computation
disp(' observation data construction ')
d_obs = data2Dpartest(nx,nz,h,v_true,ks,kr,Source,t,dt,0);   % observed data
d_cal = data2Dpartest(nx,nz,h,v,ks,kr,Source,t,dt,0);        % calculated data

r_initial = d_obs - d_cal;                                   % initial residuals
[n,m] = size(d_obs);

%% adding noise
noise = randn(n,m);
d_obsn = d_obs + 0.0*noise;
SNR = snr(d_obs,0.0*noise)

%% lowpass filter
fdlp = 5;            % lowpass filter dominant frequency
if fd ~= fdlp
    targetwavelet = 100*( 1-2.*pi.^2.*fdlp.^2.*(k.*dt-(0.5)).^2 ).*exp( -pi.^2.*fdlp.^2.*(k.*dt-(0.5)).^2 );
    targetwavelet =  targetwavelet';
    Source2= wienerfilter(Source,targetwavelet,Source);
    d_obs2= wienerfilter(Source,targetwavelet,d_obsn);
else
    d_obs2 = d_obsn;
    Source2 = Source;
end
Sourcelp = Source2;



%% FWI loop initialization
e_old= [];
v_old = []'
e_normal = 1;
iter = 1
e_new = 1;
eta = 0.9;
grad_old_vec = zeros(nz*nx,1);
p_old_vec = zeros(nz*nx,1);
alp = [];
e_normal_old =[];
Model_rms_initial = sqrt(sum(((v_true(:) - v(:)).^2)./(nx*nz)));
model_rms_old = [];
model_rms_normal_old = [];
[Dh,Dv] = drivop(nz,nx,2);
L = ((Dh./(h^2))'*(Dh./(h^2))) + ((Dv./(h^2))'*(Dv./(h^2)));
i = 0;

%% FWI loop
while iter <= 100
    
    iter_vec = 1:iter;
    
    disp(' calculated data construction ')
    [d_cal,sw] = data2Dpartest(nx,nz,h,v,ks,kr,Source2,t,dt,0);  % sw: forward wavefield sampled at each model grid for all shots
    
    
    %% residual and misfit
    r = d_obs2 - d_cal;

    e_new = r(:)'*r(:) %+ lam*v(:)'*L*v(:); 
    
    
    
    model_rms_new = sqrt(sum(((v_true(:) - v(:)).^2)./(nx*nz)));
    
    %% multi-scaling
    i_old = i;
    if iter > 1
        if  i == 0 && ((iter > 25 && fd ~= fdlp) || (abs(e_new - e(iter-1)) < 0.0001 ||(e_new > 1.1*e(iter-1))))
            disp ( ' frequency change ' )
          
            i = i + 1
            fdlp = 10;
            targetwavelet = 100*( 1-2.*pi.^2.*fdlp.^2.*(k.*dt-(0.5)).^2 ).*exp( -pi.^2.*fdlp.^2.*(k.*dt-(0.5)).^2 );
            targetwavelet =  targetwavelet';
            Source2= wienerfilter(Source,targetwavelet,Source);
            d_obs2= wienerfilter(Source,targetwavelet,d_obs);
            [d_cal,sw] = data2Dpartest(nx,nz,h,v,ks,kr,Source2,t,dt,0);
            r = d_obs2 - d_cal;
            e_new = r(:)'* r(:) %+ lam*v(:)'*L*v(:);
        elseif i == 1 && ((iter > 50 && fd ~= fdlp) || (abs(e_new - e(iter-1)) < 0.0001 ||(e_new > 1.1*e(iter-1))))
            i = i + 1
            fdlp = 15;
            targetwavelet = 100*( 1-2.*pi.^2.*fdlp.^2.*(k.*dt-(0.2)).^2 ).*exp( -pi.^2.*fdlp.^2.*(k.*dt-(0.2)).^2 );
            targetwavelet =  targetwavelet';
            Source2= wienerfilter(Source,targetwavelet,Source);
            d_obs2= wienerfilter(Source,targetwavelet,d_obs);
            [d_cal,sw] = data2Dpartest(nx,nz,h,v,ks,kr,Source2,t,dt,0);
            r = d_obs2 - d_cal;
            e_new = r(:)'* r(:) %+ lam*v(:)'*L*v(:);
        elseif i == 2 && ((iter > 75 && fd ~= fdlp) || (abs(e_new - e(iter-1)) < 0.0001 ||(e_new > 1.1*e(iter-1))))
            i = i+1
            disp ( ' main frequency band ' )
            d_obs2 = d_obsn;
            Source2 = Source;
            [d_cal,sw] = data2Dpartest(nx,nz,h,v,ks,kr,Source2,t,dt,0);
            r = d_obs2 - d_cal;
            e_new = r(:)'* r(:) %+ lam*v(:)'*L*v(:);
        elseif fd == fdlp && i == 0
            i = i+3;
        elseif (abs(e_new - e(iter-1)) < 0.0001 ||(e_new > 1.1*e(iter-1) && e_new < 1.2*e(iter-1))) && i == 3
            disp( ' minimum reached ' )
            break;
            
        elseif  isnan(e_new)
            disp ( 'NaN' )
            v = v_olds;
            break;
        elseif e_new >= 1.2*e(iter-1)
            disp ( 'High' )
            break;
        elseif e_new < 0.01*e_initial
            disp ( 'Low' )
            break;
        end
    end
    
    e = [e_old e_new];
    e_old = e;
    
    if iter == 1
        e_initial = e(1);
    elseif i ~= i_old
        e_initial = e_new;
    end
    
    
    e_normal_new = e(iter)./e_initial;
    e_normal = [e_normal_old e_normal_new];
    e_normal_old = e_normal;
    
    Model_rms = [model_rms_old model_rms_new];
    model_rms_old = Model_rms;
    model_rms_normal_new = Model_rms(iter)./Model_rms_initial;
    Model_rms_normal = [model_rms_normal_old model_rms_normal_new];
    model_rms_normal_old = Model_rms_normal;
    
    %% trace camparison
    figure(1)
    subplot(1,2,1)
    plot(e_normal,'--x')
    title('normalized cost function')
    ylabel( ' normalized L2 norm of residuals ' )
    xlabel( ' Iteration ')
    axis square
    pause(1)
    
    subplot(1,2,2)
    plot(Model_rms_normal,'--o')
    title('normalized Model RMS')
    xlabel( ' Iteration ')
    axis square
    pause(1)
    
    figure(3)
    subplot(3,2,1)
    plot(k,d_obs2(:,200),k,d_cal(:,200))
    title('200^{th} trace')
    axis square
    pause(1)
    
    subplot(3,2,2)
    plot(k,d_obs2(:,900),k,d_cal(:,900))
    title('900^{th} trace')
    axis square
    pause(1)
    
    subplot(3,2,3)
    plot(k,d_obs2(:,1500),k,d_cal(:,1500))
    title('1500^{th} trace')
    axis square
    pause(1)
    
    subplot(3,2,4)
    plot(k,d_obs2(:,2600),k,d_cal(:,2600))
    title('2600^{th} trace')
    axis square
    pause(1)
    
    subplot(3,2,5)
    plot(k,d_obs2(:,3400),k,d_cal(:,3400))
    title('3400^{th} trace')
    axis square
    pause(1)
    
    subplot(3,2,6)
    plot(k,d_obs2(:,500),k,d_cal(:,500))
    title('500^{th} trace')
    axis square
    pause(1)
    
    disp(' Back Propagation ')
    rsw = BackPropagationPartest(nx,nz,h,v,ks,kr,r,t,dt,0);
    
    %% gradiant of misfit function
    disp(' update direction calculation ')
    tic
    grad_pershot_vec = rsw(1:end-2,:).*((1/dt^2)*diff(sw,2,1));
    grad_pershot_vec = (2./(repmat(v(:)',[1 length(ks)]).^3)).*sum(grad_pershot_vec);
    grad = reshape(grad_pershot_vec,[nx*nz,length(ks)]);
    grad_vec = sum(grad,2);
    grad = reshape(grad_vec,[nz,nx]);
    
    %% CG algorithm
    if iter == 1
        beta = 0;
    elseif i ~= i_old
        beta = 0;
    else
        beta = (grad_vec'*(grad_vec - grad_old_vec))/((grad_vec - grad_old_vec)'*p_old_vec); % Hestenes - Steifel
        % beta = (norm(grad_vec).^2)/(norm(grad_old_vec).^2);                                % Fletcher - Reeves
        % beta = (grad_vec'*(grad_vec - grad_old_vec))/(norm(grad_vec).^2);                  % Polak - Rebier
        % beta = (norm(grad_vec)^2) / (p_old_vec'*(grad_vec - grad_old_vec));                % Dai - Yuan
    end
    
    p_vec = -grad_vec + beta*p_old_vec;
    p = reshape(p_vec,[nz,nx]);
    %     p(1:10,:) = 0;
    
    grad_old_vec = grad_vec;
    p_old_vec = p_vec;
    toc
    
    disp(' evaluating step length ')
    tic
    
    alpha_t = 5e12;
    while max(max(abs(alpha_t*p))) > max(max(v))/100
        alpha_t = 0.9 * alpha_t;
    end
    vv = v + alpha_t*p;
    d_cal_new = data2Dpartest(nx,nz,h,vv,ks,kr,Source2,t,dt,0);
    r_d_cal = d_cal_new(:) - d_cal(:);
    alpha = ((r_d_cal'*r(:))/(r_d_cal'*r_d_cal))*alpha_t;
    toc
    alp = [alp alpha];
    
    v_olds = v;
    v = v + alpha*p;
    
    if mod (iter,25) == 0
        v_new = v;
        v_total = [v_old v_new];
        v_old = v_total;
    end
    
    %% visualization
    figure(2);%'position',[2 340 1320 600])
    subplot(2,2,1);
    imagesc((1:nx)*h,(1:nz)*h,alpha*p);
    title('Update Direction');
    xlabel( ' horizontal distance (m) ');
    ylabel( ' depth (m)');
    colorbar;
    colormap;
    axis image;
    pause(1);
    
    
    subplot(2,2,2);
    imagesc((1:nx)*h,(1:nz)*h,v_initial);
    title('initial model');
    xlabel( ' horizontal distance (m)');
    ylabel( ' depth (m)');
    colorbar;
    colormap;
    axis image;
    caxis([min(min(v_true)) max(max(v_true))]);
    pause(1);
    
    subplot(2,2,3);
    imagesc((1:nx)*h,(1:nz)*h,v_true)
    title('true model')
    xlabel( ' horizontal distance (m)');
    ylabel( ' depth (m)');
    colorbar;
    colormap;
    axis image;
    caxis([min(min(v_true)) max(max(v_true))]);
    pause(1);
    
    subplot(2,3,4);
    imagesc((1:nx)*h,(1:nz)*h,v)
    title('updated model')
    xlabel( ' horizontal distance (m)');
    ylabel( ' depth (m)');
    colorbar;
    colormap;
    axis image;
    caxis([min(min(v_true)) max(max(v_true))]);
    pause(1);
    
    
    
    clear rsw sw d_cal
    
    iter = iter+1
    save('CG-marmousi-nomulti-thesis_final2','e','alp','Model_rms','SNR','fd','fdlp','d_obs','d_obsn', 'r','r_initial','kr','ks','Source','Sourcelp','v','v_initial','h','v_true')
end
save('CG-marmousi-nomulti-thesis_final2','e','alp','v_total','Model_rms','SNR','fd','fdlp','d_obs','d_obsn', 'r','r_initial','kr','ks','Source','Sourcelp','v','v_initial','h','v_true')
toc
delete(gcp('nocreate'))