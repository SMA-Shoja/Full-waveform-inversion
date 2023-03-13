clc
clear
close all
tic
%% True Model
% v_true = zeros(100,100);
% v_true(:,:) = 2500;
% v_true(33:end,1:33) = 3500;
% v_true(50:end,34:70) = 4000;
load('C_true.mat');

v_true = v_true + 2000;
[nz,nx] = size(v_true);
h = 10; % dx = dz = h

%% time sampling
t = 1;
dt = .001;
nt = floor(t/dt);
k = 1:nt;
T = k.*dt;
f = (1/dt)*(0:nt/2)/nt;
%% sources locations
ns_up = 3:5:nx-3;
ms_up = 2*ones(1,length(ns_up));
ms_right = 3:5:nz-3;
ns_right = (nx-2)*ones(1,length(ms_right));
ns_down = 10:10:nx-3;
ms_down = nz-5*ones(1,length(ns_down));
ms_left = 3:5:nz-3;
ns_left = 2*ones(1,length(ms_left));
ms = [ms_up];
ns = [ns_up];
ks = ((ns-1)*nz)+ms;

%% recievers locations
G = 3;
nr_1 = G:5:nx-3;
mr_1 = 2*ones(1,length(nr_1));
mr_2 = G:1:nz-3;
nr_2 = (nx-2)*ones(1,length(mr_2));
nr_3 = nx-3:-1:G;
mr_3 = nz-3*ones(1,length(nr_3));
mr_4 = G:1:nz-3;
nr_4 = 2*ones(1,length(mr_4));
mr = [mr_1];
nr = [nr_1];
kr = ((nr-1)*nz)+mr;


%% forward modeling
fd = 30;
fd_old = fd;
Source = 100*( 1-2.*pi.^2.*fd.^2.*(k.*dt-(0.01)).^2 ).*exp( -pi.^2.*fd.^2.*(k.*dt-(0.01)).^2 );
Source= Source';
[d_obs,sw] = data2Dpartest(nx,nz,h,v_true,ks,kr,Source,t,dt,0);

    
    
 
 rsw = BackPropagationPartest(nx,nz,h,v_true,ks,kr,d_obs,t,dt,0); % residuals backpropagation
    
    %% gradiant of misfit function
    disp(' update direction calculation ')
    tic
    grad_pershot_vec = rsw.*sw;
    grad_pershot_vec = sum(grad_pershot_vec);
    grad = reshape(grad_pershot_vec,[nx*nz,length(ks)]);
    grad_vec = sum(grad,2);
    grad = reshape(grad_vec,[nz,nx]);
   
    image = grad;
    %% visualization
    figure(2);
    
    subplot(1,2,1);
    imagesc((1:nx)*h,(1:nz)*h,v_true)
    title('velocity model')
    xlabel( ' horizontal distance (m)');
    ylabel( ' depth (m)');
    colorbar;
    colormap;
    axis image;
    caxis([min(min(v_true)) max(max(v_true))]);
    pause(1);
    
    subplot(1,2,2);
    imagesc((1:nx)*h,(1:nz)*h,image(10:end,:))
    title('RTM image')
    xlabel( ' horizontal distance (m)');
    ylabel( ' depth (m)');
    colorbar;
    colormap;
    axis image;
    pause(1);
toc