function[x] = BackPropagationPartest(nx,nz,h,v,ks,kr,source,t,dt,Anim)
tic
%% residuals back propagation loop
nt = t/dt;

[~,M] = size(source);
source = flipud(source);
    parfor shotnumber = 1:length(ks)
        rsw_new = zeros(nt,nz*nx,length(ks));
        p1 = zeros(nz,nx);
        p2 = zeros(nz,nx);
        p3 = zeros(nz,nx);
        for k = 1:nt
            %% source
            if M > 1
                if k <= length(source)
                    p2(kr) = source(k,(shotnumber-1)*length(kr)+1:shotnumber*length(kr));
                end
            elseif M == 1
                if k <= length(source)
                    p2(kr) = source(k);
                end 
            end
            %% main domain- 2nd order temporal and 4th order spatial FD
            
            j = 3:nz-2; i = 3:nx-2;
            
            p3(j,i) = ( (v(j,i)*dt).^2 ).*(...
                -60*p2(j,i)+...
                16*( p2(j,i-1)+p2(j,i+1)+p2(j-1,i)+p2(j+1,i) )-...
                ( p2(j,i-2)+p2(j,i+2)+p2(j+2,i)+p2(j-2,i) ) )/(12*h^2)+...
                2*p2(j,i)-p1(j,i);
            
            j = 2:nz-1;
            %left
            p3(j,2) = ((v(j,2).^2*dt^2)/(h^2)).*...
                (p2(j,3)-4*p2(j,2)+p2(j,1)+p2(3:nz,2)+p2(j-1,2))+...
                2*p2(j,2)-p1(j,2);
            
            % right
            p3(j,nx-1) = ((v(j,nx-1).^2*dt^2)/(h^2)).*...
                (p2(j,nx)-4*p2(j,nx-1)+p2(j,nx-2)+p2(j+1,nx-1)+p2(j-1,nx-1))+...
                2*p2(j,nx-1)-p1(j,nx-1);
            % top
            i = 2:nx-1;
            p3(2,i) = ((v(2,i).^2*dt^2)/(h^2)).*...
                (p2(2,i+1)-4*p2(2,i)+p2(2,i-1)+p2(3,i)+p2(1,i))+...
                2*p2(2,i)-p1(2,i);
            
            % bottom
            p3(nz-1,i) = ((v(nz-1,i).^2*dt^2)/(h^2)).*...
                (p2(nz-1,i+1)-4*p2(nz-1,i)+p2(nz-1,i-1)+p2(nz,i)+p2(nz-2,i))+...
                2*p2(nz-1,i)-p1(nz-1,i);
            
            %% boundry condition ABC#1
            
            %     %top
            %     p3(1,1:nz) = p2(1,1:nz)+(c(1,1:nz)*dt/h).*(p2(2,1:nz)-p2(1,1:nz));
            %     %left
            %     p3(1:nz,1) = p2(1:nz,1)+(c(1:nz,1)*dt/h).*(p2(1:nz,2)-p2(1:nz,1));
            %     %right
            %     p3(1:nz,nx) = p2(1:nz,nx)+(c(1:nz,nx)*dt/h).*(p2(1:nz,nx-1)-p2(1:nz,nx));
            %     %bottom
            %     p3(nz,1:nz) = p2(nz,1:nz)+(c(nz,1:nx)*dt/h).*(p2(nz-1,1:nz)-p2(nz,1:nz));
            
            %% boundry condition ABC#2
            
            i = 2:nx-1; j = 2:nz-1;
            % top
            p3(1,i) = 2*p2(1,i)-p1(1,i)+...
                ( v(1,i)*dt/h ).*( p2(2,i)-p2(1,i)-p1(2,i)+p1(1,i) )+...
                0.5*( (v(1,i)*(dt/h )).^2 ).*( p2(1,i+1)-2*p2(1,i)+p2(1,i-1) );
            % left
            p3(j,1) = 2*p2(j,1)-p1(j,1)+...
                ( v(j,1)*dt/h ).*( p2(j,2)-p2(j,1)-p1(j,2)+p1(j,1) )+...
                0.5*( (v(j,1)*(dt/h )).^2 ).*( p2(j+1,1)-2*p2(j,1)+p2(j-1,1) );
            %right
            p3(j,nx) = 2*p2(j,nx)-p1(j,nx)-...
                ( v(j,nx)*dt/h ).*( p2(j,nx)-p2(j,nx-1)-p1(j,nx)+p1(j,nx-1) )+...
                0.5*( (v(j,nx)*(dt/h )).^2 ).*( p2(j+1,nx)-2*p2(j,nx)+p2(j-1,nx) );
            %bottom
            p3(nz,i) = 2*p2(nz,i)-p1(nz,i)-...
                ( v(nz,i)*dt/h ).*( p2(nz,i)-p2(nz-1,i)-p1(nz,i)+p1(nz-1,i) )+...
                0.5*( (v(nz,i)*(dt/h )).^2 ).*( p2(nz,i+1)-2*p2(nz,i)+p2(nz,i-1) );
            %top-left corner
            p3(1,1) = p2(1,1)+(((v(1,1)*dt)/(sqrt(2)*h))*(p2(2,1)-p2(1,1)))+...
                (((v(1,1)*dt)/(sqrt(2)*h))*(p2(1,2)-p2(1,1)));
            p3(1,2) = p2(1,2)+(((v(1,2)*dt)/(sqrt(2)*h))*(p2(2,2)-p2(1,2)))+...
                (((v(1,2)*dt)/(sqrt(2)*h))*(p2(1,3)-p2(1,2)));
            p3(2,1) = p2(2,1)+(((v(2,1)*dt)/(sqrt(2)*h))*(p2(3,1)-p2(2,1)))+...
                (((v(2,1)*dt)/(sqrt(2)*h))*(p2(2,2)-p2(2,1)));
            %top-right corner
            p3(1,nx) = p2(1,nx)+(((v(1,nx)*dt)/(sqrt(2)*h))*(p2(2,nx)-p2(1,nx)))-...
                (((v(1,nx)*dt)/(sqrt(2)*h))*(p2(1,nx)-p2(1,nx-1)));
            p3(2,nx) = p2(2,nx)+(((v(2,nx)*dt)/(sqrt(2)*h))*(p2(3,nx)-p2(2,nx)))-...
                (((v(2,nx)*dt)/(sqrt(2)*h))*(p2(2,nx)-p2(2,nx-1)));
            p3(1,nx-1) = p2(1,nx-1)+(((v(1,nx-1)*dt)/(sqrt(2)*h))*(p2(2,nx-1)-p2(1,nx-1)))-...
                (((v(1,nx-1)*dt)/(sqrt(2)*h))*(p2(1,nx-1)-p2(1,nx-2)));
            %bottom_left corner
            p3(nz,1) = p2(nz,1)-(((v(nz,1)*dt)/(sqrt(2)*h))*(p2(nz,1)-p2(nz-1,1)))+...
                (((v(nz,1)*dt)/(sqrt(2)*h))*(p2(nz,2)-p2(nz,1)));
            p3(nz,2) = p2(nz,2)-(((v(nz,2)*dt)/(sqrt(2)*h))*(p2(nz,2)-p2(nz-1,2)))+...
                (((v(nz,2)*dt)/(sqrt(2)*h))*(p2(nz,3)-p2(nz,2)));
            p3(nz-1,1)= p2(nz-1,1)-(((v(nz-1,1)*dt)/(sqrt(2)*h))*(p2(nz-1,1)-p2(nz-2,1)))+...
                (((v(nz-1,1)*dt)/(sqrt(2)*h))*(p2(nz-1,2)-p2(nz-1,1)));
            %bottom-right corner
            p3(nz,nx) = p2(nz,nx)-(((v(nz,nx)*dt)/(sqrt(2)*h))*(p2(nz,nx)-p2(nz-1,nx)))-...
                (((v(nz,nx)*dt)/(sqrt(2)*h))*(p2(nz,nx)-p2(nz,nx-1)));
            p3(nz,nx-1) = p2(nz,nx-1)-(((v(nz,nx-1)*dt)/(sqrt(2)*h))*(p2(nz,nx-1)-p2(nz-1,nx-1)))-...
                (((v(nz,nx-1)*dt)/(sqrt(2)*h))*(p2(nz,nx-1)-p2(nz,nx-2)));
            p3(nz-1,nx) = p2(nz-1,nx)-(((v(nz-1,nx)*dt)/(sqrt(2)*h))*(p2(nz-1,nx)-p2(nz-2,nx)))-...
                (((v(nz-1,nx)*dt)/(sqrt(2)*h))*(p2(nz-1,nx)-p2(nz-1,nx-1)));
            %% time marching
            p1 = p2;
            p2 = p3;
            
            %% animation
            if Anim == 1
                if mod(k,50)==0
                    figure (1)
                    imagesc(p3);
                    axis equal tight
                    colormap gray
                    colorbar
                    %caxis([-.1 .1])
                    pause(0.001);
                end
            end
            
            %% data
            rsw_new(k,:,shotnumber) = p3(1:end); % rsw_new: residuals wavefield sampled at each model parameter per shot
            
        end
        
           rsw(:,:,shotnumber) = rsw_new(:,:,shotnumber);    % IMPORTANR NOTE : rsw is reversed in time
        
    end
    rsw = reshape(rsw,nt,nz*nx*length(ks));
    x = flipud(rsw);
    toc