function [LP1,LP2] = LowpassKaiser(pb,sb,spa,pbr,pb2,sb2,spa2,pbr2,dt)
LP1 = designfilt('lowpassfir','PassbandFrequency',pb,'StopbandFrequency',sb,...
    'StopbandAttenuation',spa,'PassbandRipple',pbr,'SampleRate', 1/dt,...
    'DesignMethod','kaiserwin');
LP2 = designfilt('lowpassfir','PassbandFrequency',pb2,'StopbandFrequency',sb2,...
    'StopbandAttenuation',spa2,'PassbandRipple',pbr2,'SampleRate', 1/dt,...
    'DesignMethod','kaiserwin');
% fvtool(lpkai);
% %
% %   lpham = designfilt('lowpassfir','CutoffFrequency', 7,...
% % %       'SampleRate', 1000,'FilterOrder',200,'DesignMethod','Hamming');
% %
% t = 1;
% dt = .001;
% nt = t/dt;
% k = 1:nt;
% fd =30;
% f = (1/dt)*(0:nt/2)/nt;
% Source = ( 1-2.*pi.^2.*fd.^2.*(k.*dt-(0.2)).^2 ).*exp( -pi.^2.*fd.^2.*(k.*dt-(0.2)).^2 );
% 
% %
% % df = filter(Hr,Source);
% % [n,m] = size(d_obs);
% % d_obs2 = [d_obs;zeros(3*length(lpkai.Coefficients),m)];
% %  Source2 = [Source zeros(1,3*length(lpkai.Coefficients))];
%             Source2 = filtfilt(Hd,Source);
% %             Source2 = Source2(1:length(Source));
% % df = df(1:n,:);
% figure(1)
% plot(k,Source,k,Source2)
% ssp = fft(Source2);
% Ssp = fft(Source);
% figure(2)
% plot(f,abs(Ssp(1:length(Ssp)/2+1)),f,abs(ssp(1:length(ssp)/2+1)))
% figure(3)
% subplot(1,2,1)
% imagesc(d_obs2)
% subplot(1,2,2)
% imagesc(df)