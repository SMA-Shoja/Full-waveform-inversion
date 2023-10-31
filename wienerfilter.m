function output= wienerfilter(originialwavelet,targetwavelet,datatobefiltered)

[~,n] = size(datatobefiltered);
wo = fft(originialwavelet);
wt = fft(targetwavelet);

winerlp = (wt.*conj(wo))./(abs((wo.*wo))+1e-5);

wf = repmat(winerlp(:),1,n).*fft(datatobefiltered);

output = real(ifft(wf));
