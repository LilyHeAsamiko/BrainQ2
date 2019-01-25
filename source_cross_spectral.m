% source cross spectrum computation
function [E, SNR, SCS] = source_cross_spectral(L,X,f,e)
%X = zef.measurements(e,:);
%L = zef_MEG_L;
L = L(:,1:(length(L)/9));
[CS, E, SNR] = cross_spectral(X,f,e);
SCS = L'*CS*L;
showm = SCS(:);
ind = 1:size(L,2):size(L,2)^2;
figure,
spectrogram(showm(ind));
title('Source cross sprctral')
end