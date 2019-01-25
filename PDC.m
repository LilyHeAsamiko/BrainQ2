function [j, PDCX] = PDC(X, f)
    A = corr(X);
    Af = fft(A);
    j = min(abs(Af-f),1);
    a = Af(:,j(1,2));
    denom = sqrt(abs(a).^2); 
    PDCX = Af./denom;
    figure,
    spectrogram(PDCX(:))
    title([{'the frequency '},{int2str(f)},{' Hz is about at time '},{int2str(j(1,1))},{' ms with PDC approximation'}])
end