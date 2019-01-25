function [PDC,A] = PDC(X,i,j,f)

Xf = fft(X);
[r,c] = find(min(abs(Xf-f)));
for m = 1:size(X,1)
    for n = 1:size(X,2)
        Af(m,n) = corr(X(m,n),r(j));
    end
end

for j = 1:size(X,2)
    for i = 1:j
    PDCX{i,j} = sum(Af(:,i:j),2)./sqrt((abs(Af(:,j))).^2);
    end
end

area(PDCX{i,j});

end