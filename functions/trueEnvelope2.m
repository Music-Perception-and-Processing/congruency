function [SE,CC] = trueEnvelope2(X,C,nCC,theta)

A = 20.*log10(X);
ii = 2;

while sum(A > C(:,end)+theta) > 0

    A = max(20.*log10(X),C(:,ii-1));
    c = dct(A);
    c(nCC+1:end) = 0;
    C(:,ii) = idct(c);

    ii = ii + 1;

end

SE = C(:,end);
CC = c;

end