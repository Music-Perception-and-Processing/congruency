function trueSE = cepstralSmoothing(S,C,fb,cf)

A = zeros(10,length(S));
A(1,:) = log10(S);
theta = 2;

for i = 2:10

    A(i,:) = max([20.*log10(S); C(i-1,:)]);

    erbA(i,:) = fb*10.^(A(i,:)./20)';

    erbC = 20.*dct(log10(erbA(i,:)),[],2);
    erbC = idct(erbC,[],2);
    C(i,:) = interp1(cf,erbC,1:22050,'spline');

end

trueSE = C(end,:);

end

