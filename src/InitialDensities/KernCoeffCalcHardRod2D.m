% KernCoeffCalcHardRod2D.m
% Michael Stefferson
%
% Calculate the coefficients of the Legendre expansion of the kernal
% |sin \gamma| for a 2d hard rod interactions.

function d2nVec = KernCoeffCalcHardRod2D(Nc)
% Build a vector of the kernal's coefficients.
d2nVec        = zeros(1,Nc);

for n = 1:Nc;
    d2nVec(n) = 4 / ( pi * ( 4*n^2 - 1) );
end

% keyboard
end