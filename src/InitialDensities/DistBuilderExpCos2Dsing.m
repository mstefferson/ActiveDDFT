% DistBuilderExpCos2Dsing
% Michael Stefferson
%
% Builds the normalized 3D distributions of the form
% f(\theta) = exp( \sum a_2n P_2n (cos(\theta) ) / Z
%
% Inputs:
% Nc: Max number of coeff
% x : spatial coordin. vector. x = cos(\theta)
% CoeffMat: Matrix of coefficients
%
% Outputs:
% f_mat =  [ f(a_2) ; f(a_2, a_4); ...; f(a_2, a_4,...a_{2Nc})

function f = DistBuilderExpCos2Dsing(Nc,phi,Coeff)
% Distribution

CosSum = zeros(1, length(phi) );
for i = 1:Nc      
          CosSum     = CosSum + Coeff(i) * cos(2*i*phi);    
end
%        Normalize it again to be safe.
 f = exp( CosSum ) ./ (trapz_periodic( phi, exp( CosSum ) ) );
% keyboard