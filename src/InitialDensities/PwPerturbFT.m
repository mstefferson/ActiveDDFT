% PwPerturbFT.m
% Perturbs input rhp
% Perturbations of the form
%  \sum epsilon_k exp( i (k_x x + k_y y + k_phi phi)
% epsilon_k is an input parameter
% Does everything in Fourier space

function [rho] = PwPerturbFT(rho,systemObj,rhoInit)

% N3 can be commonly used. Declare it
N3 = systemObj.n1 * systemObj.n2 * systemObj.n3;
% Perturb coeff is the weight times equilbrium
% concentration. Make sure it isn't too large
% Find isotropic density
isoDen = systemObj.c / ( systemObj.l3);

maxPerturb = isoDen * rhoInit.WeightPert * ...
  (2*rhoInit.NumModesX) * (2*rhoInit.NumModesY) * (2*rhoInit.NumModesM);
if min(min(min(rho))) < maxPerturb
  coeffMax = isoDen .* ...
    rhoInit.WeightPert * min(min(min(rho))) / maxPerturb;
else
  coeffMax = isoDen .* rhoInit.WeightPert;
end
% If it's not random, set Coeff ourside the loop
% scale by N3 b/c of FT factor
if rhoInit.RandomAmp == 0
  coeff = coeffMax * N3 * ( 1 + sqrt(-1) );
else
  coeffTemp = coeffMax * N3;
end
% Handle perturbations in Fourier space
rhoFT = fftshift( fftn( rho ) );

% Give index of k = 0
if systemObj.n1 ~= 1
  kx0   = systemObj.n1 / 2 + 1;
else
  kx0 = 1;
end
if systemObj.n2 ~= 1
  ky0   = systemObj.n2 / 2 + 1;
else
  ky0 = 1;
end
if systemObj.n3 ~= 1
  km0   = systemObj.n3 / 2 + 1;
else
  km0 = 1;
end

% Build Perturbation Matrix

m1 = rhoInit.NumModesX;
m2 = rhoInit.NumModesY;
m3 = rhoInit.NumModesM;
perturb =  zeros( 2*m1 + 1, 2*m2 + 1, 2*m3 + 1 );

if rhoInit.RandomAmp == 1;
  perturb( 1 : m1+1, 1 : m2+1, 1 : m3+1 ) = ...
    coeffTemp .* ( ( -1 + 2 .* rand( m1+1, m2+1, m3+1 ) ) + ...
    + sqrt(-1) .*  ( -1 + 2 .* rand( m1+1, m2+1, m3+1 ) ) );
  perturb( m1+1 : 2*m1+1, m2+1 : 2*m2+1, 1 : m3 ) = ...
    coeffTemp .* ( ( -1 + 2 .* rand( m1+1, m2+1, m3 ) ) + ...
    + sqrt(-1) .* ( -1 + 2 .* rand( m1+1, m2+1, m3 ) ) ) ;
  perturb( 1 : m1, m2+2 : 2*m2+1, 1 : m3+1 ) = ...
    coeffTemp .* ( ( -1 + 2 .* rand( m1, m2, m3+1 ) ) + ...
    sqrt(-1) .* ( -1 + 2 .* rand( m1, m2, m3+1 ) ) );
  perturb( m1+2 : 2*m1+1, 1 : m2 , 1 : m3 ) = ...
    coeffTemp .* ( ( -1 + 2 .* rand( m1, m2, m3 ) ) + ...
    + sqrt(-1) .* ( -1 + 2 .* rand( m1, m2, m3 ) ) );
else
  perturb( 1 : m1+1, 1 : m2+1, 1 : m3+1 ) = coeff;
  perturb( m1+1 : 2*m1+1, m2+1 : 2*m2+1, 1 : m3+1 ) = coeff;
  perturb( 1 : m1, m2+2 : 2*m2+1, 1 : m3+1 ) = coeff;
  perturb( m1+2 : 2*m1+1, 1 : m2 , 1 : m3 ) = coeff;
end

% Reflection
perturb( m1+1 : 2*m1+1, m2+1 : 2*m2+1, m3+1 : 2*m3+1 ) = ...
  conj(flip(flip(flip( perturb( 1 : m1+1, 1 : m2+1, 1 : m3+1  ),1),2),3) );
perturb( 1 : m1+1, 1 : m2+1,  m3+1 : 2*m3+1 ) = ...
  conj(flip(flip(flip( perturb(  m1+1 : 2*m1+1, m2+1 : 2*m2+1, 1 : m3+1 ),1),2),3) ) ;
perturb( m1+2 : 2*m1+1, 1 : m2,  m3+1 : 2*m3+1 ) = ...
  conj(flip(flip(flip( perturb(  1 : m1, m2+2 : 2*m2+1, 1 : m3+1  ),1),2),3) );
perturb( 1 : m1, m2+2 : 2*m2+1,  m3+2 : 2*m3+1 ) = ...
  conj(flip(flip(flip( perturb(  m1+2 : 2*m1+1, 1 : m2 , 1 : m3  ),1),2),3) );

perturb(m1+1,m2+1,m3+1) = 0;


% Add perturbation to rhoFT
rhoFT( kx0 - m1 : kx0 + m1, ky0 - m2 : ky0 + m2, km0 - m3 : km0 + m3) = ...
  rhoFT( kx0 - m1 : kx0 + m1, ky0 - m2 : ky0 + m2, km0 - m3 : km0 + m3) + perturb;
% Inverse transform and Take real part
rho = real( ifftn( ifftshift( rhoFT ) ) );

% keyboard
% Fix negative rho if that happened.
[rho] = FixNegDenFnc(rho);

end
