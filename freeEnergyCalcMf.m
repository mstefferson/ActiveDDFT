% f = freeEnergyCalcMf( rho, systemObj )
%
% Calculates the mean field (or 2nd virial) contribution to the free energy
%
function f = freeEnergyCalcMf( rho, systemObj, particleObj )
% grid stuff
dx1 = systemObj.l1/systemObj.n1;
x1 = dx1 .* (0:systemObj.n1-1);
dx2 = systemObj.l2/systemObj.n2;
x2 = dx2 .* (0:systemObj.n2-1);
% fft rho
rhoFt = fftshift( fftn( rho ) );
% initialize mu
muEx = zeros( systemObj.n1, systemObj.n2, systemObj.n3 );
% get interobj
gridObj.k3ind0 = floor( systemObj.n3 / 2 ) + 1;
[interObj] = interObjMaker( particleObj, systemObj, gridObj );
% get chemical potential
% check for softshoulder
if strcmp( particleObj.interLr, 'softshoulder')
  % chemical potential
  [muExFt] =  muExCalcMfFt( rhoFt(:,:,interObj.k3ind0), interObj.vFt, systemObj, interObj.muMfScale );
  % get it in real space
  muEx = muEx + real( ifftn( ifftshift( muExFt ) ) );
end
% check for mayer
if strcmp( particleObj.interHb, 'mayer')
  % chemical potential mayer
  muExFt = muExCalcMayerLF(rhoFt, interObj.FmFt, systemObj,...
      interObj.muMayerScale, interObj.muMayerInds, interObj.muMayerMinusInds);
  % get it in real space
  muEx = muEx + real( ifftn( ifftshift( muExFt ) ) );
end
% free energy for mean field (or onsager)
integrand = muEx .* rho;
if systemObj.n3 == 1
  f = trapz_periodic( x1, trapz_periodic( x2, integrand, 2 ), 1 );
else
  dx3 = systemObj.l3/systemObj.n3;
  x3 = dx3 .* (0:systemObj.n3-1);
  f = trapz_periodic( x1, trapz_periodic( x2, trapz_periodic( ...
    x3, integrand, 3 ), 2 ), 1 );
end
end
