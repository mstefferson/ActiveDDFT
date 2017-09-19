% f = freeEnergyCalcId( rho, systemObj )
%
% Calculates the ideal gas contribution to the free energy
%
function f = freeEnergyCalcId( rho, systemObj )
% grid stuff
dx1 = systemObj.l1/systemObj.n1;
x1 = dx1 .* (0:systemObj.n1-1);
dx2 = systemObj.l2/systemObj.n2;
x2 = dx2 .* (0:systemObj.n2-1);
% Free energy for a 2D non-interacting gas
integrand = rho .* ( log( rho ) - 1 );
if systemObj.n3 == 1
  f = trapz_periodic( x1, trapz_periodic( x2, integrand, 2 ), 1 );
else
  dx3 = systemObj.l3/systemObj.n3;
  x3 = dx3 .* (0:systemObj.n3-1);
  f = trapz_periodic( x1, trapz_periodic( x2, trapz_periodic( ...
    x3, integrand, 3 ), 2 ), 1 );
end
% scale by kbT
f = systemObj.tmp * f;
end
