function [JxDiff, JyDiff, JphiDiff] = ...
  fluxDiff( rho, D, dx, dy, dphi, systemObj )
% Commonly used variables
n1 = systemObj.n1;
n2 = systemObj.n2;
n3 = systemObj.n3;
% Allocate
drho_dx = zeros( n1, n2, n3 );
drho_dy = zeros( n1, n2, n3 );
drho_dphi =  zeros( n1, n2, n3 );

% Take derivatives on points. df(i) = f(i+1) - f(i-1)
% Bulk derivatives
drho_dx( 2:n1-1, :, : ) = 1 / ( 2 * dx ) .* ...
  ( rho( 3:n1, :, : ) - rho( 1:n1-2, :, : ) );
drho_dy( :, 2:n2-1, : ) = 1 / ( 2 * dy ) .* ...
  ( rho( :, 3:n2, : ) - rho( :, 1:n2-2, : ) );
drho_dphi( :, :, 2:n3-1 ) = 1 / ( 2 * dphi ) .* ...
  ( rho( :, :, 3:n3 ) - rho( :, :, 1:n3-2 ) );
% Surface derivatives
drho_dx( 1, :, : ) = 1 / ( 2 * dx ) .* ...
  ( rho( 2, :, : ) - rho( n1, :, : ) );
drho_dy( :, 1, : ) = 1 / ( 2 * dy ) .* ...
  ( rho( :, 2, : ) - rho( :, n2, : ) );
drho_dphi( :, : , 1) = 1 / ( 2 * dphi ) .* ...
  ( rho( :, : , 2 ) - rho( :, : , n3 ) );
drho_dx( n1, :, : ) = 1 / ( 2 * dx ) .* ...
  ( rho( 1, :, : ) - rho( n1-1, :, : ) );
drho_dy( :, n2, : ) = 1 / ( 2 * dy ) .* ...
  ( rho( :, 1, : ) - rho( :, n2-1, : ) );
drho_dphi( :, : , n3 ) = 1 / ( 2 * dphi ) .* ...
  ( rho( :, : , 1 ) - rho( :, : , n3-1 ) );

% Calculate fluxes
JxDiff = - ( D.xx .* drho_dx + D.xy .* drho_dy );
JyDiff = - ( D.xy .* drho_dx + D.yy .* drho_dy );
JphiDiff = - ( D.mm .* drho_dphi );

