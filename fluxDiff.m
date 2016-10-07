function [JxDiff, JyDiff, JphiDiff]= fluxDiff( rho, D, dx, dy, dphi, systemObj )
% Commonly used variables
Nx = systemObj.Nx;
Ny = systemObj.Ny;
Nm = systemObj.Nm;
% Allocate
drho_dx = zeros( Nx, Ny, Nm );
drho_dy = zeros( Nx, Ny, Nm );
drho_dphi =  zeros( Nx, Ny, Nm );

% Take derivatives on points. df(i) = f(i+1) - f(i-1)
% Bulk derivatives
drho_dx( 2:Nx-1, :, : ) = 1 / ( 2 * dx ) .* ...
  ( rho( 3:Nx, :, : ) - rho( 1:Nx-2, :, : ) );
drho_dy( :, 2:Ny-1, : ) = 1 / ( 2 * dy ) .* ...
  ( rho( :, 3:Ny, : ) - rho( :, 1:Ny-2, : ) );
drho_dphi( :, :, 2:Nm-1 ) = 1 / ( 2 * dphi ) .* ...
  ( rho( :, :, 3:Nm ) - rho( :, :, 1:Nm-2 ) );
% Surface derivatives
drho_dx( 1, :, : ) = 1 / ( 2 * dx ) .* ...
  ( rho( 2, :, : ) - rho( Nx, :, : ) );
drho_dx( Nx, :, : ) = 1 / ( 2 * dx ) .* ...
  ( rho( Nx-1, :, : ) - rho( 1, :, : ) );
drho_dy( :, 1, : ) = 1 / ( 2 * dy ) .* ...
  ( rho( :, 2, : ) - rho( :, Ny, : ) );
drho_dy( :, Ny, : ) = 1 / ( 2 * dy ) .* ...
  ( rho( :, Ny-1, : ) - rho( :, 1, : ) );
drho_dphi( :, : , 1) = 1 / ( 2 * dphi ) .* ...
  ( rho( :, : , 2 ) - rho( :, : , Nm ) );
drho_dphi( :, : , Nm ) = 1 / ( 2 * dphi ) .* ...
  ( rho( :, : , Nm-1 ) - rho( :, : , 1 ) );

% Calculate fluxes
JxDiff = - ( D(1,1) .* drho_dx + D(1,2) .* drho_dy );
JyDiff = - ( D(2,1) .* drho_dx + D(2,2) .* drho_dy );
JphiDiff = - ( D(3,3) .* drho_dphi );

