function fluxCmpr( rho, rho_FT, flags, systemObj, particleObj, gridObj )

% Build 3d trig functions
[~,~,phi3D] = meshgrid( gridObj.y, gridObj.x, gridObj.phi );
cosPhi = cos(phi3D);
sinPhi = sin(phi3D);
% Build diffusion matrix
cos2 = cosPhi .^ 2;
sin2 = sinPhi .^ 2;
cossin = cosPhi .* sinPhi;
% Grid spacing for derivatives
dx = gridObj.x(2) - gridObj.x(1);
dy = gridObj.y(2) - gridObj.y(1);
dphi = gridObj.phi(2) - gridObj.phi(1);

% Build diffusion object
if flags.AnisoDiff
    [diffObj] =  DiffMobCoupCoeffCalc( systemObj.Tmp,...
      particleObj.mobPar,particleObj.mobPerp,particleObj.mobRot,...
      gridObj.kx, gridObj.ky, gridObj.km, ...
      gridObj.kx2D, gridObj.ky2D,particleObj.vD);
    Dpar = diffObj.D_par;
    Dperp = diffObj.D_perp;
    D = [ Dpar.*cos2 + Dperp.*sin2, ( Dpar -  Dperp) .* cossin, 0;...
      ( Dpar -  Dperp) .* cossin, Dperp.*cos2 + Dpar.*sin2, 0;...
      0, 0, diffObj.D_rot ];
else
  [diffObj] = DiffMobCoupCoeffCalcIsoDiff(...
    systemObj.Tmp,particleObj.mobPos,particleObj.mobRot, ...
    gridObj.kx, gridObj.ky, gridObj.km);
  D = [ diffObj.D_pos, 0, 0; 0, diffObj.D_pos, 0; 0 0  diffObj.D_rot ];
end

keyboard

% Flux from diffusion
[JxDiff, JyDiff, JphiDiff] = fluxDiff( rho, D, dx, dy, dphi, systemObj );
JxDiffAve = trapz_period( gridObj.phi, JxDiff, 3);
JyDiffAve = trapz_period( gridObj.phi, JyDiff, 3);
JphiDiffAve = trapz_period( gridObj.phi, JphiDiff, 3);

% Flux from interactions
[JxInt, JyInt, JphiInt] = fluxInt( rho, rho_FT, D, diffObj, systemObj );
JxIntAve = trapz_period( gridObj.phi, JxInt, 3);
JyIntAve = trapz_period( gridObj.phi, JyInt, 3);
JphiIntAve = trapz_period( gridObj.phi, JphiInt, 3);

% Flux from driving
[JxDr, JyDr, JphiDr] = fluxDr( rho, particleObj.vD, systemObj, cosPhi, sinPhi );
JxDrAve = trapz_period( gridObj.phi, JxDr, 3);
JyDrAve = trapz_period( gridObj.phi, JyDr, 3);
JphiDrAve = trapz_period( gridObj.phi, JphiDr, 3);

keyboard
