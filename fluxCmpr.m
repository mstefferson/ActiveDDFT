function fluxCmpr( rho, rho_FT, C, flags, systemObj, particleObj, gridObj )

% Build 3d trig functions
[~,~,phi3D] = meshgrid( gridObj.x2, gridObj.x1, gridObj.x3 );
cosPhi = cos(phi3D);
sinPhi = sin(phi3D);
% Build diffusion matrix
cos2 = cosPhi .^ 2;
sin2 = sinPhi .^ 2;
cossin = cosPhi .* sinPhi;
% Grid spacing for derivatives
dx = gridObj.x1(2) - gridObj.x1(1);
dy = gridObj.x2(2) - gridObj.x2(1);
dphi = gridObj.x3(2) - gridObj.x3(1);

% 2D k vectors
[ky2D, kx2D] = meshgrid( gridObj.k2, gridObj.k1 );

% Build diffusion object
if flags.AnisoDiff
    [diffObj] =  DiffMobCoupCoeffCalc( systemObj.Tmp,...
      particleObj.mobPar,particleObj.mobPerp,particleObj.mobRot,...
      gridObj.k1, gridObj.k2, gridObj.km, ...
      kx2D, ky2D,particleObj.fD);
    Dpar = diffObj.D_par;
    Dperp = diffObj.D_perp;
    D.xx = Dpar.*cos2 + Dperp.*sin2;
    D.xy = ( Dpar -  Dperp) .* cossin;
    D.yy = Dperp.*cos2 + Dpar.*sin2;
    D.mm = diffObj.D_rot;
else
  [diffObj] = DiffMobCoupCoeffCalcIsoDiff(...
    systemObj.Tmp,particleObj.mobPos,particleObj.mobRot, ...
    gridObj.k1, gridObj.k2, gridObj.km);
  D.xx = diffObj.D_pos; 
  D.xy = 0;  
  D.yy = diffObj.D_pos;
  D.mm= diffObj.D_rot ;
end
% Mobility for j_int
mob.xx = D.xx ./ systemObj.Tmp;
mob.xy = D.xy ./ systemObj.Tmp;
mob.yy = D.yy ./ systemObj.Tmp;
mob.mm = D.mm ./ systemObj.Tmp;
% Flux from diffusion
%[jxDiff, jyDiff, jphiDiff] = fluxDiff( rho, D, dx, dy, dphi, systemObj );
[jxDiff, jyDiff, jphiDiff] = fluxDiffFt( rho_FT, D, diffObj);
jxDiffAve = trapz_periodic( gridObj.x3, jxDiff, 3);
jyDiffAve = trapz_periodic( gridObj.x3, jyDiff, 3);
jphiDiffAve = trapz_periodic( gridObj.x3, jphiDiff, 3);
jposMagDiff = jxDiffAve .^ 2 + jyDiffAve .^2;
% jmagDiff = jposMagDiff + jphiDiffAve.^2;
jposMagDiff = sqrt( jposMagDiff );
% jmagDiff = sqrt( jmagDiff );

% Flux from interactions
[jxInt, jyInt, jphiInt] = ...
  fluxInt( rho, rho_FT, mob, diffObj, systemObj, particleObj );
jxIntAve = trapz_periodic( gridObj.x3, jxInt, 3);
jyIntAve = trapz_periodic( gridObj.x3, jyInt, 3);
jphiIntAve = trapz_periodic( gridObj.x3, jphiInt, 3);
jposMagInt = jxIntAve .^ 2 + jyIntAve .^2;
% jmagInt = jposMagInt + jphiIntAve.^2;
jposMagInt = sqrt( jposMagInt );
% jmagInt = sqrt( jmagInt );

% Flux from driving
[jxDr, jyDr, jphiDr] = ...
  fluxDrive( rho, particleObj.fD, systemObj, cosPhi, sinPhi );
jxDrAve = trapz_periodic( gridObj.x3, jxDr, 3);
jyDrAve = trapz_periodic( gridObj.x3, jyDr, 3);
% jphiDrAve = trapz_periodic( gridObj.x3, jphiDr, 3);
jposMagDr = jxDrAve .^ 2 + jyDrAve .^2;
% jmagDr = jposMagDr + jphiDrAve.^2;
jposMagDr =  sqrt( jposMagDr );
% jmagDr = sqrt( jmagDr );

% Total
jxT = jxDiff + jxInt + jxDr;
jyT = jyDiff + jyInt + jyDr;
% jTsp  = sqrt( jxT .^ 2 + jyT .^ 2 );
jphiT = jphiDiff + jphiInt + jphiDr;
jxTAve = trapz_periodic( gridObj.x3, jxT, 3);
jyTAve = trapz_periodic( gridObj.x3, jyT, 3);
jphiTAve = trapz_periodic( gridObj.x3, jphiT, 3);
jposMagT  = sqrt( jxTAve .^ 2 + jyTAve .^ 2 );

%%
%% Make subset of vector fields to plot
x = gridObj.x1;
y = gridObj.x2;

fluxPlotSpt( x,y, systemObj, C, ...
  jxDiffAve, jyDiffAve, jxIntAve, jyIntAve,...
  jxTAve, jyTAve, jxDrAve, jyDrAve,...
  jposMagT, jposMagDiff, jposMagInt, jposMagDr);

fluxPlotPhi( x,y, jphiDiffAve, jphiIntAve, jphiTAve )
