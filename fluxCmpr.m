function fluxCmpr( rho, rho_FT, C, N, P, flags, systemObj, particleObj, gridObj )

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

% 2D k vectors
[ky2D, kx2D] = meshgrid( gridObj.ky, gridObj.kx );
keyboard
% Build diffusion object
if flags.AnisoDiff
    [diffObj] =  DiffMobCoupCoeffCalc( systemObj.Tmp,...
      particleObj.mobPar,particleObj.mobPerp,particleObj.mobRot,...
      gridObj.kx, gridObj.ky, gridObj.km, ...
      kx2D, ky2D,particleObj.vD);
    Dpar = diffObj.D_par;
    Dperp = diffObj.D_perp;
    D.xx = Dpar.*cos2 + Dperp.*sin2;
    D.xy = ( Dpar -  Dperp) .* cossin;
    D.yy = Dperp.*cos2 + Dpar.*sin2;
    D.mm = diffObj.D_rot;
else
  [diffObj] = DiffMobCoupCoeffCalcIsoDiff(...
    systemObj.Tmp,particleObj.mobPos,particleObj.mobRot, ...
    gridObj.kx, gridObj.ky, gridObj.km);
  D.xx = diffObj.D_pos; 
  D.xy = 0;  
  D.yy = diffObj.D_pos;
  D.mm= diffObj.D_rot ;
end

% Flux from diffusion
[jxDiff, jyDiff, jphiDiff] = fluxDiff( rho, D, dx, dy, dphi, systemObj );
jxDiffAve = trapz_periodic( gridObj.phi, jxDiff, 3);
jyDiffAve = trapz_periodic( gridObj.phi, jyDiff, 3);
jphiDiffAve = trapz_periodic( gridObj.phi, jphiDiff, 3);
jposMagDiff = jxDiffAve .^ 2 + jyDiffAve .^2;
jmagDiff = jposMagDiff + jphiDiffAve.^2;
jposMagDiff = sqrt( jposMagDiff );
% jmagDiff = sqrt( jmagDiff );


% Flux from interactions
[jxInt, jyInt, jphiInt] = ...
  fluxInt( rho, rho_FT, D, diffObj, systemObj, particleObj );
jxIntAve = trapz_periodic( gridObj.phi, jxInt, 3);
jyIntAve = trapz_periodic( gridObj.phi, jyInt, 3);
jphiIntAve = trapz_periodic( gridObj.phi, jphiInt, 3);
jposMagInt = jxIntAve .^ 2 + jyIntAve .^2;
jmagInt = jposMagInt + jphiIntAve.^2;
jposMagInt = sqrt( jposMagInt );
% jmagInt = sqrt( jmagInt );

% Flux from driving
[jxDr, jyDr, jphiDr] = ...
  fluxDrive( rho, particleObj.vD, systemObj, cosPhi, sinPhi );
jxDrAve = trapz_periodic( gridObj.phi, jxDr, 3);
jyDrAve = trapz_periodic( gridObj.phi, jyDr, 3);
jphiDrAve = trapz_periodic( gridObj.phi, jphiDr, 3);
jposMagDr = jxDrAve .^ 2 + jyDrAve .^2;
jmagDr = jposMagDr + jphiDrAve.^2;
jposMagDr =  sqrt( jposMagDr );
% jmagDr = sqrt( jmagDr );

%%
figure()
% Diffusion Quiver
subplot( 2, 3, 1)
imagesc( gridObj.x, gridObj.y, jposMagDiff' );
colorbar
hold on
quiver( gridObj.x, gridObj.y, jxDiffAve', jyDiffAve','color',[1,1,1])
hold off
ax = gca;
ax.YDir = 'normal';
axis square
title('Diffusion spatial flux');
xlabel('x'); ylabel('y');

% Interactions Quiver
subplot( 2, 3, 2)
imagesc( gridObj.x, gridObj.y, jposMagInt' )
colorbar
hold on
quiver( gridObj.x, gridObj.y, jxIntAve', jyIntAve','color',[1,1,1])
hold off
ax = gca;
ax.YDir = 'normal';
axis square
title('Interaction spatial flux');
xlabel('x'); ylabel('y');

% Driving Quiver
subplot( 2, 3, 3)
imagesc( gridObj.x, gridObj.y, jposMagDr' )
colorbar
hold on
quiver( gridObj.x, gridObj.y, jxDrAve', jyDrAve','color',[1,1,1])
hold off
ax = gca;
ax.YDir = 'normal';
axis square
title('Driving spatial flux');
xlabel('x'); ylabel('y');

% Diffusion angular magnitude
subplot( 2, 3, 4)
imagesc( gridObj.x, gridObj.y, jphiDiffAve' )
colorbar
ax = gca;
ax.YDir = 'normal';
axis square
title('Diffusion angular flux');
xlabel('x'); ylabel('y');

% Interaction angular magnitude
subplot( 2, 3, 5)
imagesc( gridObj.x, gridObj.y, jphiIntAve' )
colorbar
ax = gca;
ax.YDir = 'normal';
axis square
title('Interaction angular flux');
xlabel('x'); ylabel('y');

% driving angular magnitude
subplot( 2, 3, 6)
imagesc( gridObj.x, gridObj.y, jphiDrAve' )
colorbar
ax = gca;
ax.YDir = 'normal';
axis square
title('Driving angular flux');
xlabel('x'); ylabel('y');

