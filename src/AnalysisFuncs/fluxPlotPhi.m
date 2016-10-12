% Plot a quiver on top of C with the magnitude below it

function fluxPlotPhi( x,y, jphiDiffAve, jphiIntAve, jphiTAve )

%% Quiver plots
figure()

% Diffusion angular magnitude
subplot( 1, 3, 1)
imagesc( x, y, jphiDiffAve' )
colorbar
ax = gca;
ax.YDir = 'normal';
axis square
title('Diffusion angular flux');
xlabel('x'); ylabel('y');

% Interaction angular magnitude
subplot( 1, 3, 2)
imagesc( x, y, jphiIntAve' )
colorbar
ax = gca;
ax.YDir = 'normal';
axis square
title('Interaction angular flux');
xlabel('x'); ylabel('y');

% Total angular magnitude
subplot( 1, 3, 3)
imagesc( x, y, jphiTAve' )
colorbar
ax = gca;
ax.YDir = 'normal';
axis square
title('Total angular flux');
xlabel('x'); ylabel('y');

