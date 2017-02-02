% Plot a quiver on top of C with the magnitude below it

function fluxPlotSpt( x,y, systemObj, C, ...
  jxDiffAve, jyDiffAve, jxIntAve, jyIntAve,...
  jxTAve, jyTAve, jxDrAve, jyDrAve,...
  jposMagT, jposMagDiff, jposMagInt, jposMagDr)

divNumX = 8;
divNumY = 8;
deltaX  = ceil(systemObj.n1 / divNumX );
deltaY  = ceil(systemObj.n2 / divNumY);
subIndX = 1:deltaX:(systemObj.n1 + 1 - deltaX);
subIndY = 1:deltaY:(systemObj.n2 + 1 - deltaY);
xSub = x(subIndX);
ySub = y(subIndY);

figure()

% Diff Sp Flux Quiver
subplot(2,4,1);
imagesc( x, y, C' );
colorbar
hold on
quiver( xSub, ySub, jxDiffAve(subIndX,subIndY)', ...
  jyDiffAve(subIndX,subIndY)', 'color',[1,1,1])
hold off
ax = gca;
ax.YDir = 'normal';
axis square
title('diff spatial flux');
xlabel('x'); ylabel('y');

% Int Sp Flux Quiver
subplot(2,4,2);
imagesc( x, y, C' );
colorbar
hold on
quiver( xSub, ySub, jxIntAve(subIndX,subIndY)', ...
  jyIntAve(subIndX,subIndY)', 'color',[1,1,1])
hold off
ax = gca;
ax.YDir = 'normal';
axis square
title('int spatial flux');
xlabel('x'); ylabel('y');

% Int Sp Flux Quiver
subplot(2,4,3);
imagesc( x, y, C' );
colorbar
hold on
quiver( xSub, ySub, jxDrAve(subIndX,subIndY)', ...
  jyDrAve(subIndX,subIndY)', 'color',[1,1,1])
hold off
ax = gca;
ax.YDir = 'normal';
axis square
title('drive spatial flux');
xlabel('x'); ylabel('y');

% Total Sp Flux Quiver
subplot(2,4,4);
imagesc( x, y, C' );
colorbar
hold on
quiver( xSub, ySub, jxTAve(subIndX,subIndY)', ...
  jyTAve(subIndX,subIndY)', 'color',[1,1,1])
hold off
ax = gca;
ax.YDir = 'normal';
axis square
title('total spatial flux');
xlabel('x'); ylabel('y');

% Total Sp Flux Quiver
subplot(2,4,5);
imagesc( x, y, jposMagDiff' );
colorbar
ax = gca;
ax.YDir = 'normal';
axis square
title('diff spatial flux');
xlabel('x'); ylabel('y');

% Total Sp Flux Quiver
subplot(2,4,6);
imagesc( x, y, jposMagInt' );
colorbar
ax = gca;
ax.YDir = 'normal';
axis square
title('int spatial flux');
xlabel('x'); ylabel('y');

% Total Sp Flux Quiver
subplot(2,4,7);
imagesc( x, y, jposMagDr' );
colorbar
ax = gca;
ax.YDir = 'normal';
axis square
title('Drive spatial flux');
xlabel('x'); ylabel('y');

% Total Sp Flux Quiver
subplot(2,4,8);
imagesc( x, y, jposMagT' );
colorbar
ax = gca;
ax.YDir = 'normal';
axis square
title('total spatial flux');
xlabel('x'); ylabel('y');


