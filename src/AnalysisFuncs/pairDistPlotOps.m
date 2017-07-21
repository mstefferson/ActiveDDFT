% plot the three order parameter pair distribution functions
% g0, g1, g2. Make sure to pass the rotated and centered or it's
% going to look strange.
%
function pairDistPlotOps( pDist0, pDist1, pDist2, lR, lC)
[nR, nC] = size(pDist0);
% set some things
fontSize = 10;
% grid stuff
xLim = lC * ( 1 / 2 - 1/nC);
yLim = lR * ( 1 / 2 - 1/nR);
x = lC / nC * ( -nC/2:nC/2-1);
y = lR / nR * ( -nR/2:nR/2-1);
% plot the three pair distributions
figure()
% g0
[minVal, maxVal] = calcMinMaxVal( pDist0);
ax0 = subplot(1,3,1);
ax0.FontSize = fontSize;
imagesc( x, y, pDist0 );
title('$$g_0(x,y)$$')
axis square
xlabel('$$x$$'); ylabel('$$y$$')
ax0.XLim = [-xLim xLim];
ax0.YLim = [-yLim yLim];
ax0.CLim = [minVal maxVal];
ax0.YDir = 'normal';
colorbar;
% g1 (polar)
ax1 = subplot(1,3,2);
[minVal, maxVal] = calcMinMaxVal( pDist1 );
ax1.FontSize = fontSize;
imagesc( x, y, pDist1 );
title('$$ g_1(x,y) $$')
axis square
xlabel('$$x$$'); ylabel('$$y$$')
ax1.XLim = [-xLim xLim];
ax1.YLim = [-yLim yLim];
ax1.CLim = [minVal maxVal];
ax1.YDir = 'normal';
colorbar;
% g2 (nem)
[minVal, maxVal] = calcMinMaxVal( pDist2 );
ax2 = subplot(1,3,3);
ax2.FontSize = fontSize;
imagesc( x, y, pDist2 );
title('$$ g_2 (x,y) $$')
axis square
xlabel('$$x$$'); ylabel('$$y$$')
ax2.XLim = [-xLim xLim];
ax2.YLim = [-yLim yLim];
ax2.CLim = [minVal maxVal];
ax2.YDir = 'normal';
colorbar;
end

function [minVal, maxVal] = calcMinMaxVal( tempDist )
maxVal = max(tempDist (:));
minVal = min( min(tempDist (:)), 0 );
if abs( minVal-maxVal ) < 0.1
  maxVal = minVal + 0.1;
end
end