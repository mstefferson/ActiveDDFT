% pairDistPlotSingle( pDist, lR, lC, ttlstr)
% Plot a single order parameter and a slice through x = 0 and y = 0
%
function pairDistPlotSingle( pDist, lR, lC, ttlstr)
  [nR, nC] = size(pDist);
  % set some things
  fontSize = 10;
  % grud stuff
  % since rotated, x->n2, y->n1
  centerX = nC/2+1;
  centerY = nR/2+1;
  % centering limits
  xLim = lC * ( 1 / 2 - 1/nC);
  yLim = lR * ( 1 / 2 - 1/nR);
  x = lC / nC * ( -nC/2:nC/2-1);
  y = lR / nR * ( -nR/2:nR/2-1);
  figure()
  maxVal = max(pDist(:)); 
  minVal = min( min(pDist(:)), 0 );
  % full
  ax1 = subplot(1,2,1);
  ax1.FontSize = fontSize;
  imagesc( x, y, pDist );
  axis square
  xlabel('$$x$$'); ylabel('$$y$$')
  title(ttlstr)
  ax1.XLim = [-xLim xLim];
  ax1.YLim = [-yLim yLim];
  ax1.CLim = [minVal maxVal];
  ax1.YDir = 'normal';
  colorbar
  % not, x are columns, y rows,
  ax2 = subplot(1,2,2);
  plot( x(2:end), pDist(centerY,2:end),...
    y(2:end), pDist(2:end,centerX) )
  axis square
  xlabel('position'); ylabel('g(r)');
  title('slice')
  currLim = min(xLim,yLim);
  ax2.XLim = [-currLim currLim];
  leg = legend('$$ x $$','$$ y $$');
  leg.Interpreter = 'latex';
  leg.Position = [ 0.85 0.45 0.0393 0.05 ];
  leg.FontSize = fontSize -2 ;
  ax2.FontSize = fontSize - 2;
end
