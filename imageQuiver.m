function imageQuiver( amps, xDir, yDir, method, ax)
% Set-up figure
if nargin < 5
  figure()
  ax = gca;
end
ax.NextPlot = 'add';
% Set-up ticks
% maxTicks = 6;
% [nr, nc] = size( amps );
% xTicks = 1:min(maxTicks,nc);
% yTicks = 1:min(maxTicks,nr);
% My method M( x,y )
if method == 1
  ax.YDir = 'rev';
  imagesc( ax, amps );
  % Flip y direction so y axis counts ( 1, 2, 3,...) to math p quiver
  p = quiver(ax, xDir, -yDir,0.5);
  p.Color = [1 1 1];
  xlabel('columns'); ylabel('rows');
  title('Method 1: exact Copy')
%   ax.XTick      = xTicks;
%   ax.YTick      = flip(yTicks);
%   ax.XTickLabel = xTicks;
%   ax.YTickLabel = yTicks; 
% MATLAB's method M( y,x )
elseif method == 2
  ax.YDir = 'normal';
  imagesc( ax, amps );
  % Flip y direction so y axis counts ( 1, 2, 3,...) to math p quiver
  p = quiver( ax, xDir, yDir, 0.5);
  p.Color = [1 1 1];
  xlabel('columns'); ylabel('rows');
  title('Method 2: flipped rows')
%   ax.XTick      = xTicks;
%   ax.YTick      = yTicks;
%   ax.XTickLabel = xTicks;
%   ax.YTickLabel = yTicks; 
else
  error('Not written')
end
ax.NextPlot = 'replace';
