n1 = 100;
n2 = 100;
n3 = 100;

l1 = 10;
l2 = 10;
l3 = 2*pi;
%
xcg1 = l1/4;
xcg2 = 0;
xcg3 = 0;
sigG1 = l1 / 10;
sigG2 = l1 / 20;
sigG3 = l3 / 20;
xcgIndShift1 = mod( round(n1 * xcg1 / l1), n1);
xcgIndShift2 = mod( round(n2 * xcg2 / l2), n2);
xcgIndShift3 = mod( round(n3 * xcg3 / l3), n3);
peakInd1 = xcgIndShift1 + n1/2 + 1;
peakInd2 = xcgIndShift2 + n2/2 + 1;
peakInd3 = xcgIndShift3 + n3/2 + 1;
%
x1 = l1 / n1 * (-n1/2:n1/2-1);
x2 = l2 / n2 * (-n2/2:n2/2-1);
x3 = l3 / n3 * (-n3/2:n3/2-1);

[x2mesh2d, x1mesh2d] = meshgrid( x2, x1 );
[x2mesh3d, x1mesh3d, x3mesh3d] = meshgrid( x2, x1, x3);

cCenter = exp( - 1 / 2 * ( ( x1mesh2d ./ sigG1 ) .^ 2 + ...
  ( x2mesh2d ./ sigG2 ) .^ 2 ) );
c = circshift( circshift( cCenter, xcgIndShift1, 1 ), xcgIndShift2, 2 );

% build full
cosphi = cos(x3);
cosphi = repmat( reshape( cosphi, [1 1 n3] ), [n1 n2] );
sinphi = sin(x3);
sinphi = repmat( reshape( sinphi, [1 1 n3] ), [n1 n2] );

% iso
if 0
  f1 = 1 / ( l3) .* ones( 1, 1, n3);
  rho1 = c .* f1;
  [c1, p1, px1, py1] = calcCP( rho1, x3, cosphi, sinphi );
  plotCP( c1, p1, px1, py1,x1,x2 )
end
% delta
if 0
  phiDelta = -pi/2;
  phiInd = mod( round( n3 * phiDelta / (l3) ) + n3/2, n3) + 1;
  f2 = zeros( 1, 1, n3);
  f2(phiInd) = phiInd;
  rho2 = c .* f2;
  [c2, p2, px2, py2] = calcCP( rho2, x3, cosphi, sinphi );
  plotCP( c2, p2, px2, py2, x1, x2 )
end
% gaussian?
if 0
  cCenter = exp( - 1 / 2 * ( ( x1mesh3d ./ sigG1 ) .^ 2 + ...
    ( x2mesh3d ./ sigG2 ) .^ 2  + ( x3mesh3d ./ sigG3 ) .^ 2 ) );
  rho3 = circshift( circshift( circshift( ...
    cCenter, xcgIndShift1, 1 ), xcgIndShift2, 2 ), xcgIndShift3, 3 );
  [c3, p3, px3, py3] = calcCP( rho3, x3, cosphi, sinphi );
  plotCP( c3, p3, px3, py3, x1, x2 )
end
%%
if 0
  cCenter = exp( - 1 / 2 * ( ( x1mesh3d ./ sigG1 ) .^ 2 ) ) ...
    + exp( - 1 / 2 * ( ( x2mesh3d ./ sigG2 ) .^ 2 ) ) ...
    + exp( - 1 / 2 * ( ( x3mesh3d ./ sigG3 ) .^ 2 ) ) ;
  rho3 = circshift( circshift( circshift( ...
    cCenter, xcgIndShift1, 1 ), xcgIndShift2, 2 ), xcgIndShift3, 3 );
  [c3, p3, px3, py3] = calcCP( rho3, x3, cosphi, sinphi );
  plotCP( c3, p3, px3, py3, x1, x2 )
end
%%
% non sep cos
if 0
  xShift1 = 0;
  xShift2 = 0;
  xShift3 = 0;
  rho4 = 1 + cos( 2*pi / l1 .* (x1mesh3d-xShift1) + ...
    2*pi / l2 .* (x2mesh3d-xShift2) + 2*pi / l3 .* (x3mesh3d-xShift3) );
  [c4, p4, px4, py4] = calcCP( rho4, x3, cosphi, sinphi );
  plotCP( c4, p4, px4, py4, x1, x2 )
end
% %% test
% r2 = x1mesh2d .^ 2 + x2mesh2d .^ 2;
% r = sqrt( r2 );
% c4 =  1 + cos( 2*pi / l1 .* r );
% figure()
% imagesc(r2')
% imagesc( c4' )

%% Lorentzian
if 1
  widthTot = 1;
  width1 = 0.5;
  width2 = 0.5;
  width3 = 0.5;
  shift1 = 0;
  shift2 = 0;
  shift3 = pi/4;
  cOffSet = widthTot ./ ( ( (x1mesh2d -shift1) / width1 ) .^ 2 + ...
    ( (x2mesh2d - shift2) / width2 ) .^ 2 ...
    + (widthTot/2).^2 );
  rho5 = 1 ./ ( ( (x1mesh3d -shift1) / width1 ) .^ 2 ...
    + ( (x2mesh3d - shift2) / width2 ) .^ 2 ...
    + ( (x3mesh3d - shift3) / width3 ) .^ 2 ...
    + 1 );
  normRho = trapz_periodic( x1, trapz_periodic( x2, ...
    trapz_periodic( x3, rho5, 3 ), 2), 1);
  rho5 = rho5 / normRho;
  [c5, p5, px5, py5] = calcCP( rho5, x3, cosphi, sinphi );
%   plotCP( c5, p5, px5, py5, x1, x2 )
  figure()
  subplot(2,2,1)
  plot( 1:n2, c5(n1/2+1,:) )
  subplot(2,2,2)
  plot( 1:n2, p5(n1/2+1,:) )
  subplot(2,2,3)
  pcolor( x1, x2, c5' );
  shading interp
  colorbar
  subplot(2,2,4)
  pcolor( x1, x2, p5' );
  shading interp
  colorbar
  
end
%%
function [c, p, px, py] = calcCP( rho, phi, cosphi, sinphi )

c = trapz_periodic( phi, rho, 3 );
px = trapz_periodic( phi, cosphi .* rho, 3 ) ./ c;
py = trapz_periodic( phi, sinphi .* rho, 3 ) ./ c;
p = sqrt( px .^ 2 + py .^ 2 ) ;

end

function plotCP( c, p, px, py, x, y )
[n1, n2] = size( c );
% Set up a index vector so quiver is too crowded
divNumX = 8;
divNumY = 8;
deltaX  = ceil(n1 / divNumX );
deltaY  = ceil(n2 / divNumY);
subInd1 = 1:deltaX:(n1 + 1 - deltaX);
subInd2 = 1:deltaY:(n2 + 1 - deltaX);

figure()
subplot(1,2,1)
pcolor( x, y, c' );
shading interp
colorbar
axis square
title('C')
subplot(1,2,2)
pcolor( x, y, p' );
shading interp
colorbar
axis square
title('P')
ax = gca;
ax.CLim = [0 1];
hold on
quiver(ax, x(subInd1)', y(subInd2)',...
  px(subInd1,subInd2)',...
  py(subInd1,subInd2)', 0,'color',[1,1,1] );
hold off
end
