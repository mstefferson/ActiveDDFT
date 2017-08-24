% quiver/pcolor plots x axis going across columns. Increasing column, increasing
% x. plots just like a matrix looks. For y axis, it flips the rows. i.e.,
% row 1 is plotted at the bottom.
%
set(0,'DefaultFigureWindowStyle','normal')

n1 = 128;
n2 = 128;
n3 = 128;
l1 = 10;
l2 = 10;
l3 = 2*pi;

numTicks = 8;

dN1 = round( n1 / numTicks );
dN2 = round( n1 / numTicks );

dNind1 = 1:dN1:n1;
dNind2 = 1:dN2:n2;

n1Vec = 1:n1;
n2Vec = 1:n2;
n3Vec = 1:n3;

[n2mesh, n1mesh, n3mesh] = meshgrid( n2Vec, n1Vec, n3Vec );

x = n1Vec * l1 / n1;
y = n2Vec * l1 / n1;
phi = n3Vec * l3 / n3;

[ymesh, xmesh, phimesh] = meshgrid( y, x, phi );
% Colors
% cMat = lines(4);
% cMat = cMat( [1 2 4], : );
cMat = copper(3);
% cMat = [ 1 0 0; 1 1 0; 0 0 0];
% px = cos( 2*pi / n1 * n1mesh) .* ones(n1,n2);
% py = cos( 2*pi / n1 * n2mesh) .* ones(n1,n2);
% px = n1Vec .* ones(n1,n2);

rho = cos( phimesh) .^ 2 .*  ( ( xmesh - 3* l1/4) .^2  +  ( ymesh -  l2 / 6 ) .^2 );
% py = zeros(n1,n2);
c = trapz_periodic( phi, rho, 3 );

%%
figure
% handaxes1 = axes('Position', [0 0 1 1]);
hold on
pcolor( x, y, c' )
shading interp
colorbar
xlabel('x');
ylabel('y');
axis square
ax1 = gca;
ax1.XLim = [x(1) x(end)];
ax1.YLim = [y(1) y(end)];
% plot dots
%%
dot1x = n1/4;
dot1y = 3 * n2/4;
dot2x = n1/2;
dot2y = n2/2;
dot3x = 3*n1/4;
dot3y = n2/4;
p = plot( x(dot1x), y(dot1y),'x','MarkerSize',20 );
p.Color = cMat(1,:);
p = plot( x(dot2x), y(dot2y),'x','MarkerSize',20 );
p.Color = cMat(2,:);
p = plot( x(dot3x), y(dot3y),'x','MarkerSize',20 );
p.Color = cMat(3,:);
%%
hold off
posVec = ax1.Position;
x0 = ax1.Position(1);
y0 = ax1.Position(2);
w0 = ax1.Position(3);
h0 = ax1.Position(4);

% Place second set of axes on same plot
%%
x1 = x0 + w0 / 2;
y1 = y0 + h0 / 2;
w1 = w0 / 4;
h1 = h0 / 4;
ax2 = axes('Position', [x1 y1 w1 h1]);
hold on
p = plot(phi, reshape( rho( dot1x, dot1y, : ), [1 n3] ) );
p.Color = cMat(1,:);
p = plot(phi, reshape( rho( dot2x, dot2y, : ), [1 n3] ) );
p.Color = cMat(2,:);
p = plot(phi, reshape( rho( dot3x, dot3y, : ), [1 n3] ) );
p.Color = cMat(3,:);
ax2.XLim = [ phi(1) phi(end) ]; 
axLabColor = [0.85 0 0];
ax2.XColor = axLabColor;
ax2.YColor = axLabColor;
ax2.XLabel.Color = axLabColor;
ax2.YLabel.Color = axLabColor;
hold off
axis square
set(ax2, 'Box', 'off')
xlabel('$$ \phi $$')
ylabel('$$ f( \phi ) $$')
%%
x1 = x0 + 0.62 * w0 ;
y1 = y0 + 0.67 * h0;
w1 = w0 / 4;
h1 = h0 / 4;
ax2.Position =  [x1 y1 w1 h1];


