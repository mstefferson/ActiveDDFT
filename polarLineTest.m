% quiver/pcolor plots x axis going across columns. Increasing column, increasing
% x. plots just like a matrix looks. For y axis, it flips the rows. i.e.,
% row 1 is plotted at the bottom.
%
n1 = 128;
n2 = 128;

numTicks = 8;

dN1 = round( n1 / numTicks );
dN2 = round( n1 / numTicks );

dNind1 = 1:dN1:n1;
dNind2 = 1:dN2:n2;

n1Vec = 1:n1;
n2Vec = 1:n2;


[n2mesh, n1mesh] = meshgrid( n2Vec, n1Vec );

% px = cos( 2*pi / n1 * n1mesh) .* ones(n1,n2);
% py = cos( 2*pi / n1 * n2mesh) .* ones(n1,n2);
% px = n1Vec .* ones(n1,n2);
px = n1Vec' .* ones(n1,n2);
py = zeros(n1,n2);
% py = zeros(n1,n2);
amp = px .^ 2  + py .^ 2;

figure(1)
ax = subplot(2,2,1);
pcolor( n2Vec, n1Vec, amp )
shading interp
ax.XLim = [n2Vec(1) n2Vec(end)];
ax.YLim = [n1Vec(1) n1Vec(end)];

ax = subplot(2,2,2);
quiver(n2Vec, n1Vec, px, py);
ax.XLim = [n2Vec(1) n2Vec(end)];
ax.YLim = [n1Vec(1) n1Vec(end)];

ax = subplot(2,2,3);
quiver( n2Vec(dNind2), n1Vec(dNind1), ...
  px(dNind1, dNind2), py(dNind1, dNind2));
ax.XLim = [n2Vec(1) n2Vec(end)];
ax.YLim = [n1Vec(1) n1Vec(end)];
ax.XLim = [n2Vec(1) n2Vec(end)];
ax.YLim = [n1Vec(1) n1Vec(end)];

ax = subplot(2,2,4);
pcolor( n2Vec, n1Vec, amp )
shading interp
hold on
quiver( n2Vec(dNind2), n1Vec(dNind1), ...
  px(dNind1, dNind2), py(dNind1, dNind2));
ax.XLim = [n2Vec(1) n2Vec(end)];
ax.YLim = [n1Vec(1) n1Vec(end)];
hold off

figure(2)
ax = subplot(2,2,1);
imagesc( n2Vec, n1Vec, amp )
ax.XLim = [n2Vec(1) n2Vec(end)];
ax.YLim = [n1Vec(1) n1Vec(end)];

ax = subplot(2,2,2);
quiver(n2Vec, n1Vec, px, py);
ax.XLim = [n2Vec(1) n2Vec(end)];
ax.YLim = [n1Vec(1) n1Vec(end)];

ax = subplot(2,2,3);
quiver( n2Vec(dNind2), n1Vec(dNind1), ...
  px(dNind1, dNind2), py(dNind1, dNind2));
ax.XLim = [n2Vec(1) n2Vec(end)];
ax.YLim = [n1Vec(1) n1Vec(end)];
ax.XLim = [n2Vec(1) n2Vec(end)];
ax.YLim = [n1Vec(1) n1Vec(end)];

ax = subplot(2,2,4);
imagesc( n2Vec, n1Vec, amp )
shading interp
hold on
quiver( n2Vec(dNind2), n1Vec(dNind1), ...
  px(dNind1, dNind2), py(dNind1, dNind2));
ax.XLim = [n2Vec(1) n2Vec(end)];
ax.YLim = [n1Vec(1) n1Vec(end)];
hold off

