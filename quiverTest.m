%%
theta11 = pi/4;
theta12 = -pi/2;
theta13 = pi/2;
theta21 = -pi/4;
theta22 = pi;
theta23 = pi/4;

Ax = [ cos(theta11) cos(theta12) cos(theta13); ...
      cos(theta21) cos(theta22) cos(theta23)];
Ay = [ sin(theta11) sin(theta12) sin(theta13); ...
      sin(theta21) sin(theta22) sin(theta23)];
Amps = [ 1 2 3; 4 5 6];

disp(Ax); disp(Ay);disp(Amps);

columns = 1:3;
rows = 1:2;
%%
imageQuiver( Amps, Ax, Ay, 1)
imageQuiver( Amps, Ax, Ay, 2)
%%
%% Test imagesc
figure()
subplot(2,1,1)
i1 = imagesc(Amps);
title('no switch (rev) exact copy')
%
subplot(2,1,2)
i2 = imagesc(Amps);
ax  = gca;
ax.YDir = 'normal';
title('switched (normal) flipped')
%% Test pquiver
figure()
subplot(2,1,1)
ax1 = gca;
p1 = quiver( Ax, Ay, 0.5);
title('no switch (normal) flipped')
%
subplot(2,1,2)
p2 = quiver( Ax, -Ay, 0.5);
ax2  = gca;
ax2.YDir = 'rev';
title('switched (rev) exact copy')

%% Method 1: exact copy of matrix

%figure()
% Matlab way M( y , x )
figure()
ax = gca;
ax.YDir = 'rev';
imagesc( Amps );
colorbar
ax = gca;
% Flip y direction so y axis counts ( 1, 2, 3,...) to math p quiver
%ax.YDir = 'normal';
hold
p = quiver( Ax, -Ay, 0.5);
p.Color = [1 1 1];
xlabel('columns'); ylabel('rows');
ax.XTick      = columns;
ax.YTick      = rows;
ax.XTickLabel = columns;
ax.YTickLabel = rows;
title('Method 1: exact Copy')
disp(Ax); disp(Ay);disp(Amps);
%%

%% Method 2: flip rows

%figure()
% Matlab way M( y , x )
figure()
imagesc( Amps );
colorbar
ax = gca;
% Flip y direction so y axis counts ( 1, 2, 3,...) to math p quiver
ax.YDir = 'normal';
hold
p = quiver( Ax, Ay, 0.5);
p.Color = [1 1 1];
xlabel('columns'); ylabel('rows');
ax.XTick      = columns;
ax.YTick      = rows;
ax.XTickLabel = columns;
ax.YTickLabel = rows;
title('Method 2: Flips rows')
disp(Ax); disp(Ay);disp(Amps);

%% Method 3: Invert rows and columns

% % Mike way M(x,y)
% figure()
% imagesc( Amps' );
% ax = gca;
% % Flip y direction so y axis counts ( 1, 2, 3,...) to math p quiver
% ax.YDir = 'normal';
% hold
% p = quiver(Ax', Ay',0.5);
% p.Color = [1 1 1];
% xlabel('rows'); ylabel('columns');
% title('Mike way: right')
% ax.XTick      = rows;
% ax.YTick      = columns;
% ax.XTickLabel = rows;
% ax.YTickLabel = columns;
% disp(Ax'); disp(Ay');disp(Amps);