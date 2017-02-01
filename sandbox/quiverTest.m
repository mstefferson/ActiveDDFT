%%
% theta11 = pi/4;
% theta12 = -pi/2;
% theta13 = pi/2;
% theta21 = -pi/4;
% theta22 = pi;
% theta23 = pi/4;
% % 

thetaAll = 0;
theta11 = thetaAll;
theta12 = thetaAll;
theta13 = thetaAll;
theta21 = thetaAll;
theta22 = thetaAll;
theta23 = thetaAll;

Ax = [ cos(theta11) cos(theta12) cos(theta13); ...
      cos(theta21) cos(theta22) cos(theta23)];
Ay = [ sin(theta11) sin(theta12) sin(theta13); ...
      sin(theta21) sin(theta22) sin(theta23)];
    
theta11actMat = theta11 - pi/2 ;
theta12actMat = theta12 - pi/2 ;
theta13actMat = theta13 - pi/2 ;
theta21actMat = theta21 - pi/2 ;
theta22actMat = theta22 - pi/2 ;
theta23actMat = theta23 - pi/2 ;

AxAct = [ cos(theta11actMat) cos(theta12actMat) cos(theta13actMat); ...
      cos(theta21actMat) cos(theta22actMat) cos(theta23actMat)];
AyAct = [ sin(theta11actMat) sin(theta12actMat) sin(theta13actMat); ...
      sin(theta21actMat) sin(theta22actMat) sin(theta23actMat)];
%    
% Amps = [ 1 2 3; 1.5 2.5 3.5];
Amps = [ 3 2 1; 1 3 2 ];

disp(Ax); disp(Ay);disp(Amps);

columns = 1:3;
rows = 1:2;

%
% testing flips and rotation
figure()
subplot(2,2,1)
imagesc(columns,rows,Amps);
colorbar
ax = gca;
ax.YDir = 'rev';
ax.NextPlot = 'add';
quiver(columns,rows,AxAct,-AyAct)
axis square
xlabel('columns'); ylabel('rows');
title([ 'phi M.C.S. actual matrix phi = ' num2str(theta11) ] )

subplot(2,2,2)
imagesc(rows,columns, rot90(Amps) );
colorbar
ax = gca;
ax.YDir = 'rev';
ax.NextPlot = 'add';
quiver(rows, columns, rot90( Ax ), rot90( -Ay ) )
axis square
ax.YTick = columns;
ax.YTickLabel = flip(columns);
xlabel('rows'); ylabel('columns');
title([ 'phi M.C.S. rotated matrix phi = ' num2str(theta11) ] )

subplot(2,2,3)
imagesc(columns,rows,Amps);
colorbar
ax = gca;
ax.YDir = 'rev';
ax.NextPlot = 'add';
quiver(columns,rows,Ax,-Ay)
axis square
xlabel('columns'); ylabel('rows');
title([ 'phi N.C.S. actual matrix phi = ' num2str(theta11) ] )

subplot(2,2,4)
imagesc(rows,columns, rot90(Amps) );
colorbar
ax = gca;
ax.YDir = 'rev';
ax.NextPlot = 'add';
quiver(rows, columns, rot90( -Ay ), rot90( -Ax ) )
axis square
ax.YTick = columns;
ax.YTickLabel = flip(columns);
xlabel('rows'); ylabel('columns');
title([ 'phi N.C.S. rotated matrix phi = ' num2str(theta11) ] )

%% Test image Quiver
imageQuiver( Amps, Ax, Ay, 1)
imageQuiver( Amps, Ax, Ay, 2)
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
