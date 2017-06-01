% Make a phase diagram of the runs
% plotSoftShoulderRunPhase(T,saveMe)
function plotSoftShoulderRunPhase(T,saveMe)
% make sure you are only looking at fixed box and grid length
if length( unique( T.l1 ) ) > 1 || length( unique( T.l2 ) ) > 1
  fprintf('Box length is changing\n');
  error('Box length is changing');
end
if length( unique( T.n1 ) ) > 1 || length( unique( T.n2 ) ) > 1
  fprintf('Number of gridpoints is changing\n');
  error('Number of gridpoints is changing');
end
% set up figure
figure()
numColors = 3;
colorArray = colormap(['lines(' num2str(numColors) ')']);
titlstr = [ 'DDFT soft shoulder phase map L = ' num2str( unique( T.l1 ) ) ' N = ' num2str( unique( T.n1 ) ) ];
title( titlstr );
xlabel( '$$ \rho $$' ); ylabel( '$$ a $$' );
cRange = [0 max( 4, max(T.c) )];
aRange = [0 max( 2, max(T.longE2) )];
ax = gca;
ax.XTick = cRange(1):(cRange(2)-cRange(1))/4:cRange(2);
ax.YTick = aRange(1):(aRange(2)-aRange(1))/4:aRange(2);
ax.XLim = [cRange(1) cRange(2)];
ax.YLim = [aRange(1) aRange(2)];
legcell = {'failure' 'liquid', 'crystal A', 'crystal B' };
% plot failures
scatter( T.c( T.broken == 1 ), T.longE2( T.broken == 1 ), [], [1 0 0],'x' );
hold on
% liquid
logicalCond = T.broken == 0 & strcmp( T.spatPhase, 'Liquid' );
scatter( T.c( logicalCond ), T.longE2( logicalCond ), [], colorArray(1,:),'filled' );
% crystal A
logicalCond = T.broken == 0 & strcmp( T.spatPhase, 'Crystal A' );
scatter( T.c( logicalCond ), T.longE2( logicalCond ), [], colorArray(2,:),'filled' );
% crystal B
logicalCond = T.broken == 0 & strcmp( T.spatPhase, 'Crystal B' );
scatter( T.c( logicalCond ), T.longE2( logicalCond ), [], colorArray(3,:),'filled' );
% place legend
legend(legcell,'location', 'best')
% save it
if saveMe
  savefig( gcf, 'phaseMap' )
end
