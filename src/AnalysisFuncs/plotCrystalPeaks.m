% Plot k peaks vs time
function plotCrystalPeaks( cFt, k1, k2, timeRec, systemObj, particleObj, saveMe)
% add paths for now
addpath( genpath( './src' ) )
% run disperion
paramVec = [ systemObj.n1, systemObj.n2, systemObj.l1, systemObj.l2, ...
  particleObj.lrEs1, particleObj.lrEs2, ...
  particleObj.lrLs1, particleObj.lrLs2, systemObj.c ];
disp = dispersionSoftShoulder( paramVec, 1 );
fig0 = gcf;
% get size and center
[n1, n2, nt] = size( cFt );
k1center = n1 / 2 + 1;
% k1PosLength = n1 / 2;
k2center = n2 / 2 + 1;
% k2PosLength = n1 / 2;
% get positive kx, ky vectors
k1Pos = k1( k1center:n1 );
k2Pos = k2( k2center:n2 );
% Look at positive kx, ky
cFtPos = cFt(k1center:n1, k2center:n2,:);
% Get rid of noise
noiseEps = 1e-10;
cFtPos( cFtPos < noiseEps ) = 0;
% find peaks
ftTemp = cFtPos(:,:,end);
ftTemp( ftTemp < ftTemp(1,1) / 100 ) = 0;
[~,k1PeakInds] = findpeaks( ftTemp(:,1) );
[~,k2PeakInds] = findpeaks( ftTemp(1,:) );
% Add zero so there are no bugs
k1PeakInds = [1 k1PeakInds'];
k2PeakInds = [1 k2PeakInds];
numk1peaks = length( k1PeakInds );
numk2peaks = length( k2PeakInds );
% peaks in k-space
dPeak = 12;
fig1 = figure();
colorArray = viridis( nt );
legcell = cell( nt, 1 );
for ii = 1:nt
  legcell{ii} = [ '$$ t = ' num2str( timeRec(ii) ) ' $$'];
end
% k1 peaks
subplot(1,2,1);
inds = 2: min( max(k1PeakInds) + dPeak, n1 );
c2plot =  reshape( cFtPos( inds , 1, :), [length(inds), nt] );
x2plot = repmat( k1Pos(inds)', [1 nt] );
p = plot( x2plot, c2plot, 'Linewidth', 2 );
for ii = 1:nt
  p(ii).Color = colorArray(ii,:);
end
xlabel( ' $$ k_1 $$ '); ylabel('Amplitude');
title( '$$ k_1 $$ modes')
leg = legend(legcell);
leg.Interpreter = 'latex';
% k2 peaks
ax = subplot(1,2,2);
inds = 2: min( max(k2PeakInds) + dPeak, n2 );
c2plot =  reshape( cFtPos( 1, inds, :), [length(inds), nt] );
x2plot = repmat( k2Pos(inds)', [1 nt] );
p = plot( x2plot, c2plot, 'Linewidth', 2 );
for ii = 1:nt
  p(ii).Color = colorArray(ii,:);
end
xlabel( ' $$ k_2 $$ '); ylabel('Amplitude');
title( '$$ k_2 $$ modes')
ax.XAxis.TickLabelFormat = '%,.2f';
leg = legend(legcell);
leg.Interpreter = 'latex';
% peaks in time
fig2 = figure();
% k1 peaks vs time
colorArray = viridis( numk1peaks );
subplot(2,2,1);
legcell = cell( numk1peaks, 1 );
for ii = 1:numk1peaks
  legcell{ii} = ['$$ k_1 =  ' num2str( k1Pos( k1PeakInds(ii) ),'%.1f' ) ' $$'];
end
tempFt2plot = reshape( cFtPos( k1PeakInds , 1, : ), [numk1peaks nt] )';
p = plot( timeRec, tempFt2plot, 'LineWidth', 2 );
for ii = 1:numk1peaks
  p(ii).Color = colorArray(ii,:);
end
xlabel( ' $$ t $$ '); ylabel('Amplitude');
title( '$$ k_1 $$ modes')
leg = legend(legcell);
leg.Interpreter = 'latex';
% k2 peaks vs time
colorArray = viridis( numk2peaks );
ax = subplot(2,2,2);
legcell = cell( numk2peaks, 1 );
for ii = 1:numk2peaks
  legcell{ii} = ['$$ k_2 = ' num2str( k2Pos( k2PeakInds(ii) ),'%.1f' ) ' $$'];
end
tempFt2plot = reshape( cFtPos( 1, k2PeakInds, : ), [numk2peaks nt] )' ;
p = plot( timeRec, tempFt2plot, 'LineWidth', 2  );
for ii = 1:numk2peaks
  p(ii).Color = colorArray(ii,:);
end
xlabel( ' $$ t $$ '); ylabel('Amplitude');
title( '$$ k_2 $$ modes')
ax.YAxis.TickLabelFormat = '%,.2f';
leg = legend(legcell);
leg.Interpreter = 'latex';
% unstable peaks plotyy
maxLinDispPeaks = disp.kAllPeakInds;
if length( maxLinDispPeaks ) < 3
  maxLinDispPeaks = [ maxLinDispPeaks,...
    ones(1,3-length( maxLinDispPeaks ) ) ];
end
subplot(2,2,3);
if length(k1PeakInds) >= length(k2PeakInds)
  pickedkDir = 'k1';
else
  pickedkDir = 'k2';
end
% peak 1
if length(k1PeakInds) >= length(k2PeakInds)
  tempFt2plot1 = reshape( cFtPos( maxLinDispPeaks(2), 1, : ), [1 nt] )' ;
else
  tempFt2plot1 = reshape( cFtPos(1, maxLinDispPeaks(2), : ), [1 nt] )' ;
end
% get stabillity info
if maxLinDispPeaks(2) == disp.kPeakMaxInd
  plotk1peak = 'max unstable';
  dashtype1 = '-';
elseif any( maxLinDispPeaks(2) == disp.kPeaksInds )
  plotk1peak = 'unstable';
  dashtype1 = '--';
else
  plotk1peak = 'stable';
  dashtype1 = ':';
end
% peak 2
if length(k1PeakInds) >= length(k2PeakInds)
  tempFt2plot2 = reshape( cFtPos(maxLinDispPeaks(3), 1, : ), [1 nt] )' ;
else
  tempFt2plot2 = reshape( cFtPos(1, maxLinDispPeaks(3), : ), [1 nt] )' ;
end
% get stabillity info
if maxLinDispPeaks(3) == disp.kPeakMaxInd
  plotk2peak = 'max unstable';
  dashtype2 = '-';
elseif any( maxLinDispPeaks(3) == disp.kPeaksInds)
  plotk2peak = 'unstable';
  dashtype2 = '--';
else
  plotk2peak = 'stable';
  dashtype2 = ':';
end
[~, p1, p2 ] = plotyy( timeRec, tempFt2plot1, timeRec, tempFt2plot2);
p1.LineStyle = dashtype1;
p1.LineWidth = 2;
p2.LineStyle = dashtype2;
p2.LineWidth = 2;
titlstr = ['Lin stability peaks. k dir: ' pickedkDir ];
title(titlstr)
xlabel('t');ylabel('Amplitude');
legcell = { ['$$ k_{p,l} $$ : ' plotk1peak ], ...
  [ '$$ k_{p,u} $$ : ' plotk2peak ] };
leg = legend(legcell,'location','best');
leg.Interpreter = 'latex';
% unstable peaks plot together
subplot(2,2,4);
% normalize if possible
if tempFt2plot1(1) ~= 0 && tempFt2plot2(1) ~= 0
  tempFt2plot1 = tempFt2plot1 ./ tempFt2plot1(1);
  tempFt2plot2 = tempFt2plot2 ./ tempFt2plot2(1);
  ylab = 'Amplitude / A(0)';
else
  ylab = 'Amplitude';
end
p = plot( timeRec, tempFt2plot1, timeRec, tempFt2plot2 );
ax = gca;
limFact = 1.61;
ax.YLim = [ 0 limFact * min( [ max( tempFt2plot1 ) max( tempFt2plot2 ) ] ) ];
p(1).LineStyle = dashtype1;
p(1).LineWidth = 2;
p(2).LineStyle = dashtype2;
p(2).LineWidth = 2;
titlstr = [ ' $$ k_{p,l} = ' num2str( maxLinDispPeaks(2) )...
  '$$ $$ k_{p,u} = ' num2str( maxLinDispPeaks(3) ) '$$'];
title(titlstr)
xlabel('t');ylabel(ylab);
legcell = { ['$$ k_{p,l} $$ : ' plotk1peak ], ...
  [ '$$ k_{p,u} $$ : ' plotk2peak ] };
leg = legend(legcell,'location','best');
leg.Interpreter = 'latex';
% all peaks
fig3 = figure();
colormap(fig3, viridis )
subplot(1,2,1)
imagesc( k2, k1, cFt(:,:,end)' )
colorbar
axis square
ax = gca;
ax.YDir = 'normal';
ax.YLabel = fliplr(ax.YLabel);
subplot(1,2,2)
kSub = 1:min( max([ k1PeakInds k2PeakInds ]) + dPeak, min( n1, n2 ) );
imagesc( k2Pos(kSub), k1Pos(kSub), ...
  cFtPos(kSub,kSub,end)' );
colorbar
axis square
ax = gca;
ax.YDir = 'normal';
ax.YLabel = fliplr(ax.YLabel);
% save it
if saveMe
  fig0name = 'kAmpsDisp';
  fig1name = 'kAmpsVsK';
  fig2name = 'kAmpsVsT';
  fig3name = 'kAmpsMap';
  savefig( fig0, fig0name )
  saveas( fig0, fig0name, 'jpg' )
  savefig( fig1, fig1name )
  saveas( fig1, fig1name, 'jpg' )
  savefig( fig2, fig2name )
  saveas( fig2, fig2name, 'jpg' )
  savefig( fig3, fig3name )
  saveas( fig3, fig3name, 'jpg' )
end
