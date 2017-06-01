%%
clear
currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );
plotMe = 0;
l = 100;
lrE1 = 1;
lrL1 = 1;
lrL2 = 1.855;
params = ...
  [3.00 0.50;
   2.25 1.75;
   3.75 1.75;
   3.75 1.20;
   3.75 0.75;
   4.00 0.25];
numParams = size( params, 1 );
n = ceil( 10 * l  / pi + 2 );
n = n + mod(n,2);

for ii = 1:numParams
  % run disp
  cTemp = params(ii,1);
  aTemp = params(ii,2);
  paramVec = [n, n, l, l, lrE1, aTemp, lrL1, lrL2, cTemp];
  tic
  [disp] = dispersionSoftShoulder( paramVec, plotMe);
  toc
  % plot
  figure()
  xTemp = disp.k1pos;
  yTemp = disp.omegaSlice;
  p = plot( xTemp, yTemp );
  p.LineWidth = 2;
  hold on
  numPoints = length( xTemp  );
  yTemp = zeros( 1, numPoints );
  p = plot( xTemp, yTemp,'--' );
  p.LineWidth = 2;
  ax = gca;
  ax.XLim = [ 0 max(xTemp) ];
  scalFact  = 0.25; %for limits
  omegaMin =  min( disp.omegaSlice );
  omegaMax =  max( disp.omegaSlice );
  ax.YLim = [ omegaMin + omegaMin * scalFact, max( omegaMax + omegaMax * scalFact, 10 ) ];
  xlabel(' $$ kR $$ ' ); ylabel(' $$ \omega $$ ' );
  axis square
  titstr = [ '$$\omega(k)$$: reg' num2str(ii,'%d') ' L = ' num2str(l,'%d')...
     ' a = ' num2str( aTemp, '%.2f' ) ' c = ' num2str( cTemp, '%.2f' ) ];
  title(titstr)
  fprintf( 'omega max =  %f, kmax = %f, %s\n', ...
    disp.omegaMax, disp.kPeakMax, disp.phase );
  figname = ['dispRel_reg' num2str( ii, '%d' ) ...
    '_a' num2str( aTemp, '%.2f' ) '_c' num2str( cTemp, '%.2f' ) ...
    '_l' num2str(l,'%d') '.fig' ];
  savefig( gcf, figname ) 
end
