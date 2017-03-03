plotflag = 0;
l1 = 10;
l2 = 10;
lRod = 1;
legCell = {'Pe = 0';'Pe = 10';'Pe = 20'; 'Pe = 30'; 'Pe = 60'; 'Pe = 80'; 'Pe = 100'; 'Pe = 120'};
files = {
  'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD0_IC1_SM6_t09.04';
  'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD10_IC1_SM6_t01.01';
  'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD20_IC1_SM6_t01.01';
  'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD30_IC1_SM6_t01.01';
  'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD60_IC1_SM6_t04.01';
  'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD80_IC1_SM6_t04.01';
  'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD100_IC1_SM6_t04.06';
  'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD120_IC1_SM6_t04.04';
  };
%  files = {
%    'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD0_IC1_SM6_t09.04';
%    'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD10_IC1_SM6_t01.01';
%    };
%files = {
  %'Hr_Ani1_N646464_Lx10Ly10_bc1.65_vD0_IC1_SM6_t09.04';
  %};
numFiles = length( files );
numGrid = 64;
indsFixed = 1:numGrid;
% allocate
pDistAveSave = zeros( 64, 64, numFiles );
cSave = zeros( 64, 64, numFiles );
%rho2AveSave = zeros( 64, 64, numFiles );
%
pDistPerpAve = zeros( 64, numFiles );
cPerpAve = zeros( 64, numFiles );
%rho2PerpAve = zeros( 64, numFiles );
%
pDistPerpSlice = zeros( 64, numFiles );
cPerpSlice = zeros( 64, numFiles );
%rho2PerpSlice = zeros( 64, numFiles );
%
vD = zeros( numFiles );
bc = zeros( numFiles );
%
cMinMax = zeros(  2, numFiles );
gMinMax = zeros(  2, numFiles );
%rho2MinMax = zeros(  2, numFiles );
axH = zeros(2, numFiles);
% loop over files
for ii = 1: numFiles
  % load
  filename = files{ii};
  path = ['./analyzedfiles/summarized/vDsweep/' filename ];
  load( [path '/params_' filename '.mat'] )
  load( [path '/op_' filename '.mat'] )
  rho = denRecObj.rhoFinal;
  c = C_rec(:,:,end);
  cSlice= c( n1/2 +1, : );
  [~, cInd] = max( cSlice );
  % store there just in case
  vD(ii) = particleObj.vD;
  bc(ii) = systemObj.bc;
  cMinMax(:,ii) = [ min(min( c ) ) max(max( c ) ) ];
  % build grid vectors
  [n1,n2,n3] = size( rho );
  dx1 = l1 ./ n1;
  dx2 = l2 ./ n2;
  x1 = dx1 .* (-n1/2 + 1:n1/2 ); 
  x2 = dx2 .* (-n2/2 + 1:n2/2 ); 
  % calc correlation
  [pDistAve] = pairDistCalcRho(rho,l1,l2,lRod,plotflag);
  gMinMax(:,ii) = [ min(min( pDistAve ) ) max(max( pDistAve ) ) ];
  % two body density
  %[twoBodDenAve] = twoBodyDenCalc(rho,l1,l2,lRod,plotflag);
  %rho2MinMax(:,ii) = [ min(min( twoBodDenAve ) ) max(max( twoBodDenAve ) ) ];  
   % average over parallel direction 
  pDistPerpAve(:,ii) = mean( pDistAve, 1 );
  cPerpAve(:,ii) = circshift( mean( c, 1 ), -cInd + n1/2, 2 );
  %rho2PerpAve(:,ii) = mean( twoBodDenAve, 1 );
  % a slice
  pDistPerpSlice(:,ii) = pDistAve( n1/2 +1, : );
  cPerpSlice(:,ii) = circshift( c( n1/2 +1, : ), -cInd + n1/2, 2 );
  %rho2PerpSlice(:,ii) = twoBodDenAve( n1/2 +1, : );
  % save just in case
  pDistAveSave(:,:,ii) = pDistAve;
  cSave(:,:,ii) = circshift( c, -cInd + n1/2, 2 );
  %rho2AveSave(:,:,ii) = twoBodDenAve;
end

%% get min max
minC = min( cMinMax(1,:) );
maxC = max( cMinMax(2,:) );
minG = min( gMinMax(1,:) );
maxG = max( gMinMax(2,:) );
%minRho2 = min( rho2MinMax(1,:) );
%maxRho2 = max( rho2MinMax(2,:) );
% ave plot
numSubs = 2;
figure()
subplot( 1,numSubs,1)
hold on
ax1ave = gca;
subplot( 1,numSubs,2)
hold on
ax2ave = gca;
%subplot( 1,numSubs,3)
%hold on
%ax3ave = gca;
% slice plot
figure()
subplot( 1,numSubs,1)
hold on
ax1slice = gca;
subplot( 1,numSubs,2)
hold on
ax2slice = gca;
%subplot( 1,numSubs,3)
%hold on
%ax3slice = gca;
% plot
for ii = 1:numFiles
  figure()
  % c
  subplot(1,numSubs,1)
  imagesc( x1, x2, cSave(:,:,ii) )
  axis square
  colorbar
  title( [ 'concentration c(x,y) $$\langle c \rangle = $$ ' ...
    num2str(bc(ii)) ' $$ v_D = $$ ' num2str(vD(ii)) ] )
  xlabel('x'); ylabel('y')
  axT = gca;
  axT.CLim = [minC maxC];
  % dist
  subplot(1,numSubs,2)
  imagesc( x1, x2, pDistAveSave(:,:,ii) )
  axis square
  colorbar
  title( [ 'pair distribution g(x,y) $$\langle c \rangle = $$ '...
    num2str(bc(ii)) ' $$ v_D = $$ ' num2str(vD(ii))] )
  xlabel('x'); ylabel('y')
  axT = gca;
  axT.CLim = [0 maxG];
  % rho 2
  %subplot(1,numSubs,3)
  %imagesc( x1, x2, rho2AveSave(:,:,ii) )
  %axis square
  %colorbar
  %title( [ '$$\rho^{(2)}(x,y) $$ $$ \langle c \rangle = $$ '...
    %num2str(bc(ii)) ' $$ v_D = $$ ' num2str(vD(ii))] )
  %xlabel('x'); ylabel('y')
  %axT = gca;
  %axT.CLim = [minRho2 maxRho2];
  % Plot perpedicular ave
  plot(ax1ave, x1, cPerpAve(:,ii) )
  plot(ax2ave, x1, pDistPerpAve(:, ii) )
  %plot(ax3ave, x1, rho2PerpAve(:, ii) )
  % Plot perpedicular slice
  plot(ax1slice, x1, cPerpSlice(:,ii) )
  plot(ax2slice, x1, pDistPerpSlice(:, ii) )
  %plot(ax3slice, x1, rho2PerpSlice(:, ii) )
end
ax1ave.XLabel.String = '$$ r_{\perp} $$'; 
ax2ave.XLabel.String = '$$ r_{\perp} $$';
%ax3ave.XLabel.String = '$$ r_{\perp} $$';
ax1ave.YLabel.String = '$$ c( r_{\perp} ) $$'; 
ax2ave.YLabel.String = '$$ g( r_{\perp} ) $$';
%ax3ave.YLabel.String = '$$ \rho^{(2)} ( r_{\perp} ) $$';
ax1slice.XLabel.String = '$$ r_{\perp} $$'; 
ax2slice.XLabel.String = '$$ r_{\perp} $$';
%ax3slice.XLabel.String = '$$ r_{\perp} $$';
ax1slice.YLabel.String = '$$ c( r_{\perp} ) $$'; 
ax2slice.YLabel.String = '$$ g( r_{\perp} ) $$';
%ax3slice.YLabel.String = '$$ \rho^{(2)}( r_{\perp} ) $$';
title(ax1ave, 'ave');
title(ax1slice, 'slice');
% fix lims
ax2slice.YLim = [0 maxG]; ax2ave.YLim = [0 maxG];
% legends
hl = legend(ax1ave,legCell );
hl.Location = 'best';
hl = legend(ax1slice,legCell );
hl.Location = 'best';