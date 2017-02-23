%% All files
plotOPs = 0;
plotMax = 1;
plotMin = 0;
plotAverage = 1;
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
numFiles = length( files );
numGrid = 64;
indsFixed = 1:numGrid;
% allocate
vD = zeros( numFiles );
% c
cAtMaxVec = zeros( numGrid, numFiles );
cMax = zeros( numFiles, 1 );
cAtMinVec = zeros( numGrid, numFiles );
cMin = zeros( numFiles, 1 );
c2bodyVec = zeros( numGrid, numFiles );
cFinal = zeros( numGrid, numGrid, numFiles );
% p
pAtMaxVec = zeros( numGrid, numFiles );
pMax = zeros( numFiles, 1 );
pAtMinVec = zeros( numGrid, numFiles );
pMin = zeros( numFiles, 1 );
p2bodyVec = zeros( numGrid, numFiles );
pFinal = zeros( numGrid, numGrid, numFiles );
% n
nAtMaxVec = zeros( numGrid, numFiles );
nMax = zeros( numFiles, 1 );
nAtMinVec = zeros( numGrid, numFiles );
nMin = zeros( numFiles, 1 );
n2bodyVec = zeros( numGrid, numFiles );
nFinal = zeros( numGrid, numGrid, numFiles );
% Loop
for ii = 1: numFiles
  filename = files{ii};
  path = ['./analyzedfiles/summarized/vDsweep/' filename ];
  load( [path '/op_' filename '.mat'] )
  load( [path '/params_' filename '.mat'] )
  vD(ii) = particleObj.vD;
  % get position vector
  x = 0 : systemObj.Ly / systemObj.Ny : systemObj.Ly - systemObj.Ly / systemObj.Ny;
  % Get steady state distribution
  cTemp = C_rec(:,:,end) ./ pi;
  pTemp = POP_rec(:,:,end);
  nTemp = NOP_rec(:,:,end);
  cFinal(:,:,ii) = cTemp;
  pFinal(:,:,ii) = pTemp;
  nFinal(:,:,ii) = nTemp;
  % get a slice that is averaged over non-varying direction
  cSlice = mean( cTemp, 1 );
  pSlice = mean( pTemp, 1 );
  nSlice = mean( nTemp, 1 );
  % find max
  [cMaxTemp, maxInd] = max( cSlice );
  [pMaxTemp] = max( pSlice );
  [nMaxTemp] = max( nSlice );
  [cMinTemp] = min( cSlice );
  [pMinTemp] = min( pSlice );
  [nMinTemp] = min( nSlice );
  % circ and save
  % c
  cSliceMax = circshift( cSlice, -maxInd );
  cSliceMin = circshift( cSlice, -minInd );
  cAtMaxVec(:,ii) = cSliceMax;
  cMax(ii) = cMaxTemp;
  cAtMinVec(:,ii) = cSliceMin;
  cMin(ii) = cMinTemp;
  % p
  pSliceMax = circshift( pSlice, -maxInd );
  pSliceMin = circshift( pSlice, -minInd );
  pAtMaxVec(:,ii) = pSliceMax;
  pMax(ii) = pMaxTemp;
  pAtMinVec(:,ii) = pSliceMin;
  pMin(ii) = pMinTemp;
  % p
  nSliceMax = circshift( nSlice, -maxInd );
  nSliceMin = circshift( nSlice, -minInd );
  nAtMaxVec(:,ii) = nSliceMax;
  nMax(ii) = nMaxTemp;
  nAtMinVec(:,ii) = nSliceMin;
  nMin(ii) = nMinTemp;
  % average
  for jj = 1:numGrid
    deltaInd = jj-1;
    inds = indsFixed + deltaInd;
    inds = mod( inds - 1, numGrid ) + 1;
    c2bodyVec(jj,ii) = mean( cSliceMax( 1:numGrid ) .* cSliceMax( inds ) ); 
    p2bodyVec(jj,ii) = mean( pSliceMax( 1:numGrid ) .* pSliceMax( inds ) ); 
    n2bodyVec(jj,ii) = mean( nSliceMax( 1:numGrid ) .* nSliceMax( inds ) ); 
  end
end
%%
% plot it
% OPs
if plotOPs
  for ii = 1: numFiles
    figure()
    % c
    subplot(1,3,1)
    imagesc( x,x, cFinal(:,:,ii) )
    axis square
    title(['C(x,y) vD = ' num2str( vD(ii) )] );
    xlabel('x');
    ylabel('y');
    colorbar;
    % p
    subplot(1,3,2)
    imagesc( x,x, pFinal(:,:,ii) )
    axis square
    title(['P(x,y) vD = ' num2str( vD(ii) )] );
    xlabel('x');
    ylabel('y');
    colorbar;
    % n
    subplot(1,3,3)
    imagesc( x,x, nFinal(:,:,ii) )
    axis square
    title(['N(x,y) vD = ' num2str( vD(ii) )] );
    xlabel('x');
    ylabel('y');
    colorbar;
  end
end
%% plot max
if plotMax
  figure()
  % c
  % max
  subplot(1,3,1)
  hold on
  for ii = 1:numFiles
    plot( x, cAtMaxVec(:,ii) .* cMax(ii) )
  end
  title(' $$ c^{(2)}(r_1,r_2) $$; $$ r_1 $$ at peak');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ c (r_0) c( r_{\perp} ) $$');
  axis square
  hl= legend(legCell );
  % p
  % max
  subplot(1,3,2)
  hold on
  for ii = 1:numFiles
    plot( x, pAtMaxVec(:,ii) .* pMax(ii) )
  end
  title(' $$ p^{(2)}(r_1,r_2) $$; $$ r_1 $$ at peak');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ p (r_0) p( r_{\perp} ) $$');
  axis square
  hl= legend(legCell );
  % n
  % max
  subplot(1,3,3)
  hold on
  for ii = 1:numFiles
    plot( x, nAtMaxVec(:,ii) .* nMax(ii) )
  end
  title(' $$ n^{(2)}(r_1,r_2) $$; $$ r_1 $$ at peak');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ n (r_0) n( r_{\perp} ) $$');
  axis square
  hl= legend(legCell );
end
%% plot min 
if plotMin
  figure()
  % c
  % min
  subplot(1,3,1)
  hold on
  for ii = 1:numFiles
    plot( x, cAtMinVec(:,ii) .* cMin(ii) )
  end
  title(' $$ c^{(2)}(r_1,r_2) $$; $$ r_1 $$ at min');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ c (r_0) c ( r_{\perp} ) $$');
  axis square
  hl= legend(legCell );
  % p
  % min
  subplot(1,3,2)
  hold on
  for ii = 1:numFiles
    plot( x, pAtMinVec(:,ii) .* pMin(ii) )
  end
  title(' $$ p^{(2)}(r_1,r_2) $$; $$ r_1 $$ at min');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ p (r_0) p( r_{\perp} ) $$');
  axis square
  hl= legend(legCell );
  % n
  % min
  subplot(1,3,3)
  hold on
  for ii = 1:numFiles
    plot( x, nAtMinVec(:,ii) .* nMin(ii) )
  end
  title(' $$ n^{(2)}(r_1,r_2) $$; $$ r_1 $$ at min');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ n (r_0) n( r_{\perp} ) $$');
  axis square
  hl= legend(legCell );
end
%% plot average 
if plotAverage
  figure()
  % c
  subplot(1,3,1)
  hold on
  for ii = 1:numFiles
    plot( x, c2bodyVec(:,ii) )
  end
  axis square
  title(' $$ c^{(2)}(r_1,r_2) $$; $$ r_1 $$ averaged');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ c (r_0) c( r_{\perp} ) $$');
  hl= legend(legCell );
  % p
  subplot(1,3,2)
  hold on
  for ii = 1:numFiles
    plot( x, p2bodyVec(:,ii) )
  end
  axis square
  title(' $$ p^{(2)}(r_1,r_2) $$; $$ r_1 $$ averaged');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ p (r_0) p( r_{\perp} ) $$');
  hl= legend(legCell );
  % n
  subplot(1,3,3)
  hold on
  for ii = 1:numFiles
    plot( x, n2bodyVec(:,ii) )
  end
  axis square
  title(' $$ n^{(2)}(r_1,r_2) $$; $$ r_1 $$ averaged');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ n (r_0) n( r_{\perp} ) $$');
  hl= legend(legCell );
end
