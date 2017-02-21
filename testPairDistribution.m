%% All files
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
% c
cMaxVec = zeros( numFiles, numGrid );
cMax = zeros( numFiles, 1 );
cMinVec = zeros( numFiles, numGrid );
cMin = zeros( numFiles, 1 );
c2bodyVec = zeros( numFiles, numGrid );
% p
pMaxVec = zeros( numFiles, numGrid );
pMax = zeros( numFiles, 1 );
pMinVec = zeros( numFiles, numGrid );
pMin = zeros( numFiles, 1 );
p2bodyVec = zeros( numFiles, numGrid );
% n
nMaxVec = zeros( numFiles, numGrid );
nMax = zeros( numFiles, 1 );
nMinVec = zeros( numFiles, numGrid );
nMin = zeros( numFiles, 1 );
n2bodyAveVec = zeros( numFiles, numGrid );
% Loop
for ii = 1: numFiles
  filename = files{ii};
  path = ['./analyzedfiles/summarized/vDsweep/' filename ];
  load( [path '/op_' filename '.mat'] )
  load( [path '/params_' filename '.mat'] )
  % get position vector
  x = 0 : systemObj.Ly / systemObj.Ny : systemObj.Ly - systemObj.Ly / systemObj.Ny;
  % Get steady state distribution
  cFinal = C_rec(:,:,end) ./ pi;
  nFinal = NOP_rec(:,:,end);
  % get a slice that is averaged over non-varying direction
  cSlice = mean( cFinal, 1 );
  pSlice = mean( pFinal, 1 );
  nSlice = mean( nFinal, 1 );
  % find max
  [cMax, maxInd] = max( cSlice );
  [cMin] = min( cSlice );
  [nMax] = max( nSlice );
  [nMax] = max( nSlice );
  [nMin] = min( nSlice );
  [nMin] = min( nSlice );
  % circ and save
  % c
  cSliceMax = circshift( cSlice, -maxInd );
  cSliceMin = circshift( cSlice, -minInd );
  cMaxVec(ii,:) = cSliceMax;
  cMax(ii) = cMax;
  cMinVec(ii,:) = cSliceMin;
  cMin(ii) = cMin;
  % p
  pSliceMax = circshift( pSlice, -maxInd );
  pSliceMin = circshift( pSlice, -minInd );
  pMaxVec(ii,:) = pSliceMax;
  pMax(ii) = pMax;
  pMinVec(ii,:) = pSliceMin;
  pMin(ii) = pMin;
  % p
  nSliceMax = circshift( nSlice, -maxInd );
  nSliceMin = circshift( nSlice, -minInd );
  nMaxVec(ii,:) = nSliceMax;
  nMax(ii) = nMax;
  nMinVec(ii,:) = nSliceMin;
  nMin(ii) = nMin;
  % average
  for jj = 1:numGrid
    deltaInd = jj-1;
    inds = 1+deltaInd:deltaInd:numGrid+deltaInd;
    inds = mod( inds - 1, numGrid ) + 1;
    c2bodyVec(ii,jj) = mean( cMaxVec( 1:numGrid ) .* cMaxVec( inds ) ); 
    p2bodyVec(ii,jj) = mean( pMaxVec( 1:numGrid ) .* pMaxVec( inds ) ); 
    n2bodyVec(ii,jj) = mean( nMaxVec( 1:numGrid ) .* nMaxVec( inds ) ); 
  end
end
% plot it
for ii = 1: numFiles
  figure()
  % c
  subplot(2,2,1)
  imagesc( x,x, cFinal )
  axis square
  title(['C(x,y) vD = ' num2str(particleObj.vD)] );
  xlabel('x');
  ylabel('y');
  colorbar;
  % n
  subplot(2,2,2)
  imagesc( x,x, nFinal )
  axis square
  title(['N(x,y) vD = ' num2str(particleObj.vD)] );
  xlabel('x');
  ylabel('y');
  colorbar;
  % from peak c
  subplot(2,2,3)
  plot( x, cSliceMax .* cMax )
  title(' $$ \rho^{(2)}(r_1,r_2) $$; $$ r_1 $$ at peak');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ \rho (r_0) \rho( r_{\perp} ) $$');
  axis square
  % from min c
  subplot(2,2,4)
  plot( x, cSliceMin .* cMin )
  title(' $$ \rho^{(2)}(r_1,r_2) $$; $$ r_1 $$ at min');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ \rho (r_0) \rho( r_{\perp} ) $$');
  axis square
end
% plot it
for ii = 1:numFiles
  figure()
  % c
  subplot(1,3,1)
  plot( x, c2bodyVec )
  axis square
  title(' $$ c^{(2)}(r_1,r_2) $$; $$ r_1 $$ averaged');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ c (r_0) c( r_{\perp} ) $$');
  % p
  subplot(1,3,2)
  plot( x, p2bodyVec )
  title(' $$ p^{(2)}(r_1,r_2) $$; $$ r_1 $$ averaged');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ p(r_0) p( r_{\perp} ) $$');
  axis square
  % n
  subplot(1,3,3)
  plot( x, n2bodyVec )
  title(' $$ n^{(2)}(r_1,r_2) $$; $$ r_1 $$ averaged');
  xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
  ylabel('$$ n(r_0) n( r_{\perp} ) $$');
  axis square
end
%% plot them together
% max 
figure()
subplot(1,2,1)
hold on
for ii = 1:numFiles
  plot( x, cMaxVec(ii,:) .* cMax(ii) )
end
title(' $$ \rho^{(2)}(r_1,r_2) $$; $$ r_1 $$ at peak');
xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
ylabel('$$\rho (r_0) \rho( r_{\perp} )$$');
axis square
hl= legend(legCell );
% min
subplot(1,2,2)
hold on
for ii = 1:numFiles
  plot( x, cMinVec(ii,:) .* cMin(ii) )
end
title(' $$ \rho^{(2)}(r_1,r_2) $$; $$ r_1 $$ at min');
xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
ylabel('$$\rho (r_0) \rho( r_{\perp} )$$');
axis square
hl= legend(legCell );
% plot it
figure()
% c
subplot(1,3,1)
hold on
for ii = 1:numFiles
  figure()
  plot( x, c2bodyVec )
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
  figure()
  plot( x, p2bodyVec )
end
axis square
title(' $$ p^{(2)}(r_1,r_2) $$; $$ r_1 $$ averaged');
xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
ylabel('$$ p (r_0) p( r_{\perp} ) $$');
hl= legend(legCell );
% n
subplot(1,3,2)
hold on
for ii = 1:numFiles
  figure()
  plot( x, n2bodyVec )
end
axis square
title(' $$ n^{(2)}(r_1,r_2) $$; $$ r_1 $$ averaged');
xlabel('$$ ( r_{\perp}-r_0 ) / l_{rod} $$');
ylabel('$$ n (r_0) n( r_{\perp} ) $$');
hl= legend(legCell );
