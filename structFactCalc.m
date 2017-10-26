function structFactCalc( rhoK )

% \langle \frac{1}{N} \rho( k ) \rho( -k ) \rangle
% where \rho( k ) is the FT of \rho( r ) = \sum_i \delta ( r - r_i )


% naive calculation
% flip k -> -k
rhoMinusK = flip( flip( flip( rhoK, 1 ), 2 ), 3 );
% shift k = 0 back to center
rhoMinusK = circshift( rhoMinusK, [1 1 1] );
sk_naive = rhoK .* rhoMinusK;

n = size( rhoK, 3 );
nCenter = n / 2 + 1;
dPlot = 10;
plotInds = nCenter-dPlot:nCenter+dPlot;
% plot subplot
figure()
ind = nCenter;
subplot(1,2,1)
pcolor( real( sk_naive(plotInds,plotInds,ind) ) );
colorbar;
subplot(1,2,2)
pcolor( imag( sk_naive(plotInds,plotInds,ind) ) );
colorbar;
title( ['ind = ' num2str( ind - nCenter, '%d' ) ] )
% inds for slices
ind1 = nCenter;
ind2 = plotInds;
scanInds = 0:5;
% plot slices
figure()
ax1 = gca;
hold on
title('real')
for ii = scanInds
  indTemp = nCenter+ii;
  plot( ax1, plotInds, real( sk_naive(ind1,ind2,indTemp) ) );
end
% plot slices log
figure()
ax1sl = gca;
hold on
title('real')
ax1sl.YScale = 'log';
for ii = scanInds
  indTemp = nCenter+ii;
  plot( ax1sl, plotInds, real( sk_naive(ind1,ind2,indTemp) ) );
end