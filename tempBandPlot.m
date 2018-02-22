if exist( 'bandTable', 'var') == 0
  dir2check = 'analyzedfiles/BandAnalysis';
  [bandSum, bandTable] = bandAnalysis(dir2check);
end
%% phase diagram
figure()
markSize = 20;
p = scatter( bandTable.c, bandTable.fd, 10, bandTable.cPeak );
p.Marker = '.';
p.SizeData = 1000;
colorbar
title('Band Phase diagram');
xlabel('$$c^*$$');
ylabel('$$f^d$$');
%% band slices
% get data
cWant = 1.55;
fd2Plot = bandTable.fd( bandTable.c == cWant );
[fd2Plot, indSort] = sort( fd2Plot );
cSlices2Plot = bandTable.cSlice( bandTable.c == cWant );
pSlices2Plot = bandTable.pSlice( bandTable.c == cWant );
nSlices2Plot = bandTable.nSlice( bandTable.c == cWant );
cSlices2Plot = cSlices2Plot(indSort);
pSlices2Plot = pSlices2Plot(indSort);
nSlices2Plot = nSlices2Plot(indSort);
n = length( cSlices2Plot{1} );
l = 10;
x = l / n * (-n/2:n/2-1);
% plot it
%%
figure()
ax1 = subplot(1,3,1);
axis square
xlabel('$$ x $$')
ylabel('$$ c/c^* $$')
hold on
ax2 = subplot(1,3,2);
axis square
xlabel('$$ x $$')
ylabel('$$ p $$')
hold on
ax3 = subplot(1,3,3);
axis square
xlabel('$$ x $$')
ylabel('$$ n $$')
hold on
%%
for ii = 1:length(fd2Plot)
  %%
  % get max for shifting
  [~, maxInd] = max(cSlices2Plot{ii});
  circAmount = mod( n/2 - maxInd, n );
  % grab data
  c2plot = cSlices2Plot{ii} ./ cWant;
  p2plot = pSlices2Plot{ii};
  n2plot = nSlices2Plot{ii};
  % make sure everything is row vector
  if iscolumn(c2plot); c2plot = c2plot'; end
  if iscolumn(p2plot); p2plot = p2plot'; end
  if iscolumn(n2plot); n2plot = n2plot'; end
  c2plot = circshift( c2plot, circAmount, 2 );
  p2plot = circshift( p2plot, circAmount, 2 ); 
  n2plot = circshift( n2plot, circAmount, 2 ); 
  plot(ax1, x, c2plot )
  plot(ax2, x, p2plot )
  plot(ax3, x, n2plot )
end
legendCell = cellstr(num2str(fd2Plot, 'fd=%-d'));
legend( legendCell, 'location', 'best' )
