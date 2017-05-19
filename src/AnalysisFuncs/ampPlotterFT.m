function ampPlotterFT(FTmat2plot, FTind2plot, TimeRec, kx0, ky0, km0, maxTime)

[totMode, totk] = size(FTind2plot);
figure()
plotColumns = 4;
plotRows = totMode / plotColumns;

for i = 1:totMode
  subplot(plotRows,plotColumns,i)
  [Ax, ~, ~] = plotyy( TimeRec, abs( real( FTmat2plot(i,:) ) ), ...
    TimeRec, abs( imag( FTmat2plot(i,:) ) ) );
  hold all
  scatter( [TimeRec(end) ], ...
    [0 ] );
  Ax(1).XLim = [0 maxTime];  Ax(2).XLim = [0 maxTime];
  if totk == 3
    titstr = sprintf('(%d, %d, %d)', ...
      FTind2plot(i,1) - kx0, FTind2plot(i,2) - ky0, FTind2plot(i,3) - km0 );
  elseif totk == 2
    titstr = sprintf('(%d, %d)', ...
      FTind2plot(i,1) - kx0, FTind2plot(i,2) - ky0 );
  else
    titstr = 'no title';
  end
  title(titstr)
  xlabel('Time')
  
  if(i == 1 || i == 5 ); ylabel(Ax(1),' real Amp'); end
  if(i == 4 || i == 8 ); ylabel(Ax(2),' imag Amp'); end
  
end

end
