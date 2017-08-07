function ampPlotterFT(FTmat2plot, FTind2plot, TimeRec, ...
  kx0, ky0, km0, maxTime)

[totMode, totk] = size(FTind2plot);
figure()
plotColumns = 4;
plotRows = totMode / plotColumns;

for i = 1:totMode
  subplot(plotRows,plotColumns,i)
  ft2plot = FTmat2plot(i,:);
  realFt2plot = abs( real( ft2plot ) );
  imagFt2plot = abs( imag( ft2plot ) );
  [ax, ~, ~] = plotyy( TimeRec, realFt2plot, ...
    TimeRec, imagFt2plot );
  hold all
  scatter( TimeRec(end) , 0  );
  ax(1).XLim = [0 maxTime];  ax(2).XLim = [0 maxTime];
  fixYLimits( realFt2plot, ax(1) );
  fixYLimits( imagFt2plot, ax(2) );
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
  
  if(i == 1 || i == 5 ); ylabel(ax(1),' real Amp'); end
  if(i == 4 || i == 8 ); ylabel(ax(2),' imag Amp'); end
end
% subroutines
% fix limits
  function fixYLimits( v, ax )
    scaleEps = 0.01;
    meanVal = mean(v);
    if abs( std( v ) ./ mean(v) ) < scaleEps
      minLim = min(v) - meanVal * scaleEps;
      maxLim = max(v) + meanVal * scaleEps;
      ax.YLim = [ minLim maxLim ];
    end
  end
end




