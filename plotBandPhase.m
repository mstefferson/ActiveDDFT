function plotBandPhase( bandTable, plotTheory, plotScaled )
% plot phase
cIN = 1.5;
% set-up figure
fig = figure();
fig.WindowStyle = 'normal';
fig.Position = [360 278 560 420];
% Scaled Pe and concentration
if plotScaled
  xLabel = '$$ c^* / c_{IN} $$';
  yLabel = '$$ \sqrt{ \frac{ Pe } { 6 } }$$';
  % functional guess
  c = linspace( 1.5, 1.4*1.5 );
  cS = c / cIN;
  cTheory2plot = cS;
  fUnstable = sqrt(2 * ( cS.^2 - 1 ) ) .* ( cS + 1 ) ./ ...
    ( cS .^ 2 + cS - 1 ); % actve nematic
  %fUnstable2 = sqrt( ( cS.^2 - 1 ) ); % self-regulation
  fPlot = sqrt( bandTable.fd / 6) ;
  cPlot = bandTable.c ./ cIN;
  plotPhaseDiagram( cPlot, fPlot, bandTable.cPeak, cTheory2plot, fUnstable,...
    plotTheory, xLabel, yLabel)
  % Unscaled Pe and concentration
else
  xLabel = 'Concentration $$  C^* $$';
  yLabel = 'P\''eclet Number $$ Pe $$';
  c = linspace( 1.5, 1.4*1.5 );
  cS = c/cIN;
  cTheory2plot = c;
  fUnstable = 12 * ( cS.^2 - 1 )  .* ( cS + 1 ) .^2 ./ ...
    ( cS .^ 2 + cS - 1 ) .^ 2; % actve nematic
  %fUnstable2 = sqrt( ( cS.^2 - 1 ) ); % self-regulation
  fPlot = bandTable.fd;
  cPlot = bandTable.c;
  plotPhaseDiagram( cPlot, fPlot, bandTable.cPeak, cTheory2plot, fUnstable,...
    plotTheory, xLabel, yLabel)
end
%% functions
  function plotPhaseDiagram( cPlot, fPlot, cPeak, cTheory2plot, fUnstable,...
      plotTheory, xLabel, yLabel)
    % theory  lines
    circleSize = 75;
    lineIN = ':';
    lineTheory = '-';
    lineColorIN = [0 0 0];
    lineColorTheory = [0 0 0];
    threshold = 1.01;
    % phase diagram
    ax = gca;
    axis square
    hold on
    % plot failed
    failInd = cPeak < threshold;
    p = scatter( cPlot(failInd), fPlot(failInd), 10, cPeak(failInd) );
    p.Marker = 'o';
    p.SizeData = circleSize;
    % plot success
    p = scatter( cPlot(~failInd), fPlot(~failInd), 10, cPeak(~failInd), 'filled' );
    p.Marker = 'o';
    p.SizeData =  circleSize;
    % plot IN transtion
    bigSlope = 1000000;
    plot( cTheory2plot, bigSlope * (cTheory2plot - cTheory2plot(1)), ...
      'Color', lineColorIN, 'LineStyle', lineIN )
    if plotTheory
      plot( cTheory2plot, fUnstable, ...
        'Color', lineColorTheory, 'LineStyle', lineTheory )
      % plot IN transtion
      bigSlope = 10000;
      plot( cTheory2plot, bigSlope * (cTheory2plot - cTheory2plot(1)), ...
        'Color', lineColorIN, 'LineStyle', lineIN )
      % legend
      legCell = {'Homo.', 'Band', 'IN trans.', 'Theory',};
    else
      % legend
      legCell = {'Homo.', 'Band', 'IN trans.'};
      cTheory2plot = cPlot;
    end
    % axis
    % set-up axis properties
    h = colorbar;
    h.Label.Interpreter = 'latex';
    h.Label.String = '$$ \frac{ c_{peak} }{ c^* } $$';
    h.TickLabelInterpreter = 'latex';
    h.Label.Position = [4.5 1.75 0 ];
    h.Label.Rotation = 0;
    hold off
    xlabel(ax, xLabel);
    ylabel(ax, yLabel);
    ax.XLim = [min( [cTheory2plot(:); cPlot(:)] ), ...
      max([cTheory2plot(:); cPlot(:)]) ];
    ax.YLim = [min(fPlot) max(fPlot) ];
    ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.TickLabelInterpreter = 'latex';
    box on
    % legend
    [hl] = legend(legCell, ...
      'location', 'best');
    hl.Interpreter = 'latex';
    hl.Position = [0.5570 0.6946 0.2003 0.2195];
  end
end
