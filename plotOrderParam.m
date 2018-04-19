function plotOrderParam(bandTable, xStr, yStr, xLabel, yLabel)
% grab data
xData = bandTable(:, {xStr}).Variables; 
yData = bandTable(:, {yStr}).Variables; 
xData = xData( bandTable.fd == 0);
yData = yData( bandTable.fd == 0);
% plot it
plot(xData,yData)
ax = gca;
ax.XLim = [min(xData) max(xData)];
ax.YLim = [0 1];
xlabel(xLabel)
ylabel(yLabel)
axis square
grid on
