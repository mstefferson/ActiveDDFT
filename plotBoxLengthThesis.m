%% plotBoxSizeThesis.m
dirname = '/Users/mike/Projects/ActiveDDFT/analyzedfiles/BoxSize';
my_dirs = dir([dirname '/Hr*']);
saveName = 'box_size_compare';
saveMe = 0;
%%
% titles
fig = figure();
fig.WindowStyle = 'docked';
fig.Position = [472 -206 1036 646];
myTitle = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'};
% set some figur parameters
numRow = 2;
numCol = 2;
numT = length(my_dirs);
for ii = 1:numT
  % load params and opObj from current WD
  dir_temp = [my_dirs(ii).folder '/' my_dirs(ii).name];
  fileNameTemp = [dir_temp '/'];
  fileName = dir(fileNameTemp);
  paramsMat = dir( [ fileNameTemp '/params_*' ] );
  opMat = dir( [ fileNameTemp '/op_*' ] );
  paramLoad = load( [ fileNameTemp '/' paramsMat.name ] );
  fileObj = matfile( [ fileNameTemp   '/' opMat.name ] );
  % scale c
  b = 1 / pi;
  fprintf('Scaling C by b = %f\n', b );
  C = fileObj.C_rec(:,:,:) * b;
  % set-up figure
  % set-up grid things like ticks
  n2 = paramLoad.systemObj.n2;
  l2 = paramLoad.systemObj.l2;
  x2 =  l2/n2 * [-n2/2:n2/2-1];
  % Concentration
  axh1 = subplot(numRow,numCol,ii); % Save the handle of the subplot
  axis square
  cLim = [min(C(:)) max(C(:))];
  %fixAxis( axh1, C, xTick, yTick, xLim, yLim, cLim, myTitle{plotRowId+1} )
  % shift things
  [circAmount, shiftDim] = findShift( C, paramLoad.systemObj );
  C = shiftData( C, circAmount, shiftDim );
  nt = size(C, 3);
  C_slice = reshape( C(1,:,:), [length(x2), nt] );
  %% Concentration
  ind2plot = [1 2 5 nt];
  nPlot = length(ind2plot);
  myColors = viridis(nPlot);
  my_leg = cell(nPlot,1);
  hold on
  for nn = 1:nPlot
    plotid = ind2plot(nn);
    p = plot(x2, C_slice(:,plotid));
    p.Color = myColors(nn,:);
    my_leg{nn} = ['$$ t = $$' num2str(paramLoad.timeObj.t_rec * (plotid-1))];
  end
  axh1.XLim = [x2(1) x2(end)];
  xlabel('Position $$x$$')
  ylabel('$$C(x,t)$$')
  title(['l = ' num2str(paramLoad.systemObj.l2)])
  if ii == numT
    hl = legend(my_leg);
    hl.Interpreter = 'latex';
  end
end
if saveMe
  saveas(gcf, saveName,'png')
end
 
%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%

function [circAmount, shiftDim] = findShift( dataMap, systemObj )
% shift things
[~, rows, cols, ~, ~, ~] = ...
  sliceInhomoVar( systemObj, dataMap );
if length(rows) > length(cols)
  shiftDim =  1;
  [~,maxInd] = max( dataMap(:,1) );
  circAmount = mod( systemObj.n1/2 - maxInd, systemObj.n1 );
else
  shiftDim = 2;
  [~,maxInd] = max( dataMap(1,:) );
  circAmount = mod( systemObj.n2/2 - maxInd, systemObj.n2 );
end
end
function data2plot = shiftData( data2plot, circAmount, shiftDim )
% shift it to the center
data2plot = circshift( data2plot, circAmount, shiftDim );
end

