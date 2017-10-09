function plotPdistOps( path2dir )
% get dir
dir2analyze = dir( [path2dir '/Hr_*'] );
numDirs = length(dir2analyze);
if numDirs == 0; fprintf('No dirs found\n'); end
% allocate
bcVec = zeros(1,numDirs);
vdVec = zeros(1,numDirs);
% set-up figure()
fig1 = figure();
% c rows/columns
numRow = 2;
numCol = 4;
ax1u = subplot(numRow,numCol,1);
hold(ax1u, 'on')
xlabel('$$ x $$');
ylabel('$$ C (x,0) $$');
ax1l = subplot(numRow,numCol,5);
hold(ax1l, 'on')
xlabel('$$ y $$');
ylabel('$$ C (0,y) $$');
% p rows/columns
ax2u = subplot(numRow,numCol,2);
hold(ax2u, 'on')
xlabel('$$ x $$');
ylabel('$$ P (x,0) $$');
ax2l = subplot(numRow,numCol,6);
hold(ax2l, 'on')
xlabel('$$ y $$');
ylabel('$$ P (0,y) $$');
% n rows/columns
ax3u = subplot(numRow,numCol,3);
hold(ax3u, 'on')
xlabel('$$ x $$');
ylabel('$$ N (x,0) $$');
ax3l = subplot(numRow,numCol,7);
hold(ax3l, 'on')
xlabel('$$ y $$');
ylabel('$$ N (0,y) $$');
% pdist rows/columns
ax4u = subplot(numRow,numCol,4);
hold(ax4u, 'on')
xlabel('$$ x $$');
ylabel('$$ p_0 (x,0) $$');
ax4l = subplot(numRow,numCol,8);
hold(ax4l, 'on')
xlabel('$$ y $$');
ylabel('$$ p_0 (0,y) $$');
% loop over files
for ii = 1:numDirs
  % get paths and load
  pathTemp = dir2analyze(ii).folder;
  dirTemp = dir2analyze(ii).name;
  fprintf('Plotting for %s\n', dirTemp);
  wd = [pathTemp  '/' dirTemp ];
  opStr = [wd '/op_' dirTemp '.mat'];
  paramStr = [wd '/params_' dirTemp '.mat'];
  load( [wd '/pDist_' dirTemp '.mat'] );
  load( paramStr );
  load( opStr );
  % grab bc and vd vectors
  bcVec(ii) = systemObj.bc;
  vdVec(ii) = particleObj.vD;
  % get center ind
  center1 = systemObj.n1 / 2 + 1;
  center2 = systemObj.n2 / 2 + 1;
  x = systemObj.l1/systemObj.n1*(-systemObj.n1/2:systemObj.n1/2-1);
  y = systemObj.l2/systemObj.n2*(-systemObj.n2/2:systemObj.n2/2-1);
  % grab last value
  c = C_rec(:,:,end) ./ pi;
  pop = POP_rec(:,:,end);
  nop = NOP_rec(:,:,end);
  % get max ind
  [~, maxInd] = max( c(:) );
  [max1, max2] = ind2sub( [systemObj.n1 systemObj.n2], maxInd);
  % plot it
  % x
  plot(ax1u, x, circshift( c(:,max2), center1 - max1 ) );
  plot(ax2u, x, circshift( pop(:,max2), center1 - max1 ) );
  plot(ax3u, x, circshift( nop(:,max2), center1 - max1 ) );
  plot(ax4u, x, pDistDim1.pDist0Center(:,1) )
  % y
  plot(ax1l, y, circshift( c(max1,:), center2 - max2 ) );
  plot(ax2l, y, circshift( pop(max1,:), center2 - max2 ) );
  plot(ax3l, y, circshift( nop(max1,:), center2 - max2 ) );
  plot(ax4l, x, pDistDim2.pDist0Center(1,:) )
end
% build legend
bcVecUn = unique( bcVec );
vdVecUn = unique( vdVec );
if length(bcVecUn) == 1
  legendCell = cell(1, length( vdVec ) );
  for ii = 1:length( vdVec )
    legendCell{ii} = [' $$ v_D =  ' num2str( vdVec(ii), '%-.1f') ' $$' ];
  end
elseif length(vdVecUn) == 1
    legendCell = cell(1, length( bcVec ) );
  for ii = 1:length( bcVec )
    legendCell{ii} = [' $$ c^* =  ' num2str( bcVec(ii), '%-.1f') ' $$' ];
  end
else
  legendCell = cellstr([ num2str(bcVec', '$$ c^* = $$ %-.2f ')  ...
    num2str(vdVec', '$$ v_D = $$ %-.1f') ]);
end
leg = legend( ax4l, legendCell );
leg.Position = [0.8683 0.6466 0.0779 0.3081];
leg.Interpreter = 'latex';