function plotPdistOps( path2dir )
% get dir
dir2analyze = dir( [path2dir '/Hr_*'] );
numDirs = length(dir2analyze);
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
for ii = 1:length(dir2analyze)
  dirTemp = dir2analyze(ii).name;
  fprintf('Plotting for %s\n', dirTemp);
  wd = [path2dir  '/' dirTemp ];
  pDistSearch = dir( [wd '/pDist*'] );
  opStr = [wd '/op_' dirTemp '.mat'];
  paramStr = [wd '/params_' dirTemp '.mat'];
  load( [pDistSearch.folder '/' pDistSearch.name] );
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
  % plot it
  % x
  plot(ax1u, x, c(:,center2) )
  plot(ax2u, x, pop(:,center2) )
  plot(ax3u, x, nop(:,center2) )
  plot(ax4u, x, pDistDim1.pDist0(:,1) )
  % y
  plot(ax1l, y, c(center1,:) )
  plot(ax2l, y, pop(center1,:) )
  plot(ax3l, y, nop(center1,:) )
  plot(ax4u, x, pDistDim2.pDist0(1,:) )
end

% build legend
bcVec = unique( bcVec );
vdVec = unique( vdVec );
if length(bcVec) == 1
  legendCell = cellstr(num2str(vdVec', 'vD=%-.1f'));
elseif length(vdVec) == 1
  legendCell = cellstr(num2str(bcVec', 'bc=%-.2f'));
else
  legendCell = cellstr([ num2str(bcVec', 'bc=%-.2f ')   num2str(vdVec', 'vD=%-.1f') ]);
end
leg = legend( ax4l, legendCell );
leg.Position = [0.8615    0.8346    0.0914    0.1320];

