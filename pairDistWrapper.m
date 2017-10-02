% Wrapper for pairDistribution calculator
%
function pairDistWrapper()
% set flags
g1g2Flag = 0;
plotFlag = 0;
% add paths just in case
currentDir = pwd;
addpath( genpath( [currentDir '/src'] ) );
% grab them from a specific dir using 'path4pairdist.m'. Or, run with all 
% files in WD
if exist('path4pairdist.m','file')
  moveMe = 0;
  path2dir = path4pairdist;
  dir4move = path2dir;
else
  moveMe = 1;
  path2dir = pwd;
  dir4move = [date '_pDist' ];
  mkdir(  dir4move );
end
dir2analyze = dir( [path2dir '/Hr_*'] );
% save name stuff
saveNameRoot = 'pDist';
% loop over files

for ii = 1:length(dir2analyze)
  dirTemp = dir2analyze(ii).name;
  fprintf('Analyzing for %s\n', dirTemp);
  wd = [path2dir  '/' dirTemp ];
  rhoStr = [wd '/rhoFinal_' dirTemp '.mat'];
  paramStr = [wd '/params_' dirTemp '.mat'];
  load( paramStr );
  load( rhoStr );
  % calculate slice along dim 1
  saveName = [ saveNameRoot '1_' dirTemp '.mat' ];
  % calulate pair dist
  pDistDim1 = pairDistCalcRho( rho, systemObj.l1, systemObj.l2, ...
    particleObj.lMaj, 1:systemObj.n1, 1, g1g2Flag, plotFlag );
  % calulate pair dist
  pDistDim2 = pairDistCalcRho( rho, systemObj.l1, systemObj.l2, ...
    particleObj.lMaj, 1, 1:systemObj.n2, g1g2Flag, plotFlag );
  % calculate slice along dim 1
  saveName = [ saveNameRoot '/' dirTemp '.mat' ];
  save( saveName, 'pDistDim1', 'pDistDim2' );
  % move it
  if moveMe
    movefile( wd, dir4move )
  end
  movefile( [saveNameRoot '*'], [dir4move '/' dirTemp] )
end
