% Wrapper for pairDistribution calculator
%
function pairDistWrapper()
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
% loop over files

for ii = 1:length(dir2analyze)
  dirTemp = dir2analyze(ii).name;
  fprintf('Analyzing for %s\n', dirTemp);
  wd = [path2dir  '/' dirTemp ];
  rhoStr = [wd '/rhoFinal_' dirTemp '.mat'];
  paramStr = [wd '/params_' dirTemp '.mat'];
  load( paramStr );
  load( rhoStr );
  saveName = ['pDist_' dirTemp '.mat' ];
  % calulate pair dist
  pairDistCalcRho( rho, systemObj.l1, systemObj.l2, particleObj.lMaj, 0, saveName );
  % move it
  if moveMe
    movefile( wd, dir4move )
  end
  movefile(saveName, [dir4move '/' dirTemp] )
end
