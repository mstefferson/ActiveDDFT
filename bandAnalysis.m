function [output] = bandAnalysis( dirlist );

% Allocate for what we want
numFiles = length( dirlist );
% params and such
output.p1name = 'Ly';
output.p1 = zeros( 1, numFiles );
output.p2name = 'vD';
output.p2 = zeros( 1, numFiles );
output.p3name = 'bc';
output.p3 = zeros( 1, numFiles );
output.steady = zeros( 1, numFiles );
% c
output.cMax = zeros( 1, numFiles ); % Max C
output.cMin = zeros( 1, numFiles ); % Min C
output.cAve = zeros( 1, numFiles ); % Ave C
output.cDiff = zeros( 1, numFiles ); % C Max - C Min scaled by input bc
output.cFWHM = zeros( 1, numFiles ); % C full width half max
% p
output.pMax = zeros( 1, numFiles ); % Max P
output.pMin = zeros( 1, numFiles ); % Min P
output.pAve = zeros( 1, numFiles ); % Ave P
output.pDiff = zeros( 1, numFiles ); % P Max - P Min spaled by input bp
output.pFWHM = zeros( 1, numFiles ); % P full width half max
%n
output.nMax = zeros( 1, numFiles ); % Max N
output.nMin = zeros( 1, numFiles ); % Min N
output.nAve = zeros( 1, numFiles ); % Ave N
output.nDiff = zeros( 1, numFiles ); % N Max - N Min snaled by input bn
output.nFWHM = zeros( 1, numFiles ); % N full width half max

% find direction (x or y) of max variance
[maxVar, rows, cols, posVarLab, Lvar] = bandAnalysis( systemObj, x, y ) 

% loop over files
for ii = 1 : numFiles
  dirname = dirlist( ii );
  % cd into dir and load params and opObj
  cd( dirname );
  paramsMat = dir( 'param*' );
  opMat = dir( 'op*' );
  load(paramsMat.name);
  fileObj = matfile( opMat.name ); 
  C = fileObj.C_rec(:,:,end);
  P = fileObj.PO_rec(:,:,end);
  N = fileObj.NO_rec(:,:,end);
  % record parameters and steady state
  output.p1(ii) = systemObj.Ly;
  output.p2(ii) = systemObj.vD;
  output.p3(ii) = systemObj.bc;
  output.steady(ii) = denRecObj.SteadyState;
  % run analysis for C,N,P. Be carefule about P
  [cStats, pStats, nStats] = bandStatsCPNwrap( C(rows,cols), P(rows,cols), N(rows,cols) );
  % put them in a vector
  %c
  output.cMax(ii) = cStats.maxV;
  output.cMin(ii) = cStats.minV; 
  output.cAve(ii) = cStats.aveV; 
  output.cDiff(ii) = cStats.vdiff; 
  output.cFWHM(ii) = cStats.fwhm; 
  %p
  output.pMax(ii) = pStats.maxV;
  output.pMin(ii) = pStats.minV; 
  output.pAve(ii) = pStats.aveV; 
  output.pDiff(ii) = pStats.vdiff; 
  output.pFWHM(ii) = pStats.fwhm; 
  %n
  output.nMax(ii) = nStats.maxV;
  output.nMin(ii) = nStats.minV; 
  output.nAve(ii) = nStats.aveV; 
  output.nDiff(ii) = nStats.vdiff; 
  output.nFWHM(ii) = nStats.fwhm; 

end
