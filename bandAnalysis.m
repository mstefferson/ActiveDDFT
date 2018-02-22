function [output] = bandAnalysis( dirpath )
% get all the dirs
dirlist = dir([dirpath '/Hr*']);
% Allocate for what we want
numFiles = length( dirlist );
% params and such
fdVec = zeros( 1, numFiles );
cVec = zeros( 1, numFiles );
steady = zeros( 1, numFiles );
% c
cMax = zeros( 1, numFiles ); % Max C
cMin = zeros( 1, numFiles ); % Min C
cAve = zeros( 1, numFiles ); % Ave C
cDiff = zeros( 1, numFiles ); % C Max - C Min scaled by input bc
cFDWHD = zeros( 1, numFiles ); % C full width half (max - min)
% p
pMax = zeros( 1, numFiles ); % Max P
pMin = zeros( 1, numFiles ); % Min P
pAve = zeros( 1, numFiles ); % Ave P
pDiff = zeros( 1, numFiles ); % P Max - P Min spaled by input bp
pFDWHD = zeros( 1, numFiles ); % P full width half (max - min)
pP2P = zeros( 1, numFiles ); % P peak to peak distance
%n
nMax = zeros( 1, numFiles ); % Max N
nMin = zeros( 1, numFiles ); % Min N
nAve = zeros( 1, numFiles ); % Ave N
nDiff = zeros( 1, numFiles ); % N Max - N Min snaled by input bn
nFDWHD = zeros( 1, numFiles ); % N full width half (max - min)

% loop over files
for ii = 1 : numFiles
  dirname = dirlist( ii ).name;
  % load params and opObj from current WD
  fullpath = [ dirpath '/' dirname ];
  paramsMat = dir( [ fullpath '/params_*' ] );
  opMat = dir( [ fullpath '/op_*' ] );
  load( [ fullpath '/' paramsMat.name ] );
  fileObj = matfile( [ fullpath   '/' opMat.name ] ); 
  C = fileObj.C_rec(:,:,end);
  P = fileObj.POP_rec(:,:,end);
  N = fileObj.NOP_rec(:,:,end);
  % find direction (x or y) of max variance
  [~, rows, cols, ~, NposVar, Lvar] = ...
    sliceInhomoVar( systemObj, C );
  % grad parameters
  cVec(ii) = systemObj.bc;
  fdVec(ii) = particleObj.fD;
  steady(ii) = denRecObj.SteadyState;
  % run analysis for C,N,P. Be carefule about P
  [cStats, pStats, nStats] = ...
    bandStatsCPNwrap( C(rows,cols), P(rows,cols), N(rows,cols), Lvar, NposVar );
  % put them in a vector
  cMax(ii) = cStats.maxV;
  cMin(ii) = cStats.minV; 
  cAve(ii) = cStats.aveV; 
  cDiff(ii) = cStats.vdiff; 
  cFDWHD(ii) = cStats.fwhd;    
  %p
  pMax(ii) = pStats.maxV;
  pMin(ii) = pStats.minV; 
  pAve(ii) = pStats.aveV; 
  pDiff(ii) = pStats.vdiff; 
  pFDWHD(ii) = pStats.fwhd; 
  pP2P(ii) = pStats.p2p; 
  %n
  nMax(ii) = nStats.maxV;
  nMin(ii) = nStats.minV; 
  nAve(ii) = nStats.aveV; 
  nDiff(ii) = nStats.vdiff; 
  nFDWHD(ii) = nStats.fwhd; 
end % loop over files
% sort the parameters in case they are out of order
[fdVec, sortInd] =  sort( fdVec );
cVec = cVec(sortInd);
% get unique
num4round = 10000;
cUnique = unique( round( cVec * num4round ) / num4round );
fdUnique = unique( round( fdVec * num4round ) / num4round );
numC = length( cUnique );
numFd = length( fdUnique );
% Reshape 
reshapeVec = [numC, numFd];
steady = reshape( steady(sortInd), reshapeVec ).'; 
% c
cPeak = reshape( cMax(sortInd) ./ cVec(sortInd), reshapeVec ).'; 
cMax = reshape( cMax(sortInd), reshapeVec ).';  
cMin = reshape( cMin(sortInd), reshapeVec ).';  
cAve = reshape( cAve(sortInd), reshapeVec ).';  
cDiff = reshape( cDiff(sortInd), reshapeVec ).';  
cFDWHD = reshape( cFDWHD(sortInd), reshapeVec ).';  
% p
pPeak = reshape( pMax(sortInd), reshapeVec ).'; 
pMax = reshape( pMax(sortInd), reshapeVec ).';  
pMin = reshape( pMin(sortInd), reshapeVec ).';  
pAve = reshape( pAve(sortInd), reshapeVec ).';  
pDiff = reshape( pDiff(sortInd), reshapeVec ).';  
pFDWHD = reshape( pFDWHD(sortInd), reshapeVec ).';  
% n
nPeak = reshape( nMax(sortInd), reshapeVec ).'; 
nMax = reshape( nMax(sortInd), reshapeVec ).';  
nMin = reshape( nMin(sortInd), reshapeVec ).';  
nAve = reshape( nAve(sortInd), reshapeVec ).';  
nDiff = reshape( nDiff(sortInd), reshapeVec ).';  
nFDWHD = reshape( nFDWHD(sortInd), reshapeVec ).';  
% summary
output.fd = fdUnique;
output.c = cUnique;
output.steady = steady;
% c
output.cPeak = cPeak;
output.cMax = cMax; 
output.cMin = cMin;
output.cAve = cAve;
output.cDiff = cDiff;
output.cFDWHD = cFDWHD;
% p
output.pPeak = pPeak;
output.pMax = pMax; 
output.pMin = pMin;
output.pAve = pAve;
output.pDiff = pDiff;
output.pFDWHD = pFDWHD;
% n
output.nPeak = nPeak;
output.nMax = nMax; 
output.nMin = nMin;
output.nAve = nAve;
output.nDiff = nDiff;
output.nFDWHD = nFDWHD;
end
