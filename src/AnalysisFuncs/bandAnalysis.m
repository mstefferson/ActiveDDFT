function [output,tSum] = bandAnalysis( dirpath )
fprintf('Running %s\n',mfilename);
% add paths
addpath( genpath( pwd ) );
% get all the dirs
dirlist = dir([dirpath '/Hr*']);
% Allocate for what we want
numFiles = length( dirlist );
% scale c
b = 1 / pi;
fprintf('Scaling c by b = %f\n', b ); 
% params and such
fdVec = zeros( 1, numFiles );
cVec = zeros( 1, numFiles );
steady = zeros( 1, numFiles );
% c
cSlice = cell( 1, numFiles ); % store a slice
cMax = zeros( 1, numFiles ); % Max C
cMin = zeros( 1, numFiles ); % Min C
cAve = zeros( 1, numFiles ); % Ave C
cDiff = zeros( 1, numFiles ); % C Max - C Min scaled by input bc
cFDWHD = zeros( 1, numFiles ); % C full width half (max - min)
% p
pSlice = cell( 1, numFiles ); % store a slice
pMax = zeros( 1, numFiles ); % Max P
pMin = zeros( 1, numFiles ); % Min P
pAve = zeros( 1, numFiles ); % Ave P
pDiff = zeros( 1, numFiles ); % P Max - P Min spaled by input bp
pFDWHD = zeros( 1, numFiles ); % P full width half (max - min)
pP2P = zeros( 1, numFiles ); % P peak to peak distance
%n
nSlice = cell( 1, numFiles ); % store a slice
nMax = zeros( 1, numFiles ); % Max N
nMin = zeros( 1, numFiles ); % Min N
nAve = zeros( 1, numFiles ); % Ave N
nDiff = zeros( 1, numFiles ); % N Max - N Min snaled by input bn
nFDWHD = zeros( 1, numFiles ); % N full width half (max - min)
% Starting loop over files
fprintf('Starting loop over %d files\n', numFiles );
for ii = 1 : numFiles
  dirname = dirlist( ii ).name;
  % load params and opObj from current WD
  fullpath = [ dirpath '/' dirname ];
  paramsMat = dir( [ fullpath '/params_*' ] );
  opMat = dir( [ fullpath '/op_*' ] );
  paramLoad = load( [ fullpath '/' paramsMat.name ] );
  fileObj = matfile( [ fullpath   '/' opMat.name ] );
  C = fileObj.C_rec(:,:,end);
  P = fileObj.POP_rec(:,:,end);
  N = fileObj.NOP_rec(:,:,end);
  % find direction (x or y) of max variance
  [~, rows, cols, ~, NposVar, Lvar] = ...
    sliceInhomoVar( paramLoad.systemObj, C );
  % grad parameters
  cVec(ii) = paramLoad.systemObj.bc;
  fdVec(ii) = paramLoad.particleObj.fD;
  steady(ii) = paramLoad .denRecObj.SteadyState;
  % run analysis for C,N,P. Be carefule about P
  C =  C(rows,cols) * b;
  [cStats, pStats, nStats] = ...
    bandStatsCPNwrap(  C, P(rows,cols), N(rows,cols), Lvar, NposVar );
  % put them in a vector
  cSlice{ii} = C;
  cMax(ii) = cStats.maxV;
  cMin(ii) = cStats.minV;
  cAve(ii) = cStats.aveV;
  cDiff(ii) = cStats.vdiff;
  cFDWHD(ii) = cStats.fwhd;
  %p
  pSlice{ii} = P(rows,cols);
  pMax(ii) = pStats.maxV;
  pMin(ii) = pStats.minV;
  pAve(ii) = pStats.aveV;
  pDiff(ii) = pStats.vdiff;
  pFDWHD(ii) = pStats.fwhd;
  pP2P(ii) = pStats.p2p;
  %n
  nSlice{ii} = N(rows,cols);
  nMax(ii) = nStats.maxV;
  nMin(ii) = nStats.minV;
  nAve(ii) = nStats.aveV;
  nDiff(ii) = nStats.vdiff;
  nFDWHD(ii) = nStats.fwhd;
end % loop over files
fprintf('Finished loop over files. Storing data\n');
% sort the parameters in case they are out of order
[fdVec, sortInd] =  sort( fdVec );
cVec = cVec(sortInd);
% Reshape
steady = steady(sortInd);
% c
[cSlice, cPeak, cMax, cMin, cAve, cDiff, cFDWHD ] = ...
  sortData( cSlice, cMax, cMin, cAve, cDiff, cFDWHD, sortInd );
% scale c
cPeak = cPeak ./ cVec;
% p
[pSlice, pPeak, pMax, pMin, pAve, pDiff, pFDWHD ] = ...
  sortData( pSlice, pMax, pMin, pAve, pDiff, pFDWHD, sortInd );
% n
[nSlice, nPeak, nMax, nMin, nAve, nDiff, nFDWHD ] = ...
  sortData( nSlice, nMax, nMin, nAve, nDiff, nFDWHD, sortInd );
% summary
output.fd = fdVec;
output.c = cVec;
output.steady = steady;
output.sortInd = sortInd;
% c
output.cSlice = cSlice;
output.cPeak = cPeak;
output.cMax = cMax;
output.cMin = cMin;
output.cAve = cAve;
output.cDiff = cDiff;
output.cFDWHD = cFDWHD;
% p
output.pSlice = pSlice;
output.pPeak = pPeak;
output.pMax = pMax;
output.pMin = pMin;
output.pAve = pAve;
output.pDiff = pDiff;
output.pFDWHD = pFDWHD;
% n
output.nSlice = nSlice;
output.nPeak = nPeak;
output.nMax = nMax;
output.nMin = nMin;
output.nAve = nAve;
output.nDiff = nDiff;
output.nFDWHD = nFDWHD;
% build a table

%%
tSum = table( output.fd', output.c', output.steady', ...
  output.cPeak', output.cSlice', output.cMax', output.cMin',...
  output.cAve', output.cDiff', output.cFDWHD',...
  output.pPeak', output.pSlice', output.pMax',...
  output.pMin', output.pAve', output.pDiff', output.pFDWHD',...
  output.nPeak', output.nSlice', output.nMax', ...
  output.nMin', output.nAve', output.nDiff', output.nFDWHD');
%%

tSum.Properties.VariableNames = {'fd', 'c', 'steady', ...
  'cPeak', 'cSlice', 'cMax', 'cMin', 'cAve', 'cDiff', 'cFDWHD',...
  'pPeak', 'pSlice', 'pMax', 'pMin', 'pAve', 'pDiff', 'pFDWHD',...
  'nPeak', 'nSlice', 'nMax', 'nMin', 'nAve', 'nDiff', 'nFDWHD'};

  function [tempSlice, tempPeak, tempMax, tempMin, tempAve, tempDiff, tempFDWHD ] = ...
      sortData( tempSlice, tempMax, tempMin, tempAve, tempDiff, tempFDWHD, sortInd )
    tempSlice = tempSlice(sortInd);
    tempPeak = tempMax(sortInd);
    tempMax = tempMax(sortInd);
    tempMin = tempMin(sortInd);
    tempAve = tempAve(sortInd);
    tempDiff = tempDiff(sortInd);
    tempFDWHD = tempFDWHD(sortInd);
  end
end
