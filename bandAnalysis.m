function [output] = bandAnalysis( dirpath, plotFlag )
% Set plot flag to 1 if not specified
if nargin == 1
  plotFlag = 1;
end

% get all the dirs
dirlist = dir([dirpath '/Hr*']);
% Allocate for what we want
numFiles = length( dirlist );
% params and such
output.p1name = 'l2';
output.p1 = zeros( 1, numFiles );
output.p2name = .fD';
output.p2 = zeros( 1, numFiles );
output.p3name = 'bc';
output.p3 = zeros( 1, numFiles );
output.steady = zeros( 1, numFiles );
% c
output.cMax = zeros( 1, numFiles ); % Max C
output.cMin = zeros( 1, numFiles ); % Min C
output.cAve = zeros( 1, numFiles ); % Ave C
output.cDiff = zeros( 1, numFiles ); % C Max - C Min scaled by input bc
output.cFDWHD = zeros( 1, numFiles ); % C full width half (max - min)
% p
output.pMax = zeros( 1, numFiles ); % Max P
output.pMin = zeros( 1, numFiles ); % Min P
output.pAve = zeros( 1, numFiles ); % Ave P
output.pDiff = zeros( 1, numFiles ); % P Max - P Min spaled by input bp
output.pFDWHD = zeros( 1, numFiles ); % P full width half (max - min)
output.pP2P = zeros( 1, numFiles ); % P peak to peak distance
%n
output.nMax = zeros( 1, numFiles ); % Max N
output.nMin = zeros( 1, numFiles ); % Min N
output.nAve = zeros( 1, numFiles ); % Ave N
output.nDiff = zeros( 1, numFiles ); % N Max - N Min snaled by input bn
output.nFDWHD = zeros( 1, numFiles ); % N full width half (max - min)

% loop over files
for ii = 1 : numFiles
  dirname = dirlist( ii ).name;
  % load params and opObj from current WD
  fullpath = [ dirpath '/' dirname ];
  paramsMat = dir( [ fullpath '/param*' ] );
  opMat = dir( [ fullpath '/op*' ] );
  load( [ fullpath '/' paramsMat.name ] );
  fileObj = matfile( [ fullpath   '/' opMat.name ] ); 
  C = fileObj.C_rec(:,:,end);
  P = fileObj.POP_rec(:,:,end);
  N = fileObj.NOP_rec(:,:,end);
  % find direction (x or y) of max variance
  [~, rows, cols, ~, NposVar, Lvar] = ...
    sliceInhomoVar( systemObj, C );
  % record parameters and steady state
  output.p1(ii) = systemObj.l2;
  output.p2(ii) = particleObj.fD;
  output.p3(ii) = systemObj.bc;
  output.steady(ii) = denRecObj.SteadyState;
  % run analysis for C,N,P. Be carefule about P
  [cStats, pStats, nStats] = ...
    bandStatsCPNwrap( C(rows,cols), P(rows,cols), N(rows,cols), Lvar, NposVar );
  % put them in a vector
  %c
  output.cMax(ii) = cStats.maxV;
  output.cMin(ii) = cStats.minV; 
  output.cAve(ii) = cStats.aveV; 
  output.cDiff(ii) = cStats.vdiff; 
  output.cFDWHD(ii) = cStats.fwhd; 
    
  %p
  output.pMax(ii) = pStats.maxV;
  output.pMin(ii) = pStats.minV; 
  output.pAve(ii) = pStats.aveV; 
  output.pDiff(ii) = pStats.vdiff; 
  output.pFDWHD(ii) = pStats.fwhd; 
  output.pP2P(ii) = pStats.p2p; 

  %n
  output.nMax(ii) = nStats.maxV;
  output.nMin(ii) = nStats.minV; 
  output.nAve(ii) = nStats.aveV; 
  output.nDiff(ii) = nStats.vdiff; 
  output.nFDWHD(ii) = nStats.fwhd; 

end % loop over files

% find parameter that is varying
numVaryP = 0;
if length( unique( output.p1 ) ) > 1
  p.val = output.p1;
  p.name = output.p1name;
  numVaryP = numVaryP + 1;
  analysisFlag = 1;
end
if length( unique( output.p2 ) ) > 1
  p.val = output.p2;
  p.name = output.p2name;
  numVaryP = numVaryP + 1;
  analysisFlag = 1;
end
if length( unique( output.p3 ) ) > 1
  p.val = output.p3;
  p.name = output.p3name;
  numVaryP = numVaryP + 1;
  analysisFlag = 1;
end
if  numVaryP == 0
  fprintf('No varying parameters\n');
  plotFlag = 0;
  analysisFlag = 0;
elseif numVaryP > 1
  fprintf('Too many varying parameters\n');
  plotFlag = 0;
  analysisFlag = 0;
end


% Procede if analyzing
if analysisFlag
  % sort the parameters in case they are out of order
  [p.val, indSort] = sort( p.val ); 
  output.p1 = output.p1(indSort);
  output.p2 = output.p2(indSort); 
  output.p3 = output.p3(indSort);
  output.steady = output.steady(indSort);
  output.cMax = output.cMax(indSort); 
  output.cMin = output.cMin(indSort);
  output.cAve = output.cAve(indSort);
  output.cDiff = output.cDiff(indSort);
  output.cFDWHD = output.cFDWHD(indSort);
  %p
  output.pMax = output.pMax(indSort); 
  output.pMin = output.pMin(indSort);
  output.pAve = output.pAve(indSort);
  output.pDiff = output.pDiff(indSort);
  output.pFDWHD = output.pFDWHD(indSort);
  %n
  output.nMax = output.nMax(indSort); 
  output.nMin = output.nMin(indSort);
  output.nAve = output.nAve(indSort);
  output.nDiff = output.nDiff(indSort);
  output.nFDWHD = output.nFDWHD(indSort);
end

% plotting
if plotFlag
  bandPlotMaxDiff( p, output );
  bandPlotFWHM( p, output );
  polarP2Pplot( p, output );
end
