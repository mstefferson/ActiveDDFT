function [output] = bandAnalysis( dirlist )

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
output.cFDWHD = zeros( 1, numFiles ); % C full width half (max - min)
% p
output.pMax = zeros( 1, numFiles ); % Max P
output.pMin = zeros( 1, numFiles ); % Min P
output.pAve = zeros( 1, numFiles ); % Ave P
output.pDiff = zeros( 1, numFiles ); % P Max - P Min spaled by input bp
output.pFDWHD = zeros( 1, numFiles ); % P full width half (max - min)
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
  paramsMat = dir( [ dirname '/param*' ] );
  opMat = dir( [ dirname '/op*' ] );
  load( [ dirname '/' paramsMat.name ] );
  fileObj = matfile( [ dirname   '/' opMat.name ] ); 
  C = fileObj.C_rec(:,:,end);
  P = fileObj.POP_rec(:,:,end);
  N = fileObj.NOP_rec(:,:,end);
  % find direction (x or y) of max variance
  [~, rows, cols, ~, NposVar, Lvar] = ...
    sliceInhomoVar( systemObj, C );
  % record parameters and steady state
  output.p1(ii) = systemObj.Ly;
  output.p2(ii) = particleObj.vD;
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

% Max and diff
figure()
% C max and diff
subplot(3,1,1);
ylab = 'C';
[Ax] = plotyy( p.val, output.cMax, p.val, output.cDiff );
xlabel(p.name); ylabel(Ax(1), [ ylab ' Max']); ylabel(Ax(2), [ ylab ' diff']);
title( [ ylab ' vs ' p.name ] );
% P max and diff
subplot(3,1,2);
ylab = 'P';
[Ax] = plotyy( p.val, output.pMax, p.val, output.pDiff );
xlabel(p.name); ylabel(Ax(1), [ ylab ' Max']); ylabel(Ax(2), [ ylab ' diff']);
title( [ ylab ' vs ' p.name ] );
ylab = 'Max P';
xlabel(p.name); ylabel(ylab);
title( [ ylab ' vs ' p.name ] );
% N max and diff
subplot(3,1,3);
ylab = 'N';
[Ax] = plotyy( p.val, output.nMax, p.val, output.nDiff );
xlabel(p.name); ylabel(Ax(1), [ ylab ' Max']); ylabel(Ax(2), [ ylab ' diff']);
title( [ ylab ' vs ' p.name ] );

% fwhm
figure()
% C fwhm
subplot(3,1,1);
ylab = 'C';
plot( p.val, output.cFDWHD );
xlabel(p.name); ylabel( [ylab ' FWHM']);
title( [ ylab ' vs ' p.name ] );
% P fwhm
subplot(3,1,2);
ylab = 'P';
plot( p.val, output.pFDWHD );
xlabel(p.name); ylabel( [ylab ' FWHM']);
title( [ ylab ' vs ' p.name ] );
% N fwhm
subplot(3,1,3);
ylab = 'N';
plot( p.val, output.nFDWHD );
xlabel(p.name); ylabel( [ylab ' FWHM']);
title( [ ylab ' vs ' p.name ] );

posV =  ( 0: NposVar-1 ) ./  NposVar * Lvar;
figure()
plot(posV, C(rows,cols ) );
title('C')

figure()
plot(posV, P(rows,cols ) );
title('P')

figure()
plot(posV, N(rows,cols ) );
title('N')

keyboard
