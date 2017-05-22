% Summarize run results. Put them all in a table.
function T = summarizeResults( path2dirs, genComment, plotPeaks )
try
  % add path
  currentDir = pwd;
  addpath( genpath( [currentDir '/src'] ) );
  % get files
  dirs2analyze = dir( [path2dirs '/Hr_*'] );
  numFiles = length( dirs2analyze );
  % allocate
  counter = zeros( numFiles, 1 );
  filename = cell( numFiles, 1 );
  particleType = cell( numFiles, 1 );
  shortRange = cell( numFiles, 1 );
  longRange = cell( numFiles, 1 );
  diagOp = zeros( numFiles, 1 );
  iC = zeros( numFiles, 1 );
  sM = zeros( numFiles, 1 );
  trId = zeros( numFiles, 1 );
  rId = zeros( numFiles, 1 );
  longL1 = zeros( numFiles, 1 );
  longL2 = zeros( numFiles, 1 );
  longE1 = zeros( numFiles, 1 );
  longE2 = zeros( numFiles, 1 );
  n1 = zeros( numFiles, 1 );
  n2 = zeros( numFiles, 1 );
  n3 = zeros( numFiles, 1 );
  l1 = zeros( numFiles, 1 );
  l2 = zeros( numFiles, 1 );
  l3 = zeros( numFiles, 1 );
  dt = zeros( numFiles, 1 );
  tTot = zeros( numFiles, 1 );
  ssEpsilon = zeros( numFiles, 1 );
  bc = zeros( numFiles, 1 );
  vD = zeros( numFiles, 1 );
  steady = zeros( numFiles, 1 );
  broken = zeros( numFiles, 1 );
  whatBroke = cell( numFiles, 1 );
  kPosNumPeaks = zeros( numFiles, 1);
  kPosDist = cell( numFiles, 1);
  kPosPeakInds = cell( numFiles, 1);
  spatPhase = cell( numFiles, 1);
  angPhase = cell( numFiles, 1);
  realLatticeSpacing = zeros( numFiles, 1);
  kAngNumPeaks = zeros( numFiles, 1);
  kAngDist = cell( numFiles, 1);
  kAngPeakInds = cell( numFiles, 1);
  angLatticeSpacing = zeros( numFiles, 1);
  comments = cell( numFiles, 1);
  for ii = 1:numFiles
    % get path to params
    fileTemp = [path2dirs '/' dirs2analyze(ii).name];
    loadMe = dir( [fileTemp '/params_*'] );
    load( [fileTemp '/' loadMe.name ] );
    % store everything
    counter(ii) = ii;
    filename{ii} = dirs2analyze(ii).name;
    particleType{ii} = particleObj.type;
    shortRange{ii} = particleObj.interHb;
    longRange{ii} = particleObj.interLr;
    diagOp(ii) = flags.DiagLop;
    iC(ii) = rhoInit.IntCond;
    sM(ii) = flags.StepMeth;
    trId(ii) = runObj.trialID;
    rId(ii) = runObj.runID;
    longL1(ii) =  particleObj.lrLs1;
    longL2(ii) = particleObj.lrLs2;
    longE1(ii) = particleObj.lrEs1;
    longE2(ii) = particleObj.lrEs2;
    n1(ii) = systemObj.n1;
    n2(ii) = systemObj.n2;
    n3(ii) = systemObj.n3;
    l1(ii) = systemObj.l1;
    l2(ii) = systemObj.l2;
    l3(ii) = systemObj.l3;
    dt(ii) = timeObj.dt;
    tTot(ii) = timeObj.t_tot;
    ssEpsilon(ii) =  timeObj.ss_epsilon;
    bc(ii) = systemObj.bc;
    vD(ii) = particleObj.vD;
    steady(ii) = denRecObj.SteadyState;
    broken(ii) = denRecObj.DidIBreak;
    whatBroke{ii} = denRecObj.whatBroke;
    comments{ii} = genComment;
    % Get spatial peaks
    if systemObj.n3 == 1
      C = denRecObj.rhoFinal;
    else
      dx3 = systemObj.l3 / systemObj.n3;
      x3 = 0:dx3:systemObj.l3 - dx3;
      C = trapz_periodic( x3, denRecObj.rhoFinal, 3 );
    end
    peaks = getCrystalPeaks( C, systemObj.l1, plotPeaks );
    kPosNumPeaks(ii) = peaks.numPeaksK;
    kPosDist{ii} = peaks.distK;
    kPosPeakInds{ii} = peaks.indsK;
    realLatticeSpacing(ii) = peaks.a;
    % record phase. just two for now
    if peaks.numPeaksK == 1
      spatPhase{ii} = 'Liquid';
    elseif peaks.numPeaksK == 3
      spatPhase{ii} = 'Band';
    elseif peaks.numPeaksK == 7
      spatPhase{ii} = 'Crystal';
    else
      spatPhase{ii} = 'Unknown';
    end
    % angular
    if systemObj.n3 > 1
      dx1 = systemObj.l1 / systemObj.n1;
      dx2 = systemObj.l2 / systemObj.n2;
      x1 = 0:dx1:systemObj.l1 - dx1;
      x2 = 0:dx2:systemObj.l2 - dx2;
      f = trapz_periodic( x1, trapz_periodic( x2, denRecObj.rhoFinal, 2 ), 1 );
    else
      f = 1;
    end
    peaks = getCrystalPeaks( f, systemObj.l3, plotPeaks );
    kAngNumPeaks(ii) = peaks.numPeaksK;
    kAngDist{ii} = peaks.distK;
    kAngPeakInds{ii} = peaks.indsK;
    angLatticeSpacing(ii) = peaks.a;
    % record angular phase. just two for now
    if peaks.numPeaksK == 1
      angPhase{ii} = ['Isotropic'];
    else
      angPhase{ii} = ['Nematic'];
    end
  end
  % create a table
  T = table( counter, filename, particleType, shortRange, longRange, diagOp, ...
    iC, sM, trId, rId, longL1, longL2, longE1, longE2, n1, n2, n3,...
    l1, l2, l3, dt, tTot, ssEpsilon, bc, vD, steady, broken, whatBroke, ...
    kPosNumPeaks, kPosDist, kPosPeakInds, realLatticeSpacing, ...
    kAngNumPeaks, kAngDist, kAngPeakInds, angLatticeSpacing, ...
    spatPhase, angPhase, comments );
catch err
  fprintf('%s', err.getReport('extended')) ;
  keyboard
end
