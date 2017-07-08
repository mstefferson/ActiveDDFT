% T = summarizeResults( path2dirs, genComment, plotPeaks, saveMe, saveId )
% Summarize run results. Put them all in a table.
function T = summarizeResults( path2dirs, genComment, plotPeaks, saveMe, saveId )
try
  if nargin < 4
    saveMe = 0; saveId = [];
  else
    saveRoot = 'runSum';
    if nargin == 4
      saveId = [];
    else
      saveId = ['_' saveId];
    end
  end
  % add path
  currentDir = pwd;
  addpath( genpath( [currentDir '/src'] ) );
  % get files
  dirs2analyze = dir( [path2dirs '/Hr_*'] );
  numFiles = length( dirs2analyze );
  if numFiles > 0
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
    numModesPert1 = zeros( numFiles, 1 );
    numModesPert2 = zeros( numFiles, 1 );
    numModesPert3 = zeros( numFiles, 1 );
    bc = zeros( numFiles, 1 );
    c = zeros( numFiles, 1 );
    vD = zeros( numFiles, 1 );
    steady = zeros( numFiles, 1 );
    broken = zeros( numFiles, 1 );
    whatBroke = cell( numFiles, 1 );
    kPosNumPeaks = zeros( numFiles, 1);
    kPosDist = cell( numFiles, 1);
    kPosPeaks = cell( numFiles, 1);
    kPosKa = cell( numFiles, 1);
    realLatticeSpacing = zeros( numFiles, 1);
    spatPhase = cell( numFiles, 1);
    kAngNumPeaks = zeros( numFiles, 1);
    kAngDist = cell( numFiles, 1);
    kAngPeaks = cell( numFiles, 1);
    kAngKa = cell( numFiles, 1);
    angLatticeSpacing = zeros( numFiles, 1);
    angPhase = cell( numFiles, 1);
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
      if ~isempty( particleObj.interLr )
        longL1(ii) =  particleObj.lrLs1;
        longL2(ii) = particleObj.lrLs2;
        longE1(ii) = particleObj.lrEs1;
        longE2(ii) = particleObj.lrEs2;
      else
        longL1(ii) = 0;
        longL2(ii) = 0;
        longE1(ii) = 0;
        longE2(ii) = 0;
      end
      n1(ii) = systemObj.n1;
      n2(ii) = systemObj.n2;
      n3(ii) = systemObj.n3;
      l1(ii) = systemObj.l1;
      l2(ii) = systemObj.l2;
      l3(ii) = systemObj.l3;
      dt(ii) = timeObj.dt;
      numModesPert1(ii) = rhoInit.NumModesX;
      numModesPert2(ii) = rhoInit.NumModesY;
      numModesPert3(ii) = rhoInit.NumModesM;
      tTot(ii) = timeObj.t_tot;
      ssEpsilon(ii) =  timeObj.ss_epsilon;
      bc(ii) = systemObj.bc;
      c(ii) = systemObj.c;
      vD(ii) = particleObj.vD;
      steady(ii) = denRecObj.SteadyState;
      broken(ii) = denRecObj.DidIBreak;
      whatBroke{ii} = denRecObj.whatBroke;
      comments{ii} = genComment;
      % Get spatial peaks
      if isfield( denRecObj, 'rhoFinal' )
        rhoFinal = denRecObj.rhoFinal;
      else
        loadMe = dir( [fileTemp '/rhoFinal_*'] );
        load( [fileTemp '/' loadMe.name ] );
        rhoFinal = rho;
      end
      if systemObj.n3 == 1
        C = rhoFinal;
      else
        dx3 = systemObj.l3 / systemObj.n3;
        x3 = 0:dx3:systemObj.l3 - dx3;
        C = trapz_periodic( x3, rhoFinal, 3 );
      end
      cryPeaks = getCrystalPeaks( C, systemObj.l1, plotPeaks );
      kPosNumPeaks(ii) = cryPeaks.numPeaksK;
      kPosDist{ii} = cryPeaks.dist2PeaksPos;
      kPosPeaks{ii} = cryPeaks.kVecPos;
      kPosKa{ii} = cryPeaks.ka;
      realLatticeSpacing(ii) = cryPeaks.a;
      % record phase. just two for now
      if ~isfinite(cryPeaks.a)
        spatPhase{ii} = 'Liquid';
      elseif cryPeaks.numPeaksK == 3
        spatPhase{ii} = 'Band';
      elseif cryPeaks.numPeaksK == 7
        if cryPeaks.ka > 4
          spatPhase{ii} = 'Crystal B';
        else
          spatPhase{ii} = 'Crystal A';
        end
      else
        spatPhase{ii} = 'Unknown';
      end
      % angular
      if systemObj.n3 > 1
        dx1 = systemObj.l1 / systemObj.n1;
        dx2 = systemObj.l2 / systemObj.n2;
        x1 = 0:dx1:systemObj.l1 - dx1;
        x2 = 0:dx2:systemObj.l2 - dx2;
        f = reshape( trapz_periodic( x1, trapz_periodic( x2, rhoFinal, 2 ), 1 ),...
          [1 systemObj.n3] );
      else
        f = 1;
      end
      cryPeaks = getCrystalPeaks( f, systemObj.l3, plotPeaks );
      kAngNumPeaks(ii) = cryPeaks.numPeaksK;
      kAngDist{ii} = cryPeaks.dist2PeaksPos;
      kAngPeaks{ii} = cryPeaks.kVecPos;
      kAngKa{ii} = cryPeaks.ka;
      angLatticeSpacing(ii) = cryPeaks.a;
      % record angular phase. just two for now
      if cryPeaks.numPeaksK == 1
        angPhase{ii} = 'Isotropic';
      else
        angPhase{ii} = 'Nematic';
      end
    end
    % create a table
    kPosNumPeaks(ii) = cryPeaks.numPeaksK;
    kPosDist{ii} = cryPeaks.dist2PeaksPos;
    kPosPeaks{ii} = cryPeaks.kVecPos;
    kPosKa{ii} = cryPeaks.ka;
    T = table( counter, filename, particleType, shortRange, longRange, diagOp, ...
      iC, sM, trId, rId, longL1, longL2, longE1, longE2, n1, n2, n3,...
      l1, l2, l3, dt, tTot, ssEpsilon, bc, c, vD, steady, broken, whatBroke, ...
      kPosNumPeaks, kPosDist, kPosPeaks, realLatticeSpacing, kPosKa,...
      kAngNumPeaks, kAngDist, kAngPeaks, angLatticeSpacing, kAngKa,...
      spatPhase, angPhase, comments );
  else
    T = [];
    fprintf('No files homie\n')
  end
  % save if
  if saveMe
    saveName = [saveRoot saveId];
    save( [saveName '.mat'], 'T' );
    writetable( T, [saveName '.csv'] );
    movefile( [saveName '*'], path2dirs )
  end
catch err
  fprintf('%s\n', err.getReport('extended')) ;
  keyboard
end

