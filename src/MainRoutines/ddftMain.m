
% Program is the main() for running the diffusion of 2D hard rods with
% orientation. Program handles interactions using DDFT

function  [ denRecObj ] = ...
  ddftMain( filename, paramVec, systemObj, particleObj, runObj, timeObj, ...
  rhoInit, flags, interactLrV, externalV )
% use latex for plots
set(0,'defaulttextinterpreter','latex')
% some global
global runSave
% init some flags
movieSuccess = 0;
evolvedSucess = 0;
movStr = '';
try
  % Set up denRecObj just in case code doesn't finish
  denRecObj.didIrun = 0;
  % Move parameter vector to obj
  systemObj.n1 = paramVec(1);
  systemObj.n2 = paramVec(2);
  systemObj.n3 = paramVec(3);
  systemObj.l1 = paramVec(4);
  systemObj.l2 = paramVec(5);
  particleObj.vD = paramVec(6);
  systemObj.bc = paramVec(7);
  flags.StepMeth = paramVec(8);
  runObj.runID = paramVec(9);
  particleObj.interactLrV = interactLrV{paramVec(10)};
  particleObj.externalV = externalV{paramVec(11)};
  % set some other parameters
  systemObj.c = systemObj.bc ./ particleObj.b;
  systemObj.numPart  = systemObj.c * systemObj.l1 * systemObj.l2; % number of particles
  systemObj.lBox = [systemObj.l1 systemObj.l2];
  % set up the time object
  dtOrig = timeObj.dt ;
  if timeObj.scaleDt && particleObj.vD ~=0
    timeObj.dt = min( timeObj.dt, timeObj.dt * 30  / particleObj.vD );
  end
  % Fix the time
  ss_epsilon = timeObj.ss_epsilon;
  scaleDt = timeObj.scaleDt;
  [timeObj]= ...
    TimeStepRecMaker(timeObj.dt,timeObj.t_tot,...
    timeObj.t_rec,timeObj.t_write);
  % Scale ss_epsilon by delta_t. Equivalent to checking d rho /dt has reached
  % steady state instead of d rho
  timeObj.ss_epsilon = ss_epsilon;
  timeObj.ss_epsilon_dt = ss_epsilon .* timeObj.dt;
  timeObj.scaleDt = scaleDt;
  timeObj.dt_orig = dtOrig;
  timeObj.recStartInd = 2; % start at 2 since t = 0 is ind = 1
  % Set-up save paths, file names, and matfile
  if flags.SaveMe
    saveNameRun   = ['run_' filename];
    saveNameRhoFinal   = ['rhoFinal_' filename];
    saveNameParams = ['params_' filename];
    dirName  = filename(1:end-4) ;
    if flags.MakeOP == 0
      dirName  = ['./runfiles/' dirName ];
    else
      saveNameOP   = ['op_' filename];
      if flags.MakeMovies == 1
        dirName  = ['./analyzedfiles/' dirName ];
      else
        dirName  = ['./runOPfiles/' dirName ];
      end
      opSave = matfile(saveNameOP,'Writable',true);
    end
    if exist(dirName,'dir') == 0
      mkdir(dirName);
    end
    paramSave = matfile(saveNameParams,'Writable',true);
    rhoFinalSave = matfile(saveNameRhoFinal,'Writable',true);
    runSave = matfile(saveNameRun,'Writable',true);
    denRecObj.dirName = dirName; % Just in case
  else
    dirName = pwd;
  end
  % verbose
  if flags.Verbose
    if flags.AnisoDiff
      fprintf('Running Anisotropic Diffusion\n');
    else
      fprintf('Running Isotropic Diffusion\n');
    end
  end
  % Record how long things take
  tMainID  = tic;
  % Create a file that holds warning print statements
  locString = sprintf('Loc_%s.txt', filename(1:end-4));
  lfid      = fopen(locString,'a+');    % a+ allows to append data
  fprintf(lfid,'Starting main, current code\n');
  % Make remaining objects
  % Make all the grid stuff %
  tGridID = tic;
  [gridObj] = GridMakerPBCxk(systemObj.n1,systemObj.n2,systemObj.n3,...
    systemObj.l1,systemObj.l2,systemObj.l3);
  gridrunTime = toc(tGridID);
  if flags.Verbose
    fprintf('Made grid t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, gridrunTime);
  end
  fprintf(lfid,'Made Grid: %.3g \n', gridrunTime);
  runTime.grid = gridrunTime;
  %Make diffusion coeff
  tDiffID = tic;
  [diffObj] =  DiffMobCoupCoeffCalc( systemObj.tmp,...
    particleObj.mob,particleObj.mobPar,particleObj.mobPerp,particleObj.mobRot,...
    gridObj.k1, gridObj.k2, gridObj.k3, ...
    gridObj.k1rep2, gridObj.k2rep2,particleObj.vD);
  diffRunTime = toc(tDiffID);
  if flags.Verbose
    fprintf('Made diffusion object t%d_%d: %.3g\n', ...
      runObj.trialID, runObj.runID, diffRunTime);
  end
  fprintf(lfid,'Made diffusion object: %.3g\n', diffRunTime);
  runTime.diff = diffRunTime;
  %Initialze density
  tIntDenID = tic;
  % Find non-driving steady state
  if strcmp( particleObj.type, 'rods') && strcmp( particleObj.interHb, 'mayer' )
    % Number of coefficients
    Nc    = 20;
    % Equilib distribution. Don't let bc = 1.5
    if 1.499 < systemObj.bc && systemObj.bc < 1.501
      rhoInit.bc = 1.502;
    else
      rhoInit.bc = systemObj.bc;
    end
    if systemObj.n3 == 1
      rhoInit.feq = [];
    else
      [Coeff_best,~] = CoeffCalcExpCos2D(Nc,gridObj.x3,rhoInit.bc); % Calculate coeff
      rhoInit.feq = DistBuilderExpCos2Dsing(Nc,gridObj.x3,Coeff_best);        % Build equil distribution
    end
  else
    rhoInit.feq  = 1 / ( systemObj.l3 ) .* ones( systemObj.n3, 1 );
  end
  % Build initial density
  [rho] = MakeConc(systemObj,particleObj,rhoInit,gridObj);
  intDenRunTime = toc(tIntDenID);
  if flags.Verbose
    fprintf('Made initial density t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, intDenRunTime);
  end
  fprintf(lfid,'Made initial density: %.3g\n', intDenRunTime);
  runTime.intDen = intDenRunTime;
  % Set-up interactions and external potentials
  [interObj] =  interObjMaker( particleObj, systemObj, gridObj );
  % set-up driving
  dRhoDriveFlag = flags.Drive && flags.DiagLop;
  polarDrive  = DrhoPolarDriveClass( dRhoDriveFlag, particleObj.vD, systemObj.n1,...
    systemObj.n2, systemObj.n3, gridObj.x3, gridObj.k1rep2, gridObj.k2rep2 );
  % Save everything before running body of code
  if flags.SaveMe
    runSave.flags    = flags;
    runSave.runObj    = runObj;
    runSave.systemObj = systemObj;
    runSave.particleObj = particleObj;
    runSave.timeObj  = timeObj;
    runSave.rhoInit  = rhoInit;
    % Clean up gridobj before saving
    fields2del = {'k1rep2','k2rep2'};
    gridTemp = rmfield(gridObj,fields2del);
    runSave.gridObj  = gridTemp;
    runSave.Den_rec = zeros(systemObj.n1,systemObj.n2,systemObj.n3,2);
    runSave.DenFT_rec = complex( ...
      zeros(systemObj.n1,systemObj.n2,systemObj.n3,2), 0 );
    runSave.Den_rec(:,:,:,1) = rho;
    runSave.DenFT_rec(:,:,:,1) = fftshift(fftn(rho));
    runSave.denRecObj   = denRecObj;
    runSave.numSavedRhos = 1;
    % Clean up interObj before saving
    if isfield( interObj, 'FmFt' )
      runSave.interObj = rmfield(interObj, 'FmFt');
    end
    % Save params now
    paramSave.flags = flags;
    paramSave.particleObj = particleObj;
    paramSave.rhoInit = rhoInit;
    paramSave.systemObj = systemObj;
    paramSave.runObj = runObj;
    paramSave.timeObj = timeObj;
    paramSave.denRecObj = denRecObj;
  end
  % Run the main code
  tBodyID      = tic;
  if flags.DiagLop == 1
    [denRecObj, rho]  = denEvolverFTDiagOp( ...
      rho, systemObj, particleObj, timeObj, gridObj, diffObj, interObj, ...
      polarDrive, flags, lfid);
  else
    [denRecObj, rho]  = denEvolverFT( ...
      rho, systemObj, particleObj, timeObj, gridObj, diffObj, interObj, ...
      polarDrive, flags, lfid);
  end
  bodyRunTime  = toc(tBodyID);
  evolvedSucess = 1;
  % Save it
  if flags.SaveMe
    paramSave.denRecObj = denRecObj;
    runSave.denRecObj   = denRecObj;
    rhoFinalSave.rho   = rho;
  end
  if flags.Verbose
    fprintf('Ran Main Body t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, bodyRunTime);
  end
  fprintf(lfid,'Body Run Time = %f\n\n', bodyRunTime);
  runTime.body = bodyRunTime;
  % Run movies if you want
  if flags.MakeOP  == 1
    tOpID           = tic ;
    % Save
    if systemObj.n3 > 1
      % Commonly used trig functions
      [~,~,phi3D] = meshgrid(gridObj.x2,gridObj.x1,gridObj.x3);
      cosPhi3d = cos(phi3D);
      sinPhi3d = sin(phi3D);
      cos2Phi3d = cosPhi3d .^ 2;
      sin2Phi3d = sinPhi3d .^ 2;
      cossinPhi3d = cosPhi3d .* sinPhi3d;
    else
      cosPhi3d=0;sinPhi3d=0;cos2Phi3d=0;sin2Phi3d=0;cossinPhi3d=0;
    end
    % Build time rec vector
    if denRecObj.DidIBreak == 0
      totRec = length( denRecObj.TimeRecVec );
      opTimeRecVec = denRecObj.TimeRecVec ;
    else %Don't incldue the blowed up denesity for movies. They don't like it.
      totRec = length( denRecObj.TimeRecVec ) - 1;
      opTimeRecVec = denRecObj.TimeRecVec(1:end-1) ;
    end
    opSave.OpTimeRecVec = opTimeRecVec;
    % Set up saving
    opSave.C_rec    = zeros(systemObj.n1, systemObj.n2, 2);
    if systemObj.n3 > 1
      % Distribution slice
      % currently hardcode out inset
      sliceRho.plotInset = 0;
      dotx1 = round( systemObj.n1/4 ); sliceRho.dotx1 = dotx1;
      doty1 = round( 3 * systemObj.n2/4 ); sliceRho.doty1 = doty1;
      dotx2 = round( systemObj.n1/2 ); sliceRho.dotx2 = dotx2;
      doty2 = round( systemObj.n2/2 ); sliceRho.doty2 = doty2;
      dotx3 = round( 3*systemObj.n1/4 ); sliceRho.dotx3 = dotx3;
      doty3 = round( systemObj.n2/4 ); sliceRho.doty3 = doty3;
      nFrames = length( opTimeRecVec ); sliceRho.nFrames = nFrames;
      sliceRho.slice1 = reshape( runSave.Den_rec( dotx1, doty1, :, 1:nFrames ), ...
        [ systemObj.n3 nFrames] );
      sliceRho.slice2 = reshape( runSave.Den_rec( dotx2, doty2, :, 1:nFrames ), ...
        [ systemObj.n3 nFrames] );
      sliceRho.slice3 = reshape( runSave.Den_rec( dotx3, doty3, :, 1:nFrames ), ...
        [ systemObj.n3 nFrames] );
      sliceRho.phi = gridObj.x3;
      opSave.sliceRho = sliceRho;
      opSave.POP_rec  = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.POPx_rec = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.POPy_rec = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.NOP_rec  = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.NOPx_rec = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.NOPy_rec = zeros(systemObj.n1, systemObj.n2, 2);
      opSave.aveC_rec = zeros(1,2);
      opSave.aveP_rec = zeros(1,2);
      opSave.aveN_rec = zeros(1,2);
    end
    % Break it into chunks
    numChunks = timeObj.N_chunks;
    sizeChunk = floor( totRec/ numChunks );
    if sizeChunk > 0
      numChunks = ceil( totRec/ sizeChunk);
    else
      numChunks = 1;
    end
    for i = 1:numChunks
      if i ~= numChunks
        currInd =  (i-1) * sizeChunk + 1: i * sizeChunk;
      else
        if numChunks == 1
          currInd = 1:totRec;
        else
          currInd = (i-1) * sizeChunk:totRec;
        end
      end
      % Make the records
      [OPObjTemp] = CPNrecMaker(systemObj.n1,systemObj.n2,...
        systemObj.n3, opTimeRecVec(currInd), runSave.Den_rec(:,:,:,currInd) ,...
        gridObj.x3,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d );
      % Save it
      opSave.C_rec(:,:,currInd) = OPObjTemp.C_rec;
      opSave.aveC_rec(1,currInd) = OPObjTemp.aveC_rec;
      if systemObj.n3 > 1
        opSave.POP_rec(:,:,currInd) = OPObjTemp.POP_rec;
        opSave.POPx_rec(:,:,currInd) = OPObjTemp.POPx_rec;
        opSave.POPy_rec(:,:,currInd) = OPObjTemp.POPy_rec;
        opSave.NOP_rec(:,:,currInd) = OPObjTemp.NOP_rec;
        opSave.NOPx_rec(:,:,currInd) = OPObjTemp.NOPx_rec;
        opSave.NOPy_rec(:,:,currInd) = OPObjTemp.NOPy_rec;
        opSave.aveC_rec(1,currInd) = OPObjTemp.aveC_rec;
        opSave.aveP_rec(1,currInd) = OPObjTemp.aveP_rec;
        opSave.aveN_rec(1,currInd) = OPObjTemp.aveN_rec;
      end
    end % loop over chunks
    % fix size
    opSave.C_rec = opSave.C_rec(:,:,1:totRec);
    opSave.aveC_rec = opSave.aveC_rec(1,1:totRec);
    if systemObj.n3 > 1
      opSave.POP_rec  = opSave.POP_rec(:,:,1:totRec);
      opSave.POPx_rec = opSave.POPx_rec(:,:,1:totRec);
      opSave.POPy_rec = opSave.POPy_rec(:,:,1:totRec);
      opSave.NOP_rec  = opSave.NOP_rec(:,:,1:totRec);
      opSave.NOPx_rec = opSave.NOPx_rec(:,:,1:totRec);
      opSave.NOPy_rec = opSave.NOPy_rec(:,:,1:totRec);
      opSave.aveP_rec = opSave.aveP_rec(1,1:totRec);
      opSave.aveN_rec = opSave.aveN_rec(1,1:totRec);
    end
    if flags.MakeMovies
      OPobj.OpTimeRecVec = opTimeRecVec;
      OPobj.C_rec    = opSave.C_rec;
      if systemObj.n3 > 1
        OPobj.POP_rec  = opSave.POP_rec;
        OPobj.POPx_rec = opSave.POPx_rec;
        OPobj.POPy_rec = opSave.POPy_rec;
        OPobj.NOP_rec  = opSave.NOP_rec;
        OPobj.NOPx_rec = opSave.NOPx_rec;
        OPobj.NOPy_rec = opSave.NOPy_rec;
      end
    end
    % Now do it for steady state sol
    if systemObj.n3 > 1
      [~,~,phi3D] = meshgrid(1,1,gridObj.x3);
      cosPhi3d = cos(phi3D);
      sinPhi3d = sin(phi3D);
      cos2Phi3d = cosPhi3d .^ 2;
      sin2Phi3d = sinPhi3d .^ 2;
      cossinPhi3d = cosPhi3d .* sinPhi3d;
      % Calc CPN
      [~,~,~,~,opSave.NOPeq,~,~] = ...
        OpCPNCalc(1, 1, reshape( rhoInit.feq, [1,1,systemObj.n3] ), ...
        gridObj.x3,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d);
      % Save if movies
      if flags.MakeMovies
        OPobj.NOPeq = opSave.NOPeq;
      end
    end
    opRunTime = toc(tOpID);
    if flags.Verbose
      fprintf('Made OP object t%d_%d: %.3g \n', ...
        runObj.trialID, runObj.runID, opRunTime);
    end
    fprintf(lfid,'OrderParam Run time = %f\n', opRunTime);
    runTime.op = opRunTime;
    % Make movies
    if flags.MakeMovies == 1
      movieSuccess = 0;
      % Make matlab movies
      tMovID       = tic;
      % for n3 = 1, just concentration
      if systemObj.n3 == 1
        % Save Name
        movStr = sprintf('Cmov_bc%.2f_vD%.1f_%.2d_%.2d.avi',...
          systemObj.bc,particleObj.vD,runObj.trialID, runObj.runID);
        % Run function
        CMovieMakerAvi(movStr,...
          gridObj.x1,gridObj.x2,particleObj.b .* OPobj.C_rec,...
          OPobj.OpTimeRecVec);
      else
        % Save Name
        movStr = sprintf('OPmov_bc%.2f_vD%.1f_%.2d_%.2d.avi',...
          systemObj.bc,particleObj.vD,runObj.trialID, runObj.runID);
        % Run function
        OPMovieMakerTgtherDirAvi(movStr,...
          gridObj.x1,gridObj.x2, OPobj, ...
          OPobj.OpTimeRecVec, particleObj.b, sliceRho);
      end
      movieSuccess = 1;
      % Move it
      movRunTime   = toc(tMovID);
      if flags.Verbose
        fprintf('Made movies t%d_%d: %.3g \n', ...
          runObj.trialID, runObj.runID, movRunTime);
      end
      fprintf(lfid,'Make Mov Run Time = %f\n',  movRunTime);
      runTime.mov = movRunTime;
      % Make amplitude plot
      k1ind0 = gridObj.k1ind0;
      k2ind0 = gridObj.k2ind0;
      k3ind0 = gridObj.k3ind0;
      nRec = length( denRecObj.TimeRecVec);
      % Ft amps
      if systemObj.n3 > 1
        totModes   = 12;
      else
        totModes = 4;
      end
      FTind2plot = zeros( totModes , 3 );
      FTmat2plot = zeros( totModes , nRec );
      FTind2plot(1,:) = [k1ind0     k2ind0     k3ind0 ];
      FTind2plot(2,:) = [k1ind0 + 1 k2ind0     k3ind0 ];
      FTind2plot(3,:) = [k1ind0     k2ind0 + 1 k3ind0 ];
      FTind2plot(4,:) = [k1ind0 + 1 k2ind0 + 1 k3ind0 ];
      if systemObj.n3 > 1
        FTind2plot(5,:) = [k1ind0     k2ind0     k3ind0 + 1];
        FTind2plot(6,:) = [k1ind0 + 1 k2ind0     k3ind0 + 1];
        FTind2plot(7,:) = [k1ind0     k2ind0 + 1 k3ind0 + 1];
        FTind2plot(8,:) = [k1ind0 + 1 k2ind0 + 1 k3ind0 + 1];
        FTind2plot(9,:) = [k1ind0     k2ind0     k3ind0 + 2];
        FTind2plot(10,:) = [k1ind0 + 1 k2ind0     k3ind0 + 2];
        FTind2plot(11,:) = [k1ind0     k2ind0 + 1 k3ind0 + 2];
        FTind2plot(12,:) = [k1ind0 + 1 k2ind0 + 1 k3ind0 + 2];
      end
      % Scale by N so it's N independent
      for i = 1:totModes
        FTmat2plot(i,:) =  1 / (systemObj.n1 * systemObj.n2 * systemObj.n3) .* ...
          reshape(runSave.DenFT_rec( FTind2plot(i,1), FTind2plot(i,2), FTind2plot(i,3),1:nRec ),...
          [ 1, nRec ]  );
      end
      % Plot Amplitudes
      ampPlotterFT(FTmat2plot, FTind2plot, ...
        denRecObj.TimeRecVec(1:nRec), k1ind0, k2ind0, k3ind0, timeObj.t_tot);
      % Save it
      figtl = sprintf('AmpFT.fig');
      % savefig doesn't like decimals so save it and rename it.
      savefig(gcf,figtl)
      figtl2 = sprintf('AmpFT_bc%.2f_vD%.0f_%.2d_%.2d',...
        systemObj.bc, particleObj.vD,runObj.trialID, runObj.runID);
      movefile(figtl,[figtl2 '.fig'])
      saveas(gcf, [figtl2 '.jpg'],'jpg')
      % Plot final slices of final order parameters
      if systemObj.n3 > 1
        sliceSaveTag = sprintf('SOP_bc%.2f_vD%.0f_%.2d_%.2d',...
          systemObj.bc, particleObj.vD,runObj.trialID, runObj.runID);
        if denRecObj.DidIBreak == 0
          sliceOPplot( OPobj.C_rec(:,:,end), OPobj.POP_rec(:,:,end),...
            OPobj.NOP_rec(:,:,end), systemObj, ...
            gridObj, rhoFinalSave.rho, sliceSaveTag )
        else
          stepsNb = length(OPobj.OpTimeRecVec);
          sliceOPplot( OPobj.C_rec(:,:,end), OPobj.POP_rec(:,:,end),...
            OPobj.NOP_rec(:,:,end), systemObj, ...
            gridObj, runSave.Den_rec( :,:,:, stepsNb), sliceSaveTag );
        end
        % Move it
        movefile([sliceSaveTag '*'], dirName);
        % Plot max order parameters vs time
        maxSaveTag = sprintf('MaxOP_bc%.2f_vD%.0f_%.2d_%.2d',...
          systemObj.bc, particleObj.vD,runObj.trialID, runObj.runID);
        plotMaxOPvsTime( OPobj.C_rec, OPobj.POP_rec, OPobj.NOP_rec, ...
          particleObj.b, OPobj.OpTimeRecVec, maxSaveTag );
      else
        % Plot max order parameters vs time
        maxSaveTag = sprintf('MaxC_bc%.2f_vD%.0f_%.2d_%.2d',...
          systemObj.bc, particleObj.vD,runObj.trialID, runObj.runID);
        plotMaxCvsTime( OPobj.C_rec, particleObj.b, OPobj.OpTimeRecVec, maxSaveTag );
      end
      % Move everything else
      movefile([figtl2 '*'], dirName);
      movefile([movStr '*'], dirName);
      movefile([maxSaveTag '*'], dirName);
    end % End if movies
  end % if OP
  % Save how long everything took
  totRunTime = toc(tMainID);
  if flags.Verbose
    fprintf('Run Finished t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, totRunTime);
  end
  fprintf(lfid,'Total Run time = %f\n', totRunTime);
  runTime.tot = totRunTime;
  % Move saved things
  if flags.SaveMe
    paramSave.runTime = runTime;
    movefile(saveNameRun,dirName);
    movefile( saveNameRhoFinal,dirName);
    movefile( saveNameParams,dirName);
    if flags.MakeOP == 1
      movefile( saveNameOP,dirName);
    end
  end
catch err %Catch errors
  % write the error to file and to screen
  fprintf('%s\n', err.getReport('extended')) ;
  runSave.err = err;
  % Movies can have issues to box size. If they do, just move files
  % to ./runOPfiles
  % Move saved things
  if evolvedSucess == 1 && flags.SaveMe == 1
    if flags.MakeMovies == 1 && movieSuccess == 0
      fprintf('Movies failed\n');
      if exist(movStr,'file'); delete(movStr); end
      rmdir(dirName);
      dirName = [ './runOPfiles' filename(1:end-4) ];
      denRecObj.dirName = dirName;
    end
    movefile(saveNameRun,dirName);
    movefile( saveNameParams,dirName);
    if flags.MakeOP == 1
      movefile( saveNameOP,dirName);
    end
  end
end %End try and catch
% close files
fclose(lfid);
delete(locString);
if flags.Verbose
  fprintf('Leaving Main for t%d.%d\n', ...
    runObj.trialID, runObj.runID);
end
end % End ddftMain
