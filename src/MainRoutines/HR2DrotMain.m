%Obj HR2DrotMain.m
%
% Program is the main() for running the diffusion of 2D hard rods with
% orientation. Program handles interactions using DDFT

function  [ DenRecObj ] = ...
  HR2DrotMain( filename, paramVec, systemObj, particleObj, runObj, timeObj, rhoInit, flags )

global RunSave
DenRecObj = 0;
MovieSuccess = 0;

try
  % Move parameter vector to obj
  systemObj.Nx = paramVec(1);
  systemObj.Ny = paramVec(2);
  systemObj.Nm = paramVec(3);
  systemObj.Lx = paramVec(4);
  systemObj.Ly = paramVec(5);
  particleObj.vD = paramVec(6);
  systemObj.bc = paramVec(7);
  rhoInit.IntCond = paramVec(8);
  flags.StepMeth = paramVec(9);
  runObj.runID = paramVec(10);
  systemObj.c = systemObj.bc ./ particleObj.b;
  systemObj.numPart  = systemObj.c * systemObj.Lx * systemObj.Ly; % number of particles 
  systemObj.lBox = [systemObj.Lx systemObj.Ly];
  
  % Set-up save paths, file names, and matfile
  if flags.SaveMe
    SaveNameRun   = ['run_' filename];
    if flags.MakeOP == 0
      DirName    =  './runfiles';
    else
      SaveNameOP   = ['OP_' filename];
      DirName  = filename(1:end-4) ;
      if flags.MakeMovies == 1
        DirPath  = ['./analyzedfiles/' DirName ];
        if exist(DirPath,'dir') == 0;
          mkdir('./analyzedfiles', DirName);
        end
        DirName = DirPath;
      else
        DirPath  = ['./runOPfiles/' DirName ];
        if exist(DirPath,'dir') == 0;
          mkdir('./runOPfiles', DirName);
        end
        DirName = DirPath;
      end
      OpSave = matfile(SaveNameOP,'Writable',true);
    end
    RunSave = matfile(SaveNameRun,'Writable',true);
  end
  
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
  LocString = sprintf('Loc_%s.txt', filename(1:end-4));
  lfid      = fopen(LocString,'a+');    % a+ allows to append data
  fprintf(lfid,'Starting main, current code\n');
  
  % Make remaining objects
  % Make all the grid stuff %
  tGridID = tic;
  [gridObj] = GridMakerPBCxk(...
    systemObj.Nx,systemObj.Ny,systemObj.Nm,systemObj.Lx,systemObj.Ly);
  GridrunTime = toc(tGridID);
  if flags.Verbose
    fprintf('Made grid t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, GridrunTime);
  end
  fprintf(lfid,'Made Grid: %.3g \n', GridrunTime);
  runTime.Grid = GridrunTime;
  
  %Make diffusion coeff
  tDiffID = tic;
  if flags.AnisoDiff == 1
    SptSpc = min( gridObj.x(2) - gridObj.x(1) ,...
    gridObj.y(2) - gridObj.y(1) );
    RotSpt = gridObj.phi(2) - gridObj.phi(1);

    [diffObj] =  DiffMobCoupCoeffCalc( systemObj.Tmp,...
      particleObj.mobPar,particleObj.mobPerp,particleObj.mobRot,...
      timeObj.dt, SptSpc, RotSpt,...
      gridObj.kx, gridObj.ky, gridObj.km, ...
      gridObj.kx2D, gridObj.ky2D,particleObj.vD);
  else
    [diffObj] = DiffMobCoupCoeffCalcIsoDiff(...
      systemObj.Tmp,particleObj.mobPos,particleObj.mobRot, ...
      gridObj.kx, gridObj.ky, gridObj.km);
  end
  
  DiffrunTime = toc(tDiffID);
  if flags.Verbose
    fprintf('Made diffusion object t%d_%d: %.3g\n', ...
      runObj.trialID, runObj.runID, DiffrunTime);
  end
  fprintf(lfid,'Made diffusion object: %.3g\n', DiffrunTime);
  runTime.Diff = DiffrunTime;
  
  %Initialze density
  tIntDenID = tic;
  % Find non-driving steady state
  Nc    = 20;
  % Equilib distribution. Don't let bc = 1.5
  if 1.499 < systemObj.bc && systemObj.bc < 1.501
    rhoInit.bc = 1.502;
  else
    rhoInit.bc = systemObj.bc;
  end
  [Coeff_best,~] = CoeffCalcExpCos2D(Nc,gridObj.phi,rhoInit.bc); % Calculate coeff
  rhoInit.feq = DistBuilderExpCos2Dsing(Nc,gridObj.phi,Coeff_best);        % Build equil distribution

  % Build initial density
  [rho] = MakeConc(systemObj,rhoInit,...
    gridObj.x,gridObj.y,gridObj.phi);

  IntDenrunTime = toc(tIntDenID);
  
  if flags.Verbose
    fprintf('Made initial density t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, IntDenrunTime);
  end
  fprintf(lfid,'Made initial density: %.3g\n', IntDenrunTime);
  runTime.IntDen = IntDenrunTime;
  
  % Save everything before running body of code
  % except for Grid because it's currently too bulky

  if flags.SaveMe
    RunSave.flags    = flags;
    RunSave.systemObj = systemObj;
    RunSave.particleObj = particleObj;
    RunSave.timeObj  = timeObj;
    RunSave.rhoInit  = rhoInit;
    RunSave.Den_rec = zeros(systemObj.Nx,systemObj.Ny,systemObj.Nm,2);
    RunSave.DenFT_rec = complex( ...
      zeros(systemObj.Nx,systemObj.Ny,systemObj.Nm,2), 0 );
    RunSave.Den_rec(:,:,:,1) = rho;
    RunSave.DenFT_rec(:,:,:,1) = fftshift(fftn(rho));
  end
  
  % Run the main code
  tBodyID      = tic;
  if flags.AnisoDiff == 1
    [DenRecObj]  = HR2DrotDenEvolverFTBody(...
      rho, systemObj, particleObj, timeObj, gridObj, diffObj, flags, lfid);
  else
    [DenRecObj]  = HR2DrotDenEvolverFTBodyIdC(...
      rho, systemObj, particleObj, timeObj, gridObj, diffObj, flags, lfid);
  end

  % Save it
  if flags.SaveMe
    % Clean up gridobj
    fields2del = {'kx2D','ky2D'};
    gridObj = rmfield(gridObj,fields2del);
    RunSave.gridObj  = gridObj;
    RunSave.DenRecObj = DenRecObj;
  end
  
  BodyrunTime  = toc(tBodyID);
  if flags.Verbose
    fprintf('Ran Main Body t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, BodyrunTime);
  end
  fprintf(lfid,'Body Run Time = %f\n\n', BodyrunTime);
  runTime.Body = BodyrunTime;
  
  
  % Run movies if you want
  if flags.MakeOP  == 1
    tOpID           = tic ;

    [~,~,phi3D] = meshgrid(gridObj.x,gridObj.y,gridObj.phi); 
    cosPhi3d = cos(phi3D);
    sinPhi3d = sin(phi3D);
    cos2Phi3d = cosPhi3d .^ 2;
    sin2Phi3d = sinPhi3d .^ 2;
    cossinPhi3d = cosPhi3d .* sinPhi3d;

    if  DenRecObj.DidIBreak == 0
      totRec = length( DenRecObj.TimeRecVec );
      TimeRecVecTemp = DenRecObj.TimeRecVec ;
      OpSave.OpTimeRecVec = TimeRecVecTemp;
    else %Don't incldue the blowed up denesity for movies. They don't like it.
      totRec = length( DenRecObj.TimeRecVec ) - 1;
      TimeRecVecTemp = DenRecObj.TimeRecVec(1:end-1) ;
      OpSave.OpTimeRecVec = TimeRecVecTemp;
    end
    
    % Set up saving
    OpSave.flags    = flags;
    OpSave.systemObj = systemObj;
    OpSave.particleObj = particleObj;
    OpSave.timeObj  = timeObj;
    OpSave.C_rec    = zeros(systemObj.Nx, systemObj.Ny, 2);
    OpSave.POP_rec  = zeros(systemObj.Nx, systemObj.Ny, 2);
    OpSave.POPx_rec = zeros(systemObj.Nx, systemObj.Ny, 2);
    OpSave.POPy_rec = zeros(systemObj.Nx, systemObj.Ny, 2);
    OpSave.NOP_rec  = zeros(systemObj.Nx, systemObj.Ny, 2);
    OpSave.NOPx_rec = zeros(systemObj.Nx, systemObj.Ny, 2);
    OpSave.NOPy_rec = zeros(systemObj.Nx, systemObj.Ny, 2);
    if flags.MakeMovies
      OPobj.C_rec    = zeros(systemObj.Nx, systemObj.Ny, totRec);
      OPobj.POP_rec  = zeros(systemObj.Nx, systemObj.Ny, totRec);
      OPobj.POPx_rec = zeros(systemObj.Nx, systemObj.Ny, totRec);
      OPobj.POPy_rec = zeros(systemObj.Nx, systemObj.Ny, totRec);
      OPobj.NOP_rec  = zeros(systemObj.Nx, systemObj.Ny, totRec);
      OPobj.NOPx_rec = zeros(systemObj.Nx, systemObj.Ny, totRec);
      OPobj.NOPy_rec = zeros(systemObj.Nx, systemObj.Ny, totRec);
    end
    
    % Break it into chunks
    NumChunks = timeObj.N_chunks;
    SizeChunk = floor( totRec/ NumChunks );
    if SizeChunk > 0
      NumChunks = ceil( totRec/ SizeChunk);
    else
      NumChunks = 1;
    end
    
    for i = 1:NumChunks;
      if i ~= NumChunks
        Ind =  (i-1) * SizeChunk + 1: i * SizeChunk;
      else
        if NumChunks == 1
          Ind = 1:totRec;
        else
        Ind = (i-1) * SizeChunk:totRec;
        end
      end
      
      [OPObjTemp] = CPNrecMaker(systemObj.Nx,systemObj.Ny,...
        TimeRecVecTemp(Ind), RunSave.Den_rec(:,:,:,Ind) ,...
        gridObj.phi,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d );
 
      % Save it
      OpSave.C_rec(:,:,Ind) = OPObjTemp.C_rec;
      OpSave.POP_rec(:,:,Ind) = OPObjTemp.POP_rec;
      OpSave.POPx_rec(:,:,Ind) = OPObjTemp.POPx_rec;
      OpSave.POPy_rec(:,:,Ind) = OPObjTemp.POPy_rec;
      OpSave.NOP_rec(:,:,Ind) = OPObjTemp.NOP_rec;
      OpSave.NOPx_rec(:,:,Ind) = OPObjTemp.NOPx_rec;
      OpSave.NOPy_rec(:,:,Ind) = OPObjTemp.NOPy_rec;
      
      if flags.MakeMovies
        OPobj.C_rec(:,:,Ind)    = OPObjTemp.C_rec;
        OPobj.POP_rec(:,:,Ind)  = OPObjTemp.POP_rec;
        OPobj.POPx_rec(:,:,Ind) = OPObjTemp.POPx_rec;
        OPobj.POPy_rec(:,:,Ind) = OPObjTemp.POPy_rec;
        OPobj.NOP_rec(:,:,Ind)  = OPObjTemp.NOP_rec;
        OPobj.NOPx_rec(:,:,Ind) = OPObjTemp.NOPx_rec;
        OPobj.NOPy_rec(:,:,Ind) = OPObjTemp.NOPy_rec;
      end
      
    end % loop over chunks
    
    % Now do it for steady state sol
    [~,~,phi3D] = meshgrid(1,1,gridObj.phi); 
    cosPhi3d = cos(phi3D);
    sinPhi3d = sin(phi3D);
    cos2Phi3d = cosPhi3d .^ 2;
    sin2Phi3d = sinPhi3d .^ 2;
    cossinPhi3d = cosPhi3d .* sinPhi3d;
   
    [~,~,~,~,OpSave.NOPeq,~,~] = ...
      OpCPNCalc(1, 1, reshape( rhoInit.feq, [1,1,systemObj.Nm] ), ...
      gridObj.phi,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d);
      
    if flags.MakeMovies;
      OPobj.OpTimeRecVec = TimeRecVecTemp;
      OPobj.NOPeq = OpSave.NOPeq;
    end
    
    OprunTime = toc(tOpID);
    if flags.Verbose
      fprintf('Made OP object t%d_%d: %.3g \n', ...
        runObj.trialID, runObj.runID, OprunTime);
    end
    fprintf(lfid,'OrderParam Run time = %f\n', OprunTime);
    runTime.OP = OprunTime;
    
    if flags.MakeMovies == 1
      MovieSuccess = 0;
      
      % Make matlab movies
      tMovID       = tic;
      HoldX = systemObj.Nx /2 + 1; % spatial pos placeholders
      HoldY = systemObj.Ny /2 + 1; % spatial pos placeholders
      DistRec =  reshape( RunSave.Den_rec(HoldX, HoldY, : , :),...
        [systemObj.Nm length(DenRecObj.TimeRecVec)] );
      
      % Save Name
      MovStr = sprintf('OPmov%d.%d.avi',runObj.trialID,runObj.runID);
      
      OPMovieMakerTgtherDirAvi(MovStr,...
        gridObj.x,gridObj.y,gridObj.phi,OPobj,...
        DistRec,OPobj.OpTimeRecVec);
      
      MovieSuccess = 1;
      % Move it
      movefile( MovStr, DirName  )
      
      MovrunTime   = toc(tMovID);
     if flags.Verbose
        fprintf('Made movies t%d_%d: %.3g \n', ...
          runObj.trialID, runObj.runID, MovrunTime);
      end
      fprintf(lfid,'Make Mov Run Time = %f\n',  MovrunTime);
      runTime.Mov = MovrunTime;
      
      % Make amplitude plot
      kx0 = systemObj.Nx / 2 + 1;
      ky0 = systemObj.Ny / 2 + 1;
      km0 = systemObj.Nm / 2 + 1;
      Nrec = length( DenRecObj.TimeRecVec);
      
      FTind2plot = zeros( 8, 3 );
      FTmat2plot = zeros( 8, Nrec );
      
      FTind2plot(1,:) = [kx0     ky0     km0 + 1];
      FTind2plot(2,:) = [kx0 + 1 ky0     km0 + 1];
      FTind2plot(3,:) = [kx0     ky0 + 1 km0 + 1];
      FTind2plot(4,:) = [kx0 + 1 ky0 + 1 km0 + 1];
      FTind2plot(5,:) = [kx0     ky0     km0 + 2];
      FTind2plot(6,:) = [kx0 + 1 ky0     km0 + 2];
      FTind2plot(7,:) = [kx0     ky0 + 1 km0 + 2];
      FTind2plot(8,:) = [kx0 + 1 ky0 + 1 km0 + 2];
      
      for i = 1:8
        FTmat2plot(i,:) =  reshape(...
          RunSave.DenFT_rec( FTind2plot(i,1), FTind2plot(i,2), FTind2plot(i,3),: ),...
          [ 1, Nrec ]  );
      end
      % Plot Amplitudes
      ampPlotterFT(FTmat2plot, FTind2plot, DenRecObj.TimeRecVec, systemObj.Nx, systemObj.Ny,...
        systemObj.Nm, systemObj.bc,particleObj.vD, runObj.trialID)
      
      % Save it
      figtl = sprintf('AmpFT_%d_%d',runObj.trialID, runObj.runID);
      savefig(gcf,figtl)
      saveas(gcf, figtl,'jpg')
      
      % Move it
      movefile([figtl '.fig'], DirName  )
      movefile([figtl '.jpg'], DirName  )
      
    end % End if movies
    
  end % if OP
  
  % Save how long everything took
  TotrunTime = toc(tMainID);
  if flags.Verbose
    fprintf('Run Finished t%d_%d: %.3g \n', ...
      runObj.trialID, runObj.runID, TotrunTime);
  end
  fprintf(lfid,'Total Run time = %f\n', TotrunTime);
  runTime.Tot = TotrunTime;
  
  % Move saved things
  
  if flags.SaveMe
    RunSave.runTime = runTime;
    movefile(SaveNameRun,DirName);
    
    if flags.MakeOP == 1
      movefile( SaveNameOP,DirName);
    end
    
  end
  
catch err %Catch errors
  
  % write the error to file and to screen
  fprintf('%s', err.getReport('extended')) ;
  RunSave.err = err;
  disp(err);
  
  % Movies can have issues to box size. If they do, just move files
  % to ./runOPfiles
  % Move saved things

  if flags.SaveMe
    if flags.MakeMovies == 1
      if MovieSuccess == 0
        fprintf('Movies failed\n');
        if exist(MovStr,'file'); delete(MovStr); end
        rmdir(DirName);
        DirName    =  './runOPfiles';
      end
    end
    if flags.MakeOP == 1
      DirName  = filename(1:end-4) ;
      DirPath  = ['./runOPfiles/' DirName ];
      if exist(DirPath,'dir') == 0;
        mkdir('./runOPfiles', DirName);
      end
      DirName = DirPath;
      movefile( SaveNameOP,DirName);
    end
    movefile(SaveNameRun,DirName);
  end
  
end %End try and catch

fclose(lfid);
delete(LocString);

if flags.Verbose
  fprintf('Leaving Main for t%d.%d\n', ...
    runObj.trialID, runObj.runID);
end

end % End HR2DrotVgrExeMain.m
