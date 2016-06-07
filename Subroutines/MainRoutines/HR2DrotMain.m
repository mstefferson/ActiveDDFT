%Obj HR2DrotMain.m
%
% Program is the main() for running the diffusion of 2D hard rods with
% orientation. Program handles interactions using DDFT

function  [ DenRecObj ] = ...
  HR2DrotMain( filename, paramVec, ParamObj, TimeObj, RhoInit, Flags )

global RunSave
MasterSave = 0;
DenRecObj = 0;

try
  % Move parameter vector to obj
  ParamObj.Nx = paramVec(1);
  ParamObj.Ny = paramVec(2);
  ParamObj.Nm = paramVec(3);
  ParamObj.Lx = paramVec(4);
  ParamObj.Ly = paramVec(5);
  ParamObj.vD = paramVec(6);
  ParamObj.bc = paramVec(7);
  RhoInit.IntCond = paramVec(8);
  Flags.StepMeth = paramVec(9);
  ParamObj.runID = paramVec(10);
  ParamObj.Norm  = ParamObj.bc / ( ParamObj.L_rod ^ 2 / pi ) *... % number of particles
    ParamObj.Lx * ParamObj.Ly;
  ParamObj.c = ParamObj.bc ./ ParamObj.b;
  ParamObj.L_box = [ParamObj.Lx ParamObj.Ly];
  
  % Set-up save paths, file names, and matfile
  if Flags.SaveMe
    SaveNameRun   = ['run_' filename];
    
    if Flags.MakeOP == 0
      DirName    =  './runfiles';
    else
      SaveNameOP   = ['OP_' filename];
      DirName  = filename(1:end-4) ;
      if Flags.MakeMovies == 1
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
  
  % Set some flags to 0
  EvolvedDen = 0;DenFinal = 0;DenFTFinal = 0;GridObj = 0;
  DidIBreak = 0;SteadyState = 0;MaxReldRho = 0;
  
  if Flags.Verbose
    if Flags.AnisoDiff
      fprintf('Running Anisotropic Diffusion\n');
    else
      fprintf('Running Isotropic Diffusion\n');
    end
  end
  
  % Record how long things take
  tMainID  = tic;
  
  % Create a file that holds warning print statements
  LocString = sprintf('Location_%d.%d.txt',ParamObj.trialID,ParamObj.runID);
  lfid      = fopen(LocString,'a+');    % a+ allows to append data
  fprintf(lfid,'Starting main, current code\n');
  
  % Make remaining objects
  % Make all the grid stuff %
  tGridID = tic;
  [GridObj] = GridMakerPBCxk(...
    ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);
  GridRunTime = toc(tGridID);
  if Flags.Verbose
    fprintf('Made grid t%d_%d: %.3g \n', ...
      ParamObj.trialID, ParamObj.runID, GridRunTime);
  end
  fprintf(lfid,'Made Grid: %.3g \n', GridRunTime);
  RunTime.Grid = GridRunTime;
  
  %Make diffusion coeff
  tDiffID = tic;
  if Flags.AnisoDiff == 1
    [DiffMobObj] =  DiffMobCoupCoeffCalc( ParamObj.Tmp,...
      ParamObj.Mob_par,ParamObj.Mob_perp,ParamObj.Mob_rot,...
      TimeObj.dt, min(GridObj.dx,GridObj.dy),...
      GridObj.dphi,GridObj.kx2D, GridObj.ky2D,ParamObj.vD);
  else
    [DiffMobObj] = DiffMobCoupCoeffCalcIsoDiff(...
      ParamObj.Tmp,ParamObj.Mob_pos,ParamObj.Mob_rot);
  end
  
  DiffRunTime = toc(tDiffID);
  if Flags.Verbose
    fprintf('Made diffusion object t%d_%d: %.3g\n', ...
      ParamObj.trialID, ParamObj.runID, DiffRunTime);
  end
  fprintf(lfid,'Made diffusion object: %.3g\n', DiffRunTime);
  RunTime.Diff = DiffRunTime;
  
  %Initialze density
  tIntDenID = tic;
  [rho] = MakeConc(GridObj,ParamObj,RhoInit);
  Nc    = 20;
  % Equilib distribution. Don't let bc = 1.5
  if 1.499 < ParamObj.bc && ParamObj.bc < 1.501
    RhoInit.bc = 1.502;
  else
    RhoInit.bc = ParamObj.bc;
  end
  [Coeff_best,~] = CoeffCalcExpCos2D(Nc,GridObj.phi,RhoInit.bc); % Calculate coeff
  RhoInit.feq = DistBuilderExpCos2Dsing(Nc,GridObj.phi,Coeff_best);        % Build equil distribution
  IntDenRunTime = toc(tIntDenID);
  if Flags.Verbose
    fprintf('Made initial density t%d_%d: %.3g \n', ...
      ParamObj.trialID, ParamObj.runID, IntDenRunTime);
  end
  fprintf(lfid,'Made initial density: %.3g\n', IntDenRunTime);
  RunTime.IntDen = IntDenRunTime;
  
  % Save everything before running body of code
  if Flags.SaveMe
    RunSave.Flags    = Flags;
    RunSave.ParamObj = ParamObj;
    RunSave.TimeObj  = TimeObj;
    RunSave.GridObj  = GridObj;
    RunSave.RhoInit  = RhoInit;
    RunSave.ParamObj = ParamObj;
    RunSave.Den_rec = zeros(ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,2);
    RunSave.DenFT_rec = complex( ...
      zeros(ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,2), 0 );
    RunSave.Den_rec(:,:,:,1) = rho;
    RunSave.DenFT_rec(:,:,:,1) = fftshift(fftn(rho));
  end
  
  % Run the main code
  tBodyID      = tic;
  
  if Flags.AnisoDiff == 1
    [DenRecObj]  = HR2DrotDenEvolverFTBody(...
      rho, ParamObj, TimeObj, GridObj, DiffMobObj, Flags, lfid);
  else
    [DenRecObj]  = HR2DrotDenEvolverFTBodyIdC(...
      rho, ParamObj, TimeObj, GridObj, DiffMobObj, Flags, lfid);
  end
  
  % Save it
  if Flags.SaveMe
    RunSave.DenRecObj = DenRecObj;
  end
  
  EvolvedDen = 1;
  BodyRunTime  = toc(tBodyID);
  if Flags.Verbose
    fprintf('Ran Main Body t%d_%d: %.3g \n', ...
      ParamObj.trialID, ParamObj.runID, BodyRunTime);
  end
  fprintf(lfid,'Body Run Time = %f\n\n', BodyRunTime);
  RunTime.Body = BodyRunTime;
  
  
  % Store final density and transform
  DenFinal   = DenRecObj.rhoFinal;
  DidIBreak  = DenRecObj.DidIBreak;
  SteadyState = DenRecObj.SteadyState;
  MaxReldRho  = DenRecObj.MaxReldRho;
  
  % Run movies if you want
  if Flags.MakeOP  == 1
    tOpID           = tic ;
    
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
    OpSave.Flags    = Flags;
    OpSave.ParamObj = ParamObj;
    OpSave.TimeObj  = TimeObj;
    OpSave.C_rec    = zeros(ParamObj.Nx, ParamObj.Ny, 2);
    OpSave.POP_rec  = zeros(ParamObj.Nx, ParamObj.Ny, 2);
    OpSave.POPx_rec = zeros(ParamObj.Nx, ParamObj.Ny, 2);
    OpSave.POPy_rec = zeros(ParamObj.Nx, ParamObj.Ny, 2);
    OpSave.NOP_rec  = zeros(ParamObj.Nx, ParamObj.Ny, 2);
    OpSave.NOPx_rec = zeros(ParamObj.Nx, ParamObj.Ny, 2);
    OpSave.NOPy_rec = zeros(ParamObj.Nx, ParamObj.Ny, 2);
    if Flags.MakeMovies
      OPobj.C_rec    = zeros(ParamObj.Nx, ParamObj.Ny, totRec);
      OPobj.POP_rec  = zeros(ParamObj.Nx, ParamObj.Ny, totRec);
      OPobj.POPx_rec = zeros(ParamObj.Nx, ParamObj.Ny, totRec);
      OPobj.POPy_rec = zeros(ParamObj.Nx, ParamObj.Ny, totRec);
      OPobj.NOP_rec  = zeros(ParamObj.Nx, ParamObj.Ny, totRec);
      OPobj.NOPx_rec = zeros(ParamObj.Nx, ParamObj.Ny, totRec);
      OPobj.NOPy_rec = zeros(ParamObj.Nx, ParamObj.Ny, totRec);
    end
    
    % Break it into chunks
    NumChunks = TimeObj.N_chunks;
    SizeChunk = floor( totRec/ NumChunks );
    NumChunks = ceil( totRec/ SizeChunk);
    
    for i = 1:NumChunks;
      if i ~= NumChunks
        Ind =  (i-1) * SizeChunk + 1: i * SizeChunk;
      else
        Ind = (i-1) * SizeChunk:totRec;
      end
      
      [OPObjTemp] = CPNrecMaker(ParamObj.Nx,ParamObj.Ny,...
        TimeRecVecTemp(Ind) ,GridObj,...
        RunSave.Den_rec(:,:,:,Ind) );
      
      % Save it
      OpSave.C_rec(:,:,Ind) = OPObjTemp.C_rec;
      OpSave.POP_rec(:,:,Ind) = OPObjTemp.POP_rec;
      OpSave.POPx_rec(:,:,Ind) = OPObjTemp.POPx_rec;
      OpSave.POPy_rec(:,:,Ind) = OPObjTemp.POPy_rec;
      OpSave.NOP_rec(:,:,Ind) = OPObjTemp.NOP_rec;
      OpSave.NOPx_rec(:,:,Ind) = OPObjTemp.NOPx_rec;
      OpSave.NOPy_rec(:,:,Ind) = OPObjTemp.NOPy_rec;
      
      if Flags.MakeMovies
        OPobj.C_rec(:,:,Ind)    = OPObjTemp.C_rec;
        OPobj.POP_rec(:,:,Ind)  = OPObjTemp.POP_rec;
        OPobj.POPx_rec(:,:,Ind) = OPObjTemp.POPx_rec;
        OPobj.POPy_rec(:,:,Ind) = OPObjTemp.POPy_rec;
        OPobj.NOP_rec(:,:,Ind)  = OPObjTemp.NOP_rec;
        OPobj.NOPx_rec(:,:,Ind) = OPObjTemp.NOPx_rec;
        OPobj.NOPy_rec(:,:,Ind) = OPObjTemp.NOPy_rec;
      end
      
    end % loop over chunks
    
    [~,~,~,~,OpSave.NOPeq,~,~] = ...
      OpCPNCalc(1, 1, RhoInit.feq, GridObj.phi, 1, 1, GridObj.phi3D);
    if Flags.MakeMovies;
      OPobj.OpTimeRecVec = TimeRecVecTemp;
      OPobj.NOPeq = OpSave.NOPeq;
    end
    
    OpRunTime = toc(tOpID);
    if Flags.Verbose
      fprintf('Made OP object t%d_%d: %.3g \n', ...
        ParamObj.trialID, ParamObj.runID, OpRunTime);
    end
    fprintf(lfid,'OrderParam Run time = %f\n', OpRunTime);
    RunTime.OP = OpRunTime;
    
    if Flags.MakeMovies == 1
      MovieSuccess = 0;
      
      % Make matlab movies
      tMovID       = tic;
      %         keyboard
      HoldX = ParamObj.Nx /2 + 1; % spatial pos placeholders
      HoldY = ParamObj.Ny /2 + 1; % spatial pos placeholders
      DistRec =  reshape( RunSave.Den_rec(HoldX, HoldY, : , :),...
        [ParamObj.Nm length(DenRecObj.TimeRecVec)] );
      
      % Save Name
      MovStr = sprintf('OPmov%d.%d.avi',ParamObj.trialID,ParamObj.runID);
      
      OPMovieMakerTgtherDirAvi(MovStr,...
        GridObj.x,GridObj.y,GridObj.phi,OPobj,...
        DistRec,OPobj.OpTimeRecVec);
      
      MovieSuccess = 1;
      % Move it
      movefile( MovStr, DirName  )
      
      MovRunTime   = toc(tMovID);
      if Flags.Verbose
        fprintf('Made movies t%d_%d: %.3g \n', ...
          ParamObj.trialID, ParamObj.runID, MovRunTime);
      end
      fprintf(lfid,'Make Mov Run Time = %f\n',  MovRunTime);
      RunTime.Mov = MovRunTime;
      
      % Make amplitude plot
      kx0 = ParamObj.Nx / 2 + 1;
      ky0 = ParamObj.Ny / 2 + 1;
      km0 = ParamObj.Nm / 2 + 1;
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
      ampPlotterFT(FTmat2plot, FTind2plot, DenRecObj.TimeRecVec, ParamObj.Nx, ParamObj.Ny,...
        ParamObj.Nm, ParamObj.bc,ParamObj.vD, ParamObj.trialID)
      
      % Save it
      figtl = sprintf('AmpFT_%d_%d',ParamObj.trialID, ParamObj.runID);
      savefig(gcf,figtl)
      saveas(gcf, figtl,'jpg')
      
      % Move it
      movefile([figtl '.fig'], DirName  )
      movefile([figtl '.jpg'], DirName  )
      
    end % End if movies
    
  end % if OP
  
  % Save how long everything took
  TotRunTime = toc(tMainID);
  if Flags.Verbose
    fprintf('Run Finished t%d_%d: %.3g \n', ...
      ParamObj.trialID, ParamObj.runID, TotRunTime);
  end
  fprintf(lfid,'Total Run time = %f\n', TotRunTime);
  RunTime.Tot = TotRunTime;
  
  % Move saved things
  
  if Flags.SaveMe
    RunSave.RunTime = RunTime;
    movefile(SaveNameRun,DirName);
    
    if Flags.MakeOP == 1
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
  
  if Flags.SaveMe
    if Flags.MakeMovies == 1
      if MovieSuccess == 0
        fprintf('Movies failed\n');
        if exist(MovStr,'file'); delete(MovStr); end
        rmdir(DirName);
        DirName    =  './runOPfiles';
      end
    end
    if Flags.MakeOP == 1
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

if Flags.Verbose
  fprintf('Leaving Main for t%d.%d\n', ...
    ParamObj.trialID, ParamObj.runID);
end

end % End HR2DrotVgrExeMain.m
