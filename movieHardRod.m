% movieHardRod
%
% Takes all files in ./runOPfiles, makes movies, and moves them to analyzed
try
  tstart = tic;
  % Add Subroutine path
  CurrentDir = pwd;
  addpath( genpath( [CurrentDir '/src'] ) );

  %make output directories if they don't exist
  if exist('analyzedfiles','dir') == 0; mkdir('analyzedfiles');end;

  % see how many dirs to analyze
  Dir2Analyze = dir( './runOPfiles');
  numDirs = length(Dir2Analyze) - 2;

  if numDirs
    fprintf('Making movies for %d dirs \n', numDirs);
    Dir2Analyze = Dir2Analyze(3:end);
    for ii = 1:numDirs
      % move into a dir
      dirTemp = Dir2Analyze(ii).name;
      dirFullPath = ['./runOPfiles/' dirTemp];
      % load things
      runFileName = [dirFullPath '/run_' dirTemp '.mat'];
      opFileName = [dirFullPath '/OP_' dirTemp '.mat'];
      
      RunSave = matfile( runFileName);
      OpSave  = matfile( opFileName );
      
      OPobj.C_rec    = OpSave.C_rec;
      OPobj.POP_rec  = OpSave.POP_rec;
      OPobj.POPx_rec = OpSave.POPx_rec;
      OPobj.POPy_rec = OpSave.POPy_rec;
      OPobj.NOP_rec  = OpSave.NOP_rec;
      OPobj.NOPx_rec = OpSave.NOPx_rec;
      OPobj.NOPy_rec = OpSave.NOPy_rec;
      OPobj.OpTimeRecVec = OpSave.OpTimeRecVec;
      
      DenRecObj = RunSave.DenRecObj;
      gridObj  = RunSave.gridObj;
      
      % Make matlab movies
      HoldX = systemObj.Nx /2 + 1; % spatial pos placeholders
      HoldY = systemObj.Ny /2 + 1; % spatial pos placeholders
      
      DistRec =  reshape( RunSave.Den_rec(HoldX, HoldY, : , :),...
        [systemObj.Nm length(OPobj.OpTimeRecVec)] );
      
      % Save Name
      MovStr = sprintf('OPmov%d.%d.avi',runObj.trialID,runObj.runID);
      
      % Call movie routine
      try
        OPMovieMakerTgtherDirAvi(MovStr,...
          gridObj.x,gridObj.y,gridObj.phi,OPobj,...
          DistRec,OPobj.OpTimeRecVec);
        
        MovieSuccess = 1; 
        
        % Make amplitude plot
        kx0 = systemObj.Nx / 2 + 1;
        ky0 = systemObj.Ny / 2 + 1;
        km0 = systemObj.Nm / 2 + 1;
        Nrec = length( OPobj.OpTimeRecVec);
        
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
        
        for jj = 1:8
          FTmat2plot(jj,:) =  reshape(...
            RunSave.DenFT_rec( FTind2plot(jj,1), FTind2plot(jj,2), FTind2plot(jj,3),: ),...
            [ 1, Nrec ]  );
        end
        % Plot Amplitudes
        ampPlotterFT(FTmat2plot, FTind2plot, OPobj.OpTimeRecVec, systemObj.Nx, systemObj.Ny,...
          systemObj.Nm, systemObj.bc,particleObj.vD, runObj.trialID)
        
        % Save it
        figtl = sprintf('AmpFT_%d_%d',runObj.trialID, runObj.runID);
        savefig(gcf,figtl)
        saveas(gcf, figtl,'jpg')
        
        % move it
        movefile([MovStr '*'], dirFullPath);
        movefile([figtl '*'], dirFullPath);
        movefile(dirFullPath, './analyzedfiles');
      catch err
        fprintf('Did not make movie \n');
        disp(err.message);
      end % try catch
    end % loop over dir
  else
    fprintf('Nothing to make movies for \n');
  end
catch err
  throw(err);
end
