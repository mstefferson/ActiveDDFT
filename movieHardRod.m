% movieHardRod
%
% Takes all files in ./runOPfiles, makes movies, and moves them to analyzed
function movieHardRod()

try
  tstart = tic;
  % Add Subroutine path
  currentDir = pwd;
  addpath( genpath( [currentDir '/src'] ) );
  
  %make output directories if they don't exist
  if exist('analyzedfiles','dir') == 0; mkdir('analyzedfiles');end;
  
  % see how many dirs to analyze
  dir2Analyze = dir( './runOPfiles');
  numDirs = length(dir2Analyze) - 2;
  
  if numDirs
    fprintf('Making movies for %d dirs \n', numDirs);
    dir2Analyze = dir2Analyze(3:end);
    for ii = 1:numDirs
      % move into a dir
      dirTemp = dir2Analyze(ii).name;
      dirFullPath = ['./runOPfiles/' dirTemp];
      % load things
      runFileName = [dirFullPath '/run_' dirTemp '.mat'];
      opFileName = [dirFullPath '/op_' dirTemp '.mat'];
      
      runSave = matfile( runFileName);
      opSave  = matfile( opFileName );
      
      OPobj.C_rec    = opSave.C_rec;
      OPobj.POP_rec  = opSave.POP_rec;
      OPobj.POPx_rec = opSave.POPx_rec;
      OPobj.POPy_rec = opSave.POPy_rec;
      OPobj.NOP_rec  = opSave.NOP_rec;
      OPobj.NOPx_rec = opSave.NOPx_rec;
      OPobj.NOPy_rec = opSave.NOPy_rec;
      OPobj.distSlice_rec = opSave.distSlice_rec;
      OPobj.OpTimeRecVec = opSave.OpTimeRecVec;
      
      gridObj  = runSave.gridObj;
      particleObj  = runSave.particleObj;
      systemObj  = runSave.systemObj;
      runObj  = runSave.runObj;
      
      % Save Name
      movStr = sprintf('OPmov_bc%.2f_vD%.1f_%.2d_%.2d.avi',...
        systemObj.bc,particleObj.vD,runObj.trialID, runObj.runID);
      % Call movie routine
      OPMovieMakerTgtherDirAvi(movStr,...
        gridObj.x,gridObj.y,gridObj.phi,OPobj,...
        OPobj.distSlice_rec,OPobj.OpTimeRecVec);
      
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
          runSave.DenFT_rec( FTind2plot(jj,1), FTind2plot(jj,2), FTind2plot(jj,3),: ),...
          [ 1, Nrec ]  );
      end
      
      % Plot Amplitudes
      ampPlotterFT(FTmat2plot, FTind2plot, OPobj.OpTimeRecVec, kx0, ky0, km0);
      
      % Save it
      figtl = sprintf('AmpFT.fig');
      % savefig doesn't like decimals so save it and rename it.
      savefig(gcf,figtl)
      figtl2 = sprintf('AmpFT_bc%.2f_vD%.0f_%.2d_%.2d',...
        systemObj.bc, particleObj.vD,runObj.trialID, runObj.runID);
      movefile(figtl,[figtl2 '.fig'])
      saveas(gcf, [figtl2 '.jpg'],'jpg')
      
      % move it avi and figs into directory
      movefile([movStr '*'], dirFullPath);
      movefile([figtl2 '*'], dirFullPath);
      % move directory
      movefile(dirFullPath, ['./analyzedfiles/' dirTemp ] )
      
    end % loop over dir
  else
    fprintf('Nothing to make movies for \n');
  end
catch err
  fprintf('%s', err.getReport('extended')) ;
end
