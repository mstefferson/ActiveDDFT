% movieHardRod
%
% Takes all files in ./runOPfiles, makes movies, and moves them to analyzed
function movieHardRod()
% use latex for plots
set(0,'defaulttextinterpreter','latex')
try
  % Add Subroutine path
  currentDir = pwd;
  addpath( genpath( [currentDir '/src'] ) );
  %make output directories if they don't exist
  if exist('analyzedfiles','dir') == 0; mkdir('analyzedfiles');end
  % see how many dirs to analyze
  dir2Analyze = dir( './runOPfiles/Hr_*');
  numDirs = length(dir2Analyze);
  % run if there are things to run
  if numDirs
    fprintf('Making movies for %d dirs \n', numDirs);
    for ii = 1:numDirs
      % move into a dir
      dirTemp = dir2Analyze(ii).name;
      fprintf('Movies for %s\n', dirTemp);
      dirFullPath = ['./runOPfiles/' dirTemp];
      % load things
      runFileName = [dirFullPath '/run_' dirTemp '.mat'];
      opFileName = [dirFullPath '/op_' dirTemp '.mat'];
      rhoFinalFileName = [dirFullPath '/rhoFinal_' dirTemp '.mat'];
      runSave = matfile( runFileName);
      opSave  = matfile( opFileName );
      gridObj  = runSave.gridObj;
      particleObj  = runSave.particleObj;
      systemObj  = runSave.systemObj;
      runObj  = runSave.runObj;
      denRecObj = runSave.denRecObj;
      timeObj =  runSave.timeObj;
      % get rho Final
      if isfield( denRecObj, 'rhoFinal')
        rhoFinal = denRecObj.rhoFinal;
      else
        load(rhoFinalFileName);
        rhoFinal = rho;
      end
      % op stuff
      OPobj.OpTimeRecVec = opSave.OpTimeRecVec;
      OPobj.C_rec    = opSave.C_rec;
      if systemObj.n3 > 1
        OPobj.POP_rec  = opSave.POP_rec;
        OPobj.POPx_rec = opSave.POPx_rec;
        OPobj.POPy_rec = opSave.POPy_rec;
        OPobj.NOP_rec  = opSave.NOP_rec;
        OPobj.NOPx_rec = opSave.NOPx_rec;
        OPobj.NOPy_rec = opSave.NOPy_rec;
        OPobj.distSlice_rec = opSave.distSlice_rec;
        % Save Name
        movStr = sprintf('OPmov_bc%.2f_vD%.1f_%.2d_%.2d.avi',...
          systemObj.bc,particleObj.vD,runObj.trialID, runObj.runID);
        % make movie
        OPMovieMakerTgtherDirAvi(movStr,...
          gridObj.x1,gridObj.x2,gridObj.x3,OPobj,...
          OPobj.distSlice_rec,OPobj.OpTimeRecVec);
      else
        % Save Name
        movStr = sprintf('Cmov_bc%.2f_vD%.1f_%.2d_%.2d.avi',...
          systemObj.bc,particleObj.vD,runObj.trialID, runObj.runID);
        % make movie
        CMovieMakerAvi(movStr,...
          gridObj.x1,gridObj.x2,particleObj.b .* OPobj.C_rec,...
          OPobj.OpTimeRecVec);
      end
      % Make amplitude plot
      kx0 = systemObj.n1 / 2 + 1;
      ky0 = systemObj.n2 / 2 + 1;
      km0 = floor( systemObj.n3 / 2 ) + 1;
      nRec = length( denRecObj.TimeRecVec);
      % fill in amps
      if systemObj.n3 > 1
        totModes   = 12;
        FTind2plot = zeros( totModes , 3 );
        FTmat2plot = zeros( totModes , nRec );
        FTind2plot(1,:) = [kx0     ky0     km0 ];
        FTind2plot(2,:) = [kx0 + 1 ky0     km0 ];
        FTind2plot(3,:) = [kx0     ky0 + 1 km0 ];
        FTind2plot(4,:) = [kx0 + 1 ky0 + 1 km0 ];
        FTind2plot(5,:) = [kx0     ky0     km0 + 1];
        FTind2plot(6,:) = [kx0 + 1 ky0     km0 + 1];
        FTind2plot(7,:) = [kx0     ky0 + 1 km0 + 1];
        FTind2plot(8,:) = [kx0 + 1 ky0 + 1 km0 + 1];
        FTind2plot(9,:) = [kx0     ky0     km0 + 2];
        FTind2plot(10,:) = [kx0 + 1 ky0     km0 + 2];
        FTind2plot(11,:) = [kx0     ky0 + 1 km0 + 2];
        FTind2plot(12,:) = [kx0 + 1 ky0 + 1 km0 + 2];
      else
        totModes   = 4;
        FTind2plot = zeros( totModes , 3 );
        FTmat2plot = zeros( totModes , nRec );
        FTind2plot(1,:) = [kx0   ky0  1 ];
        FTind2plot(2,:) = [kx0 + 1 ky0  1  ];
        FTind2plot(3,:) = [kx0     ky0 + 1 1];
        FTind2plot(4,:) = [kx0 + 1 ky0 + 1 1 ];
      end
      for i = 1:totModes
        % Scale by N so it's N independent
        FTmat2plot(i,:) =  1 / (systemObj.n1 * systemObj.n2 * systemObj.n3) .* ...
          reshape(runSave.DenFT_rec( FTind2plot(i,1), FTind2plot(i,2), FTind2plot(i,3),1:nRec ),...
          [ 1, nRec ]  );
      end
      % Plot Amplitudes
      ampPlotterFT(FTmat2plot, FTind2plot, ...
        denRecObj.TimeRecVec(1:nRec), kx0, ky0, km0, timeObj.t_tot);
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
        if strcmp( particleObj.interHb, 'mayer' )
          sliceSaveTag = sprintf('SOP_bc%.2f_vD%.0f_%.2d_%.2d',...
            systemObj.bc, particleObj.vD,runObj.trialID, runObj.runID);
          sliceOPplot( OPobj.C_rec(:,:,end), OPobj.POP_rec(:,:,end),...
            OPobj.NOP_rec(:,:,end), systemObj, ...
            gridObj, rhoFinal, sliceSaveTag )
          movefile([sliceSaveTag '*'], dirFullPath);
        end
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
      % crystal peaks
      nFrames = length(OPobj.OpTimeRecVec);
      if strcmp( particleObj.interLr, 'softshoulder' ) && nFrames > 1
        
        % get concentration in k space
        if systemObj.n3 == 1
          cFt = reshape( runSave.DenFT_rec(:,:,1,1:nFrames) , ...
            [systemObj.n1, systemObj.n2, nFrames ] );
          cFinal = rhoFinal;
        else
          cFt = zeros( systemObj.n1, systemObj.n2, nFrames );
          for nn = 1:nFrames
            cFt(:,:,nn) = fftshift( fftn( opSave.C_rec(:,:,nn) ) );
          end
          cFinal = opSave.C_rec(:,:,end);
        end
        plotCrystalPeaks( cFt, cFinal, gridObj.k1, gridObj.k2, ...
          OPobj.OpTimeRecVec, systemObj, particleObj, 1 );
        movefile( 'kAmps*', dirFullPath );
      end
      % move it avi and figs into directory
      movefile([movStr '*'], dirFullPath);
      movefile([figtl2 '*'], dirFullPath);
      movefile([maxSaveTag '*'], dirFullPath);
      % move directory
      movefile(dirFullPath, ['./analyzedfiles/' dirTemp ] )
    end % loop over dir
  else
    fprintf('Nothing to make movies for \n');
  end
catch err
  fprintf('%s', err.getReport('extended')) ;
end
