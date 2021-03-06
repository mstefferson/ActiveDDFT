% movieHardRod
%
% Takes all files in ./runOPfiles, makes movies, and moves them to analyzed
function movieHardRod(desiredFilesId, destinationDir, flagsOverride)
if nargin == 0
  desiredFilesId = 'Hr_*';
  destinationDir = 'analyzedfiles/';
  override = 0;
elseif nargin == 1
  destinationDir  = 'analyzedfiles/';
  override = 0;
elseif nargin == 2
  if ~isempty( destinationDir )
    if destinationDir(end) ~= '/'
      destinationDir = [ destinationDir '/'];
    end
  end
  destinationDir = [ 'analyzedfiles/' destinationDir ];
  override = 0;
else
  if ~isempty( destinationDir )
    if destinationDir(end) ~= '/'
      destinationDir = [ destinationDir '/'];
    end
  end
  destinationDir = [ 'analyzedfiles/' destinationDir ];
  override = 1;
  movieflags = flagsOverride;
end
% use latex for plots
set(0,'defaulttextinterpreter','latex')
try
  % Add Subroutine path
  currentDir = pwd;
  addpath( genpath( [currentDir '/src'] ) );
  %make output directories if they don't exist
  if exist('analyzedfiles','dir') == 0; mkdir('analyzedfiles');end
  if exist(destinationDir,'dir') == 0; mkdir(destinationDir);end
  % see how many dirs to analyze
  if exist( ['runOPfiles/' desiredFilesId], 'dir' )
    dir2Analyze = dir( ['runOPfiles/' desiredFilesId(1:end-1) '*'] );
  else
    dir2Analyze = dir( ['runOPfiles/' desiredFilesId] );
  end
  numDirs = length(dir2Analyze);
  % run if there are things to run
  if numDirs
    fprintf('Making movies for %d dirs \n', numDirs);
    for ii = 1:numDirs
      % move into a dir
      dirTemp = dir2Analyze(ii).name;
      fprintf('Movies for %s\n', dirTemp);
      dirFullPath = ['runOPfiles/' dirTemp];
      % load things
      runFileName = ls( [dirFullPath '/run_*'] );
      opFileName = ls( [dirFullPath '/op_*' ] );
      rhoFinalFileName = ls( [dirFullPath '/rhoFinal_*'] );
      % get rid of blank
      runFileName = runFileName(1:end-1);
      opFileName = opFileName(1:end-1);
      rhoFinalFileName = rhoFinalFileName(1:end-1);
      % build mat files
      runSave = matfile( runFileName );
      opSave  = matfile( opFileName );
      gridObj  = runSave.gridObj;
      particleObj  = runSave.particleObj;
      systemObj  = runSave.systemObj;
      runObj  = runSave.runObj;
      denRecObj = runSave.denRecObj;
      timeObj =  runSave.timeObj;
      % common save str
      commonSaveTag = sprintf('_bc%.2f_fD%.0f_%.2d_%.2d',...
        systemObj.bc, particleObj.fD,runObj.trialID, runObj.runID);
      if override == 0
        flags = runSave.flags;
        if isfield( flags, 'movie')
          movieflags = flags.movie;
        else
          movieflags.makeMovie = 1;
          movieflags.plotInset = 0;
          movieflags.plotMax = 1;
          movieflags.plotAmp = 1;
          movieflags.plotCrystal = 0;
          movieflags.plotSlice = 0;
        end
      end
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
      end
      if movieflags.makeMovie
        if systemObj.n3 > 1
          if isfield(opSave, 'sliceRho')
            sliceRho = opSave.sliceRho;
            sliceRho.plotInset = movieflags.plotInset;
          else
            sliceRho.plotInset = 0;
          end
          % Save Name
          movStr = ['OPmov' commonSaveTag];
          % Make movie
          OPMovieMakerTgtherDirAvi(movStr,...
            gridObj.x1,gridObj.x2, OPobj,...
            OPobj.OpTimeRecVec, particleObj.b, sliceRho);
        else
          % Save Name
          movStr = ['Cmov' commonSaveTag]; 
          % make movie
          CMovieMakerAvi(movStr,...
            gridObj.x1,gridObj.x2,particleObj.b .* OPobj.C_rec,...
            OPobj.OpTimeRecVec);
        end
        % move it
        movefile([movStr '*'], dirFullPath);
      end % plotMov
      if movieflags.plotAmp
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
        figtl2 = ['AmpFT' commonSaveTag];
        movefile(figtl,[figtl2 '.fig'])
        saveas(gcf, [figtl2 '.jpg'],'jpg')
        movefile([figtl2 '*'], dirFullPath);
      end % plotAmp
      % Plot final slices of final order parameters
      if movieflags.plotSlice
        sliceSaveTag = ['SOP' commonSaveTag];
        sliceOPplot( OPobj.C_rec(:,:,end), OPobj.POP_rec(:,:,end),...
          OPobj.NOP_rec(:,:,end), systemObj, ...
          gridObj, rhoFinal, sliceSaveTag )
        movefile([sliceSaveTag '*'], dirFullPath);
      end % plot slice
      % plot max OPs
      if movieflags.plotMax
        if systemObj.n3 > 1
          % Plot max order parameters vs time
          maxSaveTag = ['MaxOP' commonSaveTag];
          plotMaxOPvsTime( OPobj.C_rec, OPobj.POP_rec, OPobj.NOP_rec, ...
            particleObj.b, OPobj.OpTimeRecVec, maxSaveTag );
        else
          % Plot max order parameters vs time
          maxSaveTag = ['MaxC' commonSaveTag];
          plotMaxCvsTime( OPobj.C_rec, particleObj.b, OPobj.OpTimeRecVec, maxSaveTag );
        end
        % move it
        movefile([maxSaveTag '*'], dirFullPath);
      end % plotMax
      % crystal peaks
      if movieflags.plotCrystal
        plotConfirm = 0;
        nFrames = length(OPobj.OpTimeRecVec);
        if isfield(particleObj, 'interLr')
          if strcmp( particleObj.interLr, 'softshoulder') && nFrames > 1
            plotConfirm = 1;
            lrEs1 = particleObj.lrEs1;
            lrEs2 = particleObj.lrEs2;
            lrLs1 = particleObj.lrLs1;
            lrLs2 = particleObj.lrLs2;
          end
        elseif isfield(particleObj, 'interactLrV')
          for jj = 1:length(particleObj.interactLrV)
            if strcmp( particleObj.interactLrV{jj}{1}, 'ss2d' ) && nFrames > 1
              plotConfirm = 1;
              lrEs1 = particleObj.interactLrV{jj}{2}(1);
              lrEs2 = particleObj.interactLrV{jj}{2}(2);
              lrLs1 = particleObj.interactLrV{jj}{2}(3);
              lrLs2 = particleObj.interactLrV{jj}{2}(4);
            end
          end
        end
        if plotConfirm == 1
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
            OPobj.OpTimeRecVec, systemObj, lrEs1, lrEs2, lrLs1, lrLs2, 1 );
          movefile( 'kAmps*', dirFullPath );
        end
      end % plot Crystal
      % move directory
      % dont override
      if exist( [destinationDir  dirTemp], 'dir' )
        dirNew = [  dirTemp '_' ...
          datestr(now,'yyyymmdd_HH.MM.SS') ];
        movefile( dirFullPath, ['runOPfiles/' dirNew] );
        dirFullPath = ['runOPfiles/' dirNew];
      end
      movefile(dirFullPath, destinationDir )
    end % loop over dir
  else
    fprintf('Nothing to make movies for \n');
  end
catch err
  fprintf('%s', err.getReport('extended')) ;
end
