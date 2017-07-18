% opHardRod
% Calculates the Order parameters for all the ./runfiles'
function opHardRod( numFiles2Analyze )
try
  tstart = tic;
  % Add Subroutine path
  CurrentDir = pwd;
  addpath( genpath( [CurrentDir '/src'] ) );
  if nargin == 0; numFiles2Analyze = 1; end
  %Make sure it's not a string (bash)
  if isa(numFiles2Analyze,'string')
    fprintf('You gave me a string, turning it to an int\n');
    numFiles2Analyze = str2int('NumFiles2Analyze');
  end
  %make output directories if they don't exist
  if exist('runOPfiles','dir') == 0; mkdir('runOPfiles'); end
  if exist('./runfiles/analyzing','dir') == 0; mkdir('./runfiles/analyzing');end;
  % see how many dirs to analyze
  dir2Analyze = dir( './runfiles/Hr_*');
  numDirs = length(dir2Analyze);
  %Fix issues if Numfiles is less than desired amount
  if numDirs > numFiles2Analyze; numDirs = numFiles2Analyze; end
  % Move the files you want to analyze to an analyzing folder
  if numDirs
    fprintf('Starting analysis\n');
    for ii=1:numDirs
      % Grab a file
      % get dir name
      dirTemp = dir2Analyze(ii).name;
      fprintf('OP for %s\n', dirTemp);
      dirFullPath = ['./runfiles/' dirTemp];
      runFileName = [dirFullPath '/run_' dirTemp '.mat'];
      runSave = matfile( runFileName);
      % Put all variables in a struct
      denRecObj = runSave.denRecObj;
      systemObj  = runSave.systemObj;
      timeObj  = runSave.timeObj;
      rhoInit  = runSave.rhoInit;
      gridObj  = runSave.gridObj;
      n1 = systemObj.n1; n2 = systemObj.n2; n3 = systemObj.n3;
      % Build phi3D once
      [~,~,phi3D] = meshgrid(gridObj.x2,gridObj.x1,gridObj.x3);
      cosPhi3d = cos(phi3D);
      sinPhi3d = sin(phi3D);
      cos2Phi3d = cosPhi3d .^ 2;
      sin2Phi3d = sinPhi3d .^ 2;
      cossinPhi3d = cosPhi3d .* sinPhi3d;
      phi = gridObj.x3;
      % new directory name
      dirFullPath  = ['./runfiles/' dirTemp ];
      saveNameOP   = ['op_' dirTemp '.mat' ];
      % set up new matfile
      OpSave = matfile(saveNameOP,'Writable',true);
      % Old files don't have denRecObj.didIrun. Assume it did for now
      if isfield(denRecObj, 'didIrun' ) == 0
        denRecObj.didIrun = 1;
        runSave.denRecObj = denRecObj;
      end
      % See if it ran or not
      if denRecObj.didIrun == 0
        % If it didn't finish, create time vectors from size
        [~,~,~,totRec] = size( runSave.Den_rec );
        OpTimeRecVec = (0:totRec-1) .* timeObj.t_rec;
        OpSave.OpTimeRecVec = OpTimeRecVec;
        denRecObj.TimeRecVec = OpTimeRecVec;
        denRecObj.DidIBreak = 0;
        denRecObj.SteadyState = 0;
        denRecObj.rhoFinal = runSave.Den_rec(:,:,:,end);
        runSave.denRecObj  = denRecObj;
      else
        if  denRecObj.DidIBreak == 0
          totRec = length( denRecObj.TimeRecVec );
          OpTimeRecVec = denRecObj.TimeRecVec ;
          OpSave.OpTimeRecVec = OpTimeRecVec;
          fprintf('Nothing Broke totRec = %d\n',totRec);
        else %Don't incldue the blowed up density for movies. They don't like it.
          totRec = length( denRecObj.TimeRecVec ) - 1;
          OpTimeRecVec = denRecObj.TimeRecVec(1:end-1) ;
          OpSave.OpTimeRecVec = OpTimeRecVec;
          fprintf('Density Broke totRec = %d\n',totRec);
        end % If it broke
      end % If I ran
      % if only 1 record, skip it
      if totRec == 1
        fprintf('Nothing saved other than IC!!!\n')
        numPoints = 1;
      else
        numPoints = 2;
      end % if totRec == 1
      % initialize
      OpSave.C_rec    = zeros(n1, n2, numPoints);
      OpSave.POP_rec  = zeros(n1, n2, numPoints);
      OpSave.POPx_rec = zeros(n1, n2, numPoints);
      OpSave.POPy_rec = zeros(n1, n2, numPoints);
      OpSave.NOP_rec  = zeros(n1, n2, numPoints);
      OpSave.NOPx_rec = zeros(n1, n2, numPoints);
      OpSave.NOPy_rec = zeros(n1, n2, numPoints);
      % Analyze chucks in parallel
      % Break it into chunks
      numChunks = timeObj.N_chunks;
      sizeChunk = max( floor( totRec/ numChunks ), 1 );
      numChunks = ceil( totRec/ sizeChunk);
      %OpSave.NOPy_rec = zeros(systemObj.n1, systemObj.n2, 2);
      C_rec    = zeros(n1, n2,  sizeChunk);
      POP_rec  = zeros(n1, n2,  sizeChunk);
      POPx_rec = zeros(n1, n2,  sizeChunk);
      POPy_rec = zeros(n1, n2,  sizeChunk);
      NOP_rec  = zeros(n1, n2,  sizeChunk);
      NOPx_rec = zeros(n1, n2,  sizeChunk);
      NOPy_rec = zeros(n1, n2,  sizeChunk);
      % print some things
      fprintf('totPoints = %d, numChunks = %d, sizeChunk = %d\n',...
        totRec, numChunks, sizeChunk);
      for jj = 1:numChunks
        if jj ~= numChunks
          ind =  (jj-1) * sizeChunk + 1: jj * sizeChunk;
        else
          ind = (jj-1) * sizeChunk:totRec;
        end
        ind = ind( ind > 0 );
        % Temp variables
        DenRecTemp = runSave.Den_rec(:,:,:,ind);
        TimeRecVecTemp = OpTimeRecVec(ind);
        if length(ind) ~= sizeChunk
          C_rec    = zeros(n1, n2,  length(ind) );
          POP_rec  = zeros(n1, n2,  length(ind) );
          POPx_rec = zeros(n1, n2,  length(ind) );
          POPy_rec = zeros(n1, n2,  length(ind) );
          NOP_rec  = zeros(n1, n2,  length(ind) );
          NOPx_rec = zeros(n1, n2,  length(ind) );
          NOPy_rec = zeros(n1, n2,  length(ind));
        end
        % Make sure it isn't zero of infinite
        if isfinite(DenRecTemp); notInf = 1; else; notInf = 0; end
        if DenRecTemp ~= 0; notZero = 1; else; notZero = 0;end
        if notInf && notZero
          parfor kk = 1:length(ind)
            %Calculate
            [OPObjTemp] = CPNrecMaker(n1,n2,...
              TimeRecVecTemp(kk), DenRecTemp(:,:,:,kk),...
              phi,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d );
            % Store
            C_rec(:,:,kk) = OPObjTemp.C_rec;
            POP_rec(:,:,kk) = OPObjTemp.POP_rec;
            POPx_rec(:,:,kk) = OPObjTemp.POPx_rec;
            POPy_rec(:,:,kk) = OPObjTemp.POPy_rec;
            NOP_rec(:,:,kk) = OPObjTemp.NOP_rec;
            NOPx_rec(:,:,kk) = OPObjTemp.NOPx_rec;
            NOPy_rec(:,:,kk) = OPObjTemp.NOPy_rec;
          end % parloop
        else
          fprintf('Density contains NAN, INF, or zeros. Not analyzing \n')
        end
        % store it
        if totRec > 1
          OpSave.C_rec(:,:,ind)    = C_rec;
          OpSave.POP_rec(:,:,ind)  = POP_rec;
          OpSave.POPx_rec(:,:,ind) = POPx_rec;
          OpSave.POPy_rec(:,:,ind) = POPy_rec;
          OpSave.NOP_rec(:,:,ind)  = NOP_rec;
          OpSave.NOPx_rec(:,:,ind) = NOPx_rec;
          OpSave.NOPy_rec(:,:,ind) = NOPy_rec;
        else
          OpSave.C_rec    = C_rec;
          OpSave.POP_rec  = POP_rec;
          OpSave.POPx_rec = POPx_rec;
          OpSave.POPy_rec = POPy_rec;
          OpSave.NOP_rec  = NOP_rec;
          OpSave.NOPx_rec = NOPx_rec;
          OpSave.NOPy_rec = NOPy_rec;
        end
      end %loop over chunks
      % Distribution slice
      holdX = systemObj.n1 /2 + 1; % spatial pos placeholders
      holdY = systemObj.n2 /2 + 1; % spatial pos placeholders
      OpSave.distSlice_rec = reshape( ...
        runSave.Den_rec(holdX, holdY, : , 1:length(OpTimeRecVec)),...
        [systemObj.n3 length(OpTimeRecVec)] );
      % Now do it for steady state sol
      [~,~,phi3D] = meshgrid(1,1,phi);
      cosPhi3d = cos(phi3D);
      sinPhi3d = sin(phi3D);
      cos2Phi3d = cosPhi3d .^ 2;
      sin2Phi3d = sinPhi3d .^ 2;
      cossinPhi3d = cosPhi3d .* sinPhi3d;
      [~,~,~,~,OpSave.NOPeq,~,~] = ...
        OpCPNCalc(1, 1, reshape( rhoInit.feq, [1,1,n3] ), ...
        phi,cosPhi3d,sinPhi3d,cos2Phi3d,sin2Phi3d,cossinPhi3d);
      [~, ~, o] = size(OpSave.C_rec);
      % Move it
      movefile( saveNameOP,dirFullPath );
      % move directory
      movefile(dirFullPath, ['./runOPfiles/' dirTemp] )
      % use your words
      fprintf('Finished %s\n', dirTemp);
      fprintf('Rec points for C_rec = %d vs totRec = %d\n',o,totRec);
    end %loop over files
    fprintf('Looped over files\n');
  end %if analyzing
  
  end_time = toc(tstart);
  fprintf('Finished making OPs. OP made for %d files in %.2g min\n', ...
    numDirs, end_time / 60);
  
catch err
  fprintf('%s', err.getReport('extended')) ;
  keyboard
end

% end
