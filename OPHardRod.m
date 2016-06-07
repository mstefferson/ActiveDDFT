% OPHardRod
%
% Calculates the Order parameters for all the ./runfiles'
function OPHardRod( NumFiles2Analyze )
tstart = tic;
% Add Subroutine path
CurrentDir = pwd;
addpath( genpath( [CurrentDir '/Subroutines'] ) );

if nargin == 0; NumFiles2Analyze = 1; end;

%Make sure it's not a string (bash)
if isa(NumFiles2Analyze,'string');
  fprintf('You gave me a string, turning it to an int\n');
  NumFiles2Analyze = str2int('NumFiles2Analyze');
end;

%make output directories if they don't exist
if exist('runOPfiles','dir') == 0; mkdir('runOPfiles');end;
if exist('./runfiles/analyzing','dir') == 0; mkdir('./runfiles/analyzing');end;

%grab files
Files2Analyze = filelist( '.mat', './runfiles');
NumFilesTot = size(Files2Analyze,1);

%Fix issues if Numfiles is less than desired amount
if NumFiles2Analyze > NumFilesTot;
  NumFiles2Analyze = NumFilesTot;
end;

% Move the files you want to analyze to an analyzing folder
if NumFiles2Analyze;
  fprintf('Moving files to analyzing directory\n');
  %
  for ii=1:NumFiles2Analyze
    % Grab a file
    filename = Files2Analyze{ii};
    movefile( ['./runfiles/' filename], ['./runfiles/analyzing/' filename] );
  end
  
  fprintf('Starting analysis\n');
  
  
  for ii=1:NumFiles2Analyze
    
    % Grab a file
    SaveNameRun = Files2Analyze{ii};
    fprintf('Analyzing %s\n',SaveNameRun);
    
    % Put all variables in a struct
    RunSave = matfile( ['./runfiles/analyzing/' SaveNameRun] );
    %     RunSave = matfile( ['./runfiles/' SaveNameRun] );
    DenRecObj = RunSave.DenRecObj;
    ParamObj  = RunSave.ParamObj;
    TimeObj  = RunSave.TimeObj;
    Flags  = RunSave.Flags;
    RhoInit  = RunSave.RhoInit;
    GridObj  = RunSave.GridObj;
    
    DirName = SaveNameRun(5:end-4);
    SaveNameOP   = ['OP_' DirName '.mat' ];
    OpSave = matfile(SaveNameOP,'Writable',true);
    
    DirName  = ['./runOPfiles/' DirName ];
    
    if exist(DirName ,'dir') == 0;
      mkdir(DirName);
    end
    
    if  DenRecObj.DidIBreak == 0
      totRec = length( DenRecObj.TimeRecVec );
      OpTimeRecVec = DenRecObj.TimeRecVec ;
      OpSave.OpTimeRecVec = OpTimeRecVec;
    else %Don't incldue the blowed up denesity for movies. They don't like it.
      totRec = length( DenRecObj.TimeRecVec ) - 1;
      OpTimeRecVec = DenRecObj.TimeRecVec(1:end-1) ;
      OpSave.OpTimeRecVec = OpTimeRecVec;
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
    
    % Analyze chucks in parallel
    Nx = ParamObj.Nx; Ny = ParamObj.Ny;
    
    % Break it into chunks
    NumChunks = TimeObj.N_chunks;
    SizeChunk = floor( totRec/ NumChunks );
    NumChunks = ceil( totRec/ SizeChunk);
    
    %OpSave.NOPy_rec = zeros(ParamObj.Nx, ParamObj.Ny, 2);
    C_rec    = zeros(ParamObj.Nx, ParamObj.Ny,  SizeChunk);
    POP_rec  = zeros(ParamObj.Nx, ParamObj.Ny,  SizeChunk);
    POPx_rec = zeros(ParamObj.Nx, ParamObj.Ny,  SizeChunk);
    POPy_rec = zeros(ParamObj.Nx, ParamObj.Ny,  SizeChunk);
    NOP_rec  = zeros(ParamObj.Nx, ParamObj.Ny,  SizeChunk);
    NOPx_rec = zeros(ParamObj.Nx, ParamObj.Ny,  SizeChunk);
    NOPy_rec = zeros(ParamObj.Nx, ParamObj.Ny,  SizeChunk);
    
    for jj = 1:NumChunks;
      if jj ~= NumChunks
        Ind =  (jj-1) * SizeChunk + 1: jj * SizeChunk;
      else
        Ind = (jj-1) * SizeChunk:totRec;
      end
      
      DenRecTemp = RunSave.Den_rec(:,:,:,Ind);
      TimeRecVecTemp = OpTimeRecVec(Ind);
      
      if length(Ind) ~= SizeChunk;
        C_rec    = zeros(ParamObj.Nx, ParamObj.Ny,  length(Ind) );
        POP_rec  = zeros(ParamObj.Nx, ParamObj.Ny,  length(Ind) );
        POPx_rec = zeros(ParamObj.Nx, ParamObj.Ny,  length(Ind) );
        POPy_rec = zeros(ParamObj.Nx, ParamObj.Ny,  length(Ind) );
        NOP_rec  = zeros(ParamObj.Nx, ParamObj.Ny,  length(Ind) );
        NOPx_rec = zeros(ParamObj.Nx, ParamObj.Ny,  length(Ind) );
        NOPy_rec = zeros(ParamObj.Nx, ParamObj.Ny,  length(Ind));
      end
      
      
      parfor kk = 1:length(Ind);
        
        [OPObjTemp] = CPNrecMaker(Nx,Ny,...
          TimeRecVecTemp(kk) ,GridObj,...
          DenRecTemp(:,:,:,kk) );
        
        C_rec(:,:,kk) = OPObjTemp.C_rec;
        POP_rec(:,:,kk) = OPObjTemp.POP_rec;
        POPx_rec(:,:,kk) = OPObjTemp.POPx_rec;
        POPy_rec(:,:,kk) = OPObjTemp.POPy_rec;
        NOP_rec(:,:,kk) = OPObjTemp.NOP_rec;
        NOPx_rec(:,:,kk) = OPObjTemp.NOPx_rec;
        NOPy_rec(:,:,kk) = OPObjTemp.NOPy_rec;
      end % loop over chunks
      
      OpSave.C_rec(:,:,Ind)    = C_rec;
      OpSave.POP_rec(:,:,Ind)  = POP_rec;
      OpSave.POPx_rec(:,:,Ind) = POPx_rec;
      OpSave.POPy_rec(:,:,Ind) = POPy_rec;
      OpSave.NOP_rec(:,:,Ind)  = NOP_rec;
      OpSave.NOPx_rec(:,:,Ind) = NOPx_rec;
      OpSave.NOPy_rec(:,:,Ind) = NOPy_rec;
      
    end
    [~,~,~,~,OpSave.NOPeq,~,~] = ...
      OpCPNCalc(1, 1, RhoInit.feq, ...
      GridObj.phi, 1, 1, GridObj.phi3D);
    
    movefile( ['./runfiles/analyzing/' SaveNameRun], DirName );
    movefile( SaveNameOP,DirName );
    
  end %loop over files
end %if analyzing

end_time = toc(tstart);
fprintf('Finished making OPs. OP made for %d files in %.2g min\n', ...
  NumFiles2Analyze, end_time / 60);

% end




