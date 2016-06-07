
% Date
DateTime =  datestr(now);
fprintf('Starting RunHardRod: %s\n', DateTime);

% Add Subroutine path
CurrentDir = pwd;
addpath( genpath( [CurrentDir '/Subroutines'] ) );

% Make Output Directories
if ~exist('runfiles', 'dir'); mkdir('./runfiles'); end;
if ~exist('runOPfiles', 'dir'); mkdir('./runOPfiles'); end;
if ~exist('analyzedfiles', 'dir'); mkdir('./analyzedfiles'); end;

% Grab initial parameters
if exist('Params.mat','file') == 0;
  if exist('InitParams.m','file') == 0;
    cpParams
  end;
  InitParams
end
load Params.mat;

% Copy the master parameter list to ParamObj
ParamObj = ParamMaster;
RhoInit  = RhoInitMaster;
Flags    = FlagMaster;
AnisoDiffFlag = Flags.AnisoDiff;
trial = ParamObj.trial;

% Print what you are doing
if AnisoDiffFlag  == 1;
  fprintf('Anisotropic Hard Rod \n')
else
  fprintf('Isotropic Hard Rod \n')
end

% Fix the time
[TimeObj]= ...
  TimeStepRecMaker(TimeMaster.dt,TimeMaster.t_tot,...
  TimeMaster.t_record,TimeMaster.t_write);
TimeObj.ss_epsilon = TimeMaster.ss_epsilon;

% Display everythin
disp(Flags); disp(ParamObj); disp(TimeObj); disp(RhoInit);

% Make paramMat
[paramMat, numRuns] = MakeParamMat( ParamObj, RhoInit, Flags );
paramvec = zeros(numRuns,1);
% For some reason, param_mat gets "sliced". Create vectors to get arround
paramNx  = paramMat(:,1); paramNy  = paramMat(:,2);
paramNm  = paramMat(:,3); paramLx  = paramMat(:,4);
paramLy  = paramMat(:,5); paramvD  = paramMat(:,6);
parambc  = paramMat(:,7); paramIC  = paramMat(:,8);
paramSM  = paramMat(:,9); paramrun = paramMat(:,10);

% Loops over all run
if numRuns > 1
  parfor ii = 1:numRuns
    % Assign parameters
    paramvec = [ paramNx(ii) paramNy(ii) paramNm(ii) paramLx(ii) ...
      paramLy(ii) paramvD(ii) parambc(ii) paramIC(ii)...
      paramSM(ii) paramrun(ii)];
    
    % Name the file
    filename = [ 'Hr_Ani' num2str( AnisoDiffFlag ) ...
      '_N' num2str( paramNx(ii) ) num2str( paramNy(ii) ) num2str( paramNm(ii) )  ...
      '_Lx' num2str( paramLx(ii) ) 'Ly' num2str( paramLy(ii) )...
      '_vD' num2str( paramvD(ii) ) '_bc' num2str( parambc(ii) ) ...
      '_IC' num2str( paramIC(ii) ) '_SM' num2str( paramSM(ii) ) ...
      '_t' num2str( trial ) '.' num2str( paramrun(ii) ) '.mat' ];
    
    disp(filename);
    
    [DenRecObj] = ...
      HR2DrotMain( filename, paramvec, ParamObj, TimeObj, RhoInit, Flags );
  end
else
  paramvec = [ paramNx(1) paramNy(1) paramNm(1) paramLx(1) ...
    paramLy(1) paramvD(1) parambc(1) paramIC(1)...
    paramSM(1) paramrun(1)];
  
  % Name the file
  filename = [ 'Hr_Ani' num2str( AnisoDiffFlag ) ...
    '_N' num2str( paramNx(1) ) num2str( paramNy(1) ) num2str( paramNm(1) )  ...
    '_Lx' num2str( paramLx(1) ) 'Ly' num2str( paramLy(1) )...
    '_vD' num2str( paramvD(1) ) '_bc' num2str( parambc(1) ) ...
    '_IC' num2str( paramIC(1) ) '_SM' num2str( paramSM(1) ) ...
    '_t' num2str( trial ) '.' num2str( paramrun(1) ) '.mat' ];
  
  disp(filename);
  
  [DenRecObj] = ...
    HR2DrotMain( filename, paramvec, ParamObj, TimeObj, RhoInit, Flags );
end


DateTime =  datestr(now);
fprintf('Finished RunHardRod: %s\n', DateTime);
delete Params.mat

