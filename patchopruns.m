% patchruns( dir1, dir2 )
% patch together runs dir1, dir2 in runOPfiles
function patchopruns(dir1,dir2)
% patch the dirs together
dirPatched = [dir1(1:end-1) dir2( strfind(dir2,'_t'):end-1) '_patched/'];
% patched files
path2dir = [ 'runOPfiles/' dirPatched ];
if ~exist(path2dir,'dir'); mkdir(path2dir); end
% filenames
filename = [dirPatched(1:end-1) '.mat'];
paramFilePatchName = ['params_' filename];
runFilePatchName = ['run_' filename];
opFilePatchName = ['op_' filename];
rhoFinalFilePatchName = ['rhoFinal_' filename];
paramFilePatch = matfile( paramFilePatchName, 'Writable', true);
opFilePatch = matfile( opFilePatchName, 'Writable', true );
runFilePatch = matfile( runFilePatchName, 'Writable', true );
rhoFinalFilePatch = matfile( rhoFinalFilePatchName, 'Writable', true );
% load file 1
path2dir = [ 'runOPfiles/' dir1 ];
param2load = ls( [path2dir 'params_*'] );
run2load = ls( [path2dir 'run_*'] );
op2load = ls( [path2dir 'op_*'] );
% get rid of blank space
param2load = param2load(1:end-1);
run2load = run2load(1:end-1);
op2load = op2load(1:end-1);
load( param2load );
runFile1 = matfile( run2load, 'Writable', true);
opFile1 = matfile( op2load, 'Writable', true );
% store some parameters
denRecObj = runFile1.denRecObj;
checkParamVec1 = [systemObj.n1 systemObj.n2 systemObj.n3 ...
  systemObj.l1 systemObj.l2 systemObj.bc particleObj.vD  ];
nt1 = length( denRecObj.TimeRecVec );
% nt1op = length( opFile1.OpTimeRecVec );
dt1 = timeObj.dt;
tRec1 = timeObj.t_rec;
tTot1 = timeObj.t_tot;
timeRecVec1 = denRecObj.TimeRecVec;
simTime1 = denRecObj.simTime;
% totRunTime1 = runTime.tot;
runTime1 = denRecObj.runTime;
rhoInit1 = rhoInit;
% grab file 2
path2dir = [ 'runOPfiles/' dir2 ];
param2load = ls( [path2dir 'params_*'] );
run2load = ls( [path2dir 'run_*'] );
op2load = ls( [path2dir 'op_*'] );
rhoFinal2load = ls( [path2dir 'rhoFinal_*'] );
param2load = param2load(1:end-1);
run2load = run2load(1:end-1);
op2load = op2load(1:end-1);
rhoFinal2load = rhoFinal2load(1:end-1);
load( param2load );
% store some things
dt2 = timeObj.dt;
tRec2 = timeObj.t_rec;
tTot2 = timeObj.t_tot;
nt2 = length( denRecObj.TimeRecVec );
% check parameters
checkParamVec2 = [systemObj.n1 systemObj.n2 systemObj.n3 ...
  systemObj.l1 systemObj.l2 systemObj.bc particleObj.vD  ];
% build dir2 files
runFile2 = matfile( run2load, 'Writable', true );
opFile2 = matfile( op2load, 'Writable', true );
rhoFinalFile2 = matfile( rhoFinal2load, 'Writable', true );
proceed = checkSameParameters( checkParamVec1, checkParamVec2 );
% keyboard
if proceed
  timeObjPatch = timeObj;
  timeObjPatch.dt1 = dt1;
  timeObjPatch.t_rec1 = tRec1;
  timeObjPatch.t_tot1 = tTot1;
  timeObjPatch.dt2 = dt2;
  timeObjPatch.t_rec2 = tRec2;
  timeObjPatch.t_tot2 = tTot2;
  timeObjPatch.t_tot = tTot1 + tTot2;
  % for denRecObj, save file 2 with edits
  denRecObjPatch = denRecObj;
  denRecObjPatch.TimeRecVec = [ timeRecVec1 (tTot1 + denRecObj.TimeRecVec(2:end) ) ];
  denRecObjPatch.simTime = denRecObjPatch.simTime + simTime1;
%   runTime.tot = runTime.tot + totRunTime1;
  denRecObjPatch.runTime = denRecObjPatch.runTime + runTime1;
  % store param patched as param 2
  paramFilePatch.systemObj = systemObj;
  paramFilePatch.particleObj = particleObj;
  paramFilePatch.runObj = runObj;
  paramFilePatch.rhoInit = rhoInit1;
  paramFilePatch.flags = flags;
  paramFilePatch.timeObj = timeObjPatch;
  paramFilePatch.runTime = runTime;
  paramFilePatch.denRecObj = denRecObjPatch;
  % run and op stuff
  runFilePatch.systemObj = systemObj;
  runFilePatch.particleObj = particleObj;
  runFilePatch.runObj = runObj;
  runFilePatch.rhoInit = rhoInit1;
  runFilePatch.flags = flags;
  runFilePatch.timeObj = timeObjPatch;
%   runFilePatch.runTime = runTime;
  runFilePatch.gridObj = runFile2.gridObj;
  runFilePatch.denRecObj = denRecObjPatch;
  % Allocate
  runFilePatch.Den_rec = ...
    zeros(systemObj.n1, systemObj.n2, systemObj.n3, nt1 + nt2 -1 );
  runFilePatch.DenFT_rec = ...
    complex( zeros(systemObj.n1, systemObj.n2, systemObj.n3, nt1 + nt2 -1 ) );
  opFilePatch.C_rec = zeros(systemObj.n1, systemObj.n2, nt1 + nt2 -1 );
  opFilePatch.NOP_rec = zeros(systemObj.n1, systemObj.n2, nt1 + nt2 -1 );
  opFilePatch.NOPx_rec = zeros(systemObj.n1, systemObj.n2, nt1 + nt2 -1 );
  opFilePatch.NOPy_rec = zeros(systemObj.n1, systemObj.n2, nt1 + nt2 -1 );
  opFilePatch.POP_rec = zeros(systemObj.n1, systemObj.n2, nt1 + nt2 -1 );
  opFilePatch.POPx_rec = zeros(systemObj.n1, systemObj.n2, nt1 + nt2 -1 );
  opFilePatch.POPy_rec = zeros(systemObj.n1, systemObj.n2, nt1 + nt2 -1 );
  opFilePatch.aveC_rec = zeros(1, nt1 + nt2 -1 );
  opFilePatch.aveP_rec = zeros(1, nt1 + nt2 -1 );
  opFilePatch.aveN_rec = zeros(1, nt1 + nt2 -1 );
  opFilePatch.distSlice_rec = zeros(systemObj.n3, nt1 + nt2 -1 );
  % track which file to save from
  runFileTemp =  runFile1;
  opFileTemp =  opFile1;
  for ii = 1:nt1+nt2-1
    if ii < nt1+1
      ind = ii;
    elseif ii == nt1+1
      runFileTemp =  runFile2;
      opFileTemp =  opFile2;
      ind = ii+1-nt1;
    else
      ind = ii+1-nt1;
    end
    % grab it
    denTemp = runFileTemp.Den_rec(:,:,:,ind);
    denFtTemp = runFileTemp.DenFT_rec(:,:,:,ind);
    cTemp = opFileTemp.C_rec(:,:,ind);
    nTemp = opFileTemp.NOP_rec(:,:,ind);
    nxTemp = opFileTemp.NOPx_rec(:,:,ind);
    nyTemp = opFileTemp.NOPy_rec(:,:,ind);
    pTemp = opFileTemp.POP_rec(:,:,ind);
    pxTemp = opFileTemp.POPx_rec(:,:,ind);
    pyTemp = opFileTemp.POPy_rec(:,:,ind);
    cAveTemp = opFileTemp.aveC_rec(1,ind);
    pAveTemp = opFileTemp.aveP_rec(1,ind);
    nAveTemp = opFileTemp.aveN_rec(1,ind);
    distSliceTemp = opFileTemp.distSlice_rec(:,ind);
    % save it
    runFilePatch.Den_rec(:,:,:,ii) = denTemp;
    runFilePatch.DenFT_rec(:,:,:,ii) = denFtTemp;
    opFilePatch.C_rec(:,:,ii) = cTemp;
    opFilePatch.NOP_rec(:,:,ii) = nTemp;
    opFilePatch.NOPx_rec(:,:,ii) = nxTemp;
    opFilePatch.NOPy_rec(:,:,ii) = nyTemp;
    opFilePatch.POP_rec(:,:,ii) = pTemp;
    opFilePatch.POPx_rec(:,:,ii) = pxTemp;
    opFilePatch.POPy_rec(:,:,ii) = pyTemp;
    opFilePatch.aveC_rec(1,ii) = cAveTemp;
    opFilePatch.aveP_rec(1,ii) = pAveTemp;
    opFilePatch.aveN_rec(1,ii) = nAveTemp;
    opFilePatch.distSlice_rec(:,ii) = distSliceTemp;
    fprintf('%.1f done', ii / (nt1+nt2-1) );
  end %loop over time
  opFilePatch.OpTimeRecVec = denRecObjPatch.TimeRecVec;
  % rhoFinal
  rhoFinalFilePatch.rho = rhoFinalFile2.rho;
  %move it
  path2dir = [ 'runOPfiles/' dirPatched ];
  movefile(paramFilePatchName, path2dir);
  movefile(runFilePatchName, path2dir);
  movefile(opFilePatchName, path2dir);
  movefile(rhoFinalFilePatchName, path2dir);
end % if proceed
end % function

% check parameters
function proceed = checkSameParameters( checkParamVec1, checkParamVec2 )
if any( checkParamVec1(1:3) ~= checkParamVec2(1:3) )
  fprintf('System size is not the same. Cannot proceed')
  proceed = 0;
elseif any( checkParamVec1(4:end) ~= checkParamVec2(4:end) )
  printf('l1 l2 vd bc \n');
  printf('run 1: ');
  printf( '%f', checkParamVec1);
  printf( '%f', checkParamVec2);
  printf('\n run 2: ');
  printf('\n');
  proceedStr = input('Key parameters not the same. Shall I proceed (y/n)?');
  if proceedStr == 'y'
    fprintf('Key parameters are not the same, but I patching\n')
    proceed = 1;
  else
    fprintf('Key parameters are not the same. Not patching\n')
    proceed = 0;
  end
else
  fprintf('Key parameters are the same. I am patching\n')
  proceed = 1;
end
end
