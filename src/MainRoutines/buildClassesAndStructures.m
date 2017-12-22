function [timeObj, gridObj, diffObj, interObj, ...
  noise, polarDrive, densityDepDiff, runTime] = buildClassesAndStructures( ...
  systemObj, particleObj, flags, timeObj, lfid )
% set up the time object
dtOrig = timeObj.dt ;
if timeObj.scaleDt && particleObj.fD ~=0
  timeObj.dt = min( timeObj.dt, timeObj.dt * 30  / particleObj.fD );
end
% Fix the time
ss_epsilon = timeObj.ss_epsilon;
scaleDt = timeObj.scaleDt;
[timeObj]= ...
  TimeStepRecMaker(timeObj.dt,timeObj.t_tot,...
  timeObj.t_rec,timeObj.t_write);
% Scale ss_epsilon by delta_t. Equivalent to checking d rho /dt has reached
% steady state instead of d rho
timeObj.ss_epsilon = ss_epsilon;
timeObj.ss_epsilon_dt = ss_epsilon .* timeObj.dt;
timeObj.scaleDt = scaleDt;
timeObj.dt_orig = dtOrig;
timeObj.recStartInd = 2; % start at 2 since t = 0 is ind = 1
% Make all the grid stuff %
tGridID = tic;
[gridObj] = GridMakerPBCxk(systemObj.n1,systemObj.n2,systemObj.n3,...
  systemObj.l1,systemObj.l2,systemObj.l3);
gridrunTime = toc(tGridID);
if flags.Verbose
  fprintf('Made grid t%d_%d: %.3g \n', ...
    runObj.trialID, runObj.runID, gridrunTime);
end
fprintf(lfid,'Made Grid: %.3g \n', gridrunTime);
runTime.grid = gridrunTime;
%Make diffusion coeff
tDiffID = tic;
[diffObj] =  DiffMobCoupCoeffCalc( systemObj.tmp,...
  particleObj.mob,particleObj.mobPar,particleObj.mobPerp,particleObj.mobRot,...
  gridObj.k1, gridObj.k2, gridObj.k3, ...
  gridObj.k1rep2, gridObj.k2rep2,particleObj.fD);
diffRunTime = toc(tDiffID);
if flags.Verbose
  fprintf('Made diffusion object t%d_%d: %.3g\n', ...
    runObj.trialID, runObj.runID, diffRunTime);
end
fprintf(lfid,'Made diffusion object: %.3g\n', diffRunTime);
runTime.diff = diffRunTime;
% Set-up interactions and external potentials
[interObj] =  interObjMaker( particleObj, systemObj, gridObj );
% set-up noise
[noise] = DRhoNoiseClass( systemObj.noise, ...
  systemObj.n1, systemObj.n2, systemObj.n3, systemObj.l1, systemObj.l2, ...
  systemObj.l3, diffObj.D_pos, diffObj.D_rot, diffObj.ik1rep3, diffObj.ik2rep3, ...
  diffObj.ik3rep3, timeObj.dt);
% set-up driving
dRhoDriveFlag = flags.Drive && flags.DiagLop;
polarDrive  = DrhoPolarDriveClass( dRhoDriveFlag, particleObj.fD, systemObj.n1,...
  systemObj.n2, systemObj.n3, gridObj.x3, gridObj.k1rep2, gridObj.k2rep2 );
% build density dep diffusion class
densityDepDiff = densityDepDiffClassHandler( ...
  particleObj, systemObj, gridObj, diffObj, []  );

