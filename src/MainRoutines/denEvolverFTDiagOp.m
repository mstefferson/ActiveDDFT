% HR2DrotDenEvolverFTBodyIdC.M
%
% Description:
% Code is the main body  for the code that
% uses discrete Fourier transforms to solve the diffusion equation of rods in
% 2 spatial directions. These Rods are allowed to diffuse rotationally.
%
% Includes hard rod interactions, assuming in the interaction that the rods
% are infinitely thin.
%
%
%
% Density matrix is set up as rho(x,y,phi)-----> RHO(kx,ky,km)
%
% Isotropic diffusion. The propagaotr is a cube.
%
function [denRecObj,rho]  = ...
  denEvolverFTDiagOp(rho,systemObj,...
  timeObj,gridObj,diffObj,interObj,polarDrive,...
  noise, dRhoFlux, densityDepDiff, flags,lfid)
% globals
global runSave
% where you at
fprintf(lfid,'In body of code\n');
%Set N since it used so frequently
n1  = systemObj.n1;
n2  = systemObj.n2;
n3  = systemObj.n3;
% Declare dt since it's used so much
dt = timeObj.dt;
% FT initial density and max density
rho_FT = fftshift(fftn(rho));
if n3 == 1
  constConc = rho_FT( n1/2 + 1, n2/2 + 1);
else
  constConc = rho_FT( n1/2 + 1, n2/2 + 1, n3/2 + 1);
end
%Initialize matrices that change size the +1 is to include initial density
if flags.SaveMe == 1
  Density_rec = zeros( n1, n2, n3, timeObj.N_recChunk );    % Store density amplitudes
  DensityFT_rec = zeros( n1, n2, n3, timeObj.N_recChunk );   % Store k-space amplitdues
else
  Density_rec = 0;
  DensityFT_rec = 0;
end
% Recording indecies
jrectemp = 1; % Temporary holder for Density_rec
jrec     = timeObj.recStartInd; % Actual index for runSave
%Set up Diffusion operator, discrete k-space propagator, and interaction
[lop] = DiffOpBuilderIsoDiffCube(diffObj,gridObj,n1,n2,n3);
prop = exp(lop .* dt);   % Exponentiate the elements
% Interactions and driving
if flags.dRhoCalc
  rho    = real(ifftn(ifftshift(rho_FT)));
  % Calculate dRho from interactions and driving
  [GammaCube_FT, shitIsFucked, whatBroke1] = dRhoMaster( rho, rho_FT, ...
    interObj, systemObj, diffObj, polarDrive, noise, dRhoFlux, densityDepDiff );
else
  shitIsFucked = 0; shitIsFuckedTemp1 =0; shitIsFuckedTemp2 = 0;
  whatBroke1 = 0; whatBroke2 = 0; whatBroke3 = 0;
  GammaCube_FT = zeros( n1, n2, n3);
end
% Take the first step- Euler. Element by element mulitplication
if( flags.StepMeth == 0 ) % AB 1
  NlPf =  dt;
  [rho_FTnext] = DenStepperAB1cPf( prop, rho_FT, GammaCube_FT, NlPf );
elseif( flags.StepMeth == 1 ) % AB 2
  NlPf = 3 * dt / 2;
  NlPrevPf = dt / 2;
  [rho_FTnext] = DenStepperAB1cPf( prop, rho_FT, GammaCube_FT,dt);
  % Save prev Gamma if need be
  GammaCube_FTprev = GammaCube_FT;
elseif( flags.StepMeth == 2 ) % HAB1
  NlPf = dt .* prop;
  [rho_FTnext] = DenStepperHAB1cPf( prop, rho_FT, GammaCube_FT, NlPf );
elseif( flags.StepMeth == 3 ) % HAB2
  NlPf = 3 * dt / 2 .* prop;
  NlPrevPf = dt / 2 .* prop .* prop;
  [rho_FTnext] = DenStepperHAB1cPf( prop, rho_FT, GammaCube_FT, dt .* prop );
  % Save prev Gamma if need be
  GammaCube_FTprev = GammaCube_FT;
elseif( flags.StepMeth == 4 ) % BHAB1
  NlPf = dt / 2 .* ( 1 + prop);
  [rho_FTnext] = DenStepperBHAB1cPf( prop, rho_FT, GammaCube_FT, NlPf);
elseif( flags.StepMeth == 5 ) % BHAB2
  NlPf = dt / 2 * ( 2 + prop );
  NlPrevPf = dt / 2;
  [rho_FTnext] = DenStepperBHAB1cPf( ...
    prop, rho_FT, GammaCube_FT, dt / 2 .* ( 1 + prop) );
  % Save prev Gamma if need be
  GammaCube_FTprev = GammaCube_FT;
elseif( flags.StepMeth == 6 ) % Exponential Euler
  gamProp = ( prop - 1 ) ./ lop;
  gamProp( isnan( gamProp) ) = timeObj.dt;
%   gamProp( gridObj.k1ind0, gridObj.k2ind0, gridObj.k3ind0 ) = 0;
  [rho_FTnext] = DenStepperEEM1c( prop, gamProp, rho_FT,GammaCube_FT);
else
  error('No stepping method selected');
end
% Initialize some things and start time loop
tic
steadyState  = 0;
maxDrho   = 0; % Initialize this so things don't get messed up
whatBroke2 = [];
fprintf(lfid,'Starting master time loop\n');
if shitIsFucked == 0
  for t = 1:timeObj.N_time-1
    %Need to update rho!!!
    rho_FT = rho_FTnext;
    rhoPrev = rho;
    % Calculate rho if there is driving or interactions and steady
    rho    = real(ifftn(ifftshift(rho_FT)));
    if flags.dRhoCalc
      % Calculate dRho from interactions and driving
      [GammaCube_FT,shitIsFuckedTemp1, whatBroke1] = ...
        dRhoMaster( rho, rho_FT, interObj, systemObj,...
        diffObj, polarDrive, noise, dRhoFlux, densityDepDiff );
    end
    % Take a step
    if( flags.StepMeth == 0 )
      [rho_FTnext] = DenStepperAB1cPf( prop, rho_FT, GammaCube_FT, NlPf );
    elseif( flags.StepMeth == 1 )
      [rho_FTnext] = DenStepperAB2cPf( ...
        prop, rho_FT,GammaCube_FT,GammaCube_FTprev,NlPf, NlPrevPf );
      GammaCube_FTprev = GammaCube_FT;
    elseif( flags.StepMeth == 2 )
      [rho_FTnext] = DenStepperHAB1cPf( prop, rho_FT, GammaCube_FT, NlPf);
    elseif( flags.StepMeth == 3 )
      [rho_FTnext] = DenStepperHAB2cPf( ...
        prop, rho_FT, GammaCube_FT,GammaCube_FTprev, NlPf, NlPrevPf );
      GammaCube_FTprev = GammaCube_FT;
    elseif( flags.StepMeth == 4 )
      [rho_FTnext] = DenStepperBHAB1cPf( prop, rho_FT, GammaCube_FT, NlPf );
    elseif( flags.StepMeth == 5 )
      [rho_FTnext] = DenStepperBHAB2cPf( ...
        prop, rho_FT, GammaCube_FT,GammaCube_FTprev, NlPf, NlPrevPf );
      GammaCube_FTprev = GammaCube_FT;
    elseif( flags.StepMeth == 6 )
      [rho_FTnext] = DenStepperEEM1c( prop, gamProp, rho_FT,GammaCube_FT);
    end
    %Save everything
    if ( mod(t,timeObj.N_dtRec) == 0 )
      [steadyState,shitIsFuckedTemp2, whatBroke2, maxDrho] = ...
        BrokenSteadyDenTracker(rho, rhoPrev, rho_FT, constConc, timeObj, systemObj);
      if flags.SaveMe
        fprintf(lfid,'%f percent done\n',t./timeObj.N_time*100);
        DensityFT_rec(:,:,:,jrectemp)   = rho_FT;
        Density_rec(:,:,:,jrectemp)     = rho;
      end
      shitIsFucked = shitIsFuckedTemp1 + shitIsFuckedTemp2;
      % Write a chunk to disk
      if flags.SaveMe
        if ( mod(t, timeObj.N_dtChunk ) == 0 )
          % Record Density_recs to file
          jrecEnd = jrec+timeObj.N_recChunk-1;
          recIndTemp = jrec : jrecEnd;
          runSave.Den_rec(:,:,:,recIndTemp) = Density_rec;
          runSave.DenFT_rec(:,:,:,recIndTemp) = DensityFT_rec;
          runSave.numSavedRhos = recIndTemp(end);
          jrectemp = 0;
          jrec = jrecEnd + 1;
        end
      end
      jrectemp = jrectemp + 1;
      % Break out if shit is fucked or done. Write first though
      if shitIsFucked > 0 || steadyState == 1
        if flags.SaveMe
          % update record holders
          jrectemp = jrectemp - 1;
          jrecEnd = jrec + (jrectemp) - 1;
          recIndTemp = jrec:jrecEnd;
          % Save what remains
          if ~isempty(recIndTemp)
            runSave.Den_rec(:,:,:,recIndTemp) = Density_rec(:,:,:,1:jrectemp);
            runSave.DenFT_rec(:,:,:,recIndTemp) = DensityFT_rec(:,:,:,1:jrectemp);
            runSave.numSavedRhos = recIndTemp(end);
            jrec = jrecEnd + 1;
          end
        end
        break
      end
    end %end recording
  end %end time loop
end % broken from start
% Say you're done
fprintf(lfid,'Finished master time loop\n');
% Update last rho
t =  t + 1;
rho_FT    = rho_FTnext;
rho        = real(ifftn(ifftshift(rho_FT)));
trun = toc;
%Save everything
if flags.SaveMe
  if ( mod(t,timeObj.N_dtRec)== 0 )
    fprintf(lfid,'%f percent done\n',t./timeObj.N_time*100);
    % Save it is there is unsaved data
    if jrectemp > 0
      if flags.dRhoCalc
        rho    = real(ifftn(ifftshift(rho_FT)));
      end
      DensityFT_rec(:,:,:,jrectemp)   = rho_FT;
      Density_rec(:,:,:,jrectemp)     = rho;
      % Record Density_recs to file
      if ( mod(t, timeObj.N_dtChunk ) == 0 )
        jrecEnd = jrec+timeObj.N_recChunk-1;
        recIndTemp = jrec : jrecEnd;
        runSave.Den_rec(:,:,:,recIndTemp) = Density_rec;
        runSave.DenFT_rec(:,:,:,recIndTemp) = DensityFT_rec;
        runSave.numSavedRhos = recIndTemp(end);
      end
      jrec = jrecEnd + 1; % Still +1. Programs assumes this always happens
    end
  end
end %end recording
% Create vector of recorded times
jrec = jrec - 1;
TimeRecVec    = (0:jrec-1) * timeObj.t_rec;
%If something broke, return zeros. Else, return the goods
if shitIsFucked
  fprintf('Density is either negative or not conserved.\n');
  fprintf('I have done %i steps out of %i.\n',t, timeObj.N_time);
elseif steadyState
  fprintf('Things are going steady if you know what I mean.\n');
  fprintf('I have done %i steps out of %i.\n',t, timeObj.N_time);
else
  fprintf('Ran for the full time %.1f with %d jrec\n',timeObj.t_tot,jrec);
end
% Save useful info
denRecObj.didIrun      = 1;
denRecObj.DidIBreak    = shitIsFucked;
denRecObj.whatBroke    = [whatBroke1 whatBroke2];
denRecObj.SteadyState  = steadyState;
denRecObj.maxDrho      = maxDrho;
denRecObj.TimeRecVec   = TimeRecVec;
denRecObj.runTime      = trun;
denRecObj.simTime      = TimeRecVec(end);
end %function
