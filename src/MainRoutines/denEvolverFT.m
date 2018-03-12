% HR2DrotDenEvolverFTBody.m
%
% Description:
% Code is the main body  for the code that
% uses discrete Fourier transforms to solve the diffusion equation of rods in
% 2 spatial directions. These Rods are allowed to diffuse rotationally.
%
% Includes hard rod interactions, assuming in the interaction that the rods
% are infinitely thin.
%
% Anisotropic diffusion

% Density matrix is set up as rho(x,y,phi)-----> RHO(kx,ky,km)
%
% The Lopagator now includes all the terms. The all the cubes are turned
% into: a (n1*n2*n3) x (n1*n2*n3) linear operator and N^3 density vector
%
% Everything is sparsified
%
% Program never actually calculates the Lopagator, but uses expv from
% ExpoKit. Way Faster.
%
function [denRecObj, rho] = denEvolverFT(...
  rho, systemObj, timeObj, gridObj, diffObj, interObj, ...
  polarDrive, noise, dRhoFlux, flags,lfid )
% global
global runSave
% where you at
fprintf(lfid,'In body of code\n');
%Set N since it used so frequently
n1  = systemObj.n1;
n2  = systemObj.n2;
n3  = systemObj.n3;
N3 = n1*n2*n3;
N2 = n1*n2;
% Declare dt since it's used so much
dt = timeObj.dt;
% FT initial density and max density
rho_FT = fftshift(fftn(rho));
rhoVec_FT = reshape(rho_FT,N3,1);
if n3 == 1
  constConc = rho_FT( n1/2 + 1, n2/2 + 1);
else
  constConc = rho_FT( n1/2 + 1, n2/2 + 1, n3/2 + 1);
end
%Initialize matrices that change size the +1 is to include initial density
if flags.SaveMe == 1
  Density_rec       = zeros( n1, n2, n3, timeObj.N_recChunk );    % Store density amplitudes
  DensityFT_rec      = zeros( n1, n2, n3, timeObj.N_recChunk );   % Store k-space amplitdues
else
  Density_rec = 0;
  DensityFT_rec = 0;
end
% Recording indecies
jrectemp = 1; % Temporary holder for Density_rec
jrec     = timeObj.recStartInd; % Actual index for runSave
%Set up Diffusion operator, discrete k-space Lopagator, and interaction
[Lop] = DiffOpBuilderDr(diffObj,gridObj,n1,n2,n3,N2,N3);
% trig functions for driving not used here, but set it to zero
trigFnc.cosPhi3 = 0;
trigFnc.sinPhi3 = 0;
%Interactions
if flags.dRhoCalc
  [GammaEx_FT, shitIsFucked, whatBroke1] = dRhoMaster( rho, rho_FT,...     
    interObj, systemObj, diffObj, polarDrive, noise, dRhoFlux );
  GammaExVec_FT  = reshape( GammaEx_FT, N3,1);
else
  shitIsFucked = 0; shitIsFuckedTemp1 =0; shitIsFuckedTemp2 = 0;
  whatBroke1 = 0; whatBroke2 = 0; whatBroke3 = 0;
  GammaExVec_FT = zeros(N3,1);
end
% Take first step- Euler
if( flags.StepMeth == 0 ) % AB 1
  NlPf =  dt;
  [rhoVec_FTnext, ticExpInt] = DenStepperAB1Pf(...
    Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt );
elseif( flags.StepMeth == 1 ) % AB 2
  NlPf = 3 * dt / 2;
  NlPrevPf = dt / 2;
  [rhoVec_FTnext, ticExpInt] = DenStepperAB1Pf( ...
    Lop, rhoVec_FT, GammaExVec_FT, dt, dt  );
  % Save prev Gamma if need be
  GammaExVec_FTprev = GammaExVec_FT;
elseif( flags.StepMeth == 2 )  % HAB 1
  NlPf = dt;
  [rhoVec_FTnext, ticExpInt] = DenStepperHAB1Pf( ...
    Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt );
elseif( flags.StepMeth == 3 ) % HAB 2
  NlPf = 3 * dt / 2 ;
  NlPrevPf = dt / 2 ;
  [rhoVec_FTnext, ticExpInt] = DenStepperHAB1Pf( ...
    Lop, rhoVec_FT, GammaExVec_FT, dt, dt  );
  % Save prev Gamma if need be
  GammaExVec_FTprev = GammaExVec_FT;
elseif( flags.StepMeth == 4 ) % BHAB 1
  NlPf = dt / 2;
  [rhoVec_FTnext, ticExpInt] = DenStepperBHAB1Pf( ...
    Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt  );
elseif( flags.StepMeth == 5 ) % BHAB 2
  NlPf = dt ;
  NlPrevPf = dt / 2;
  NlExpPf =  dt / 2;
  [rhoVec_FTnext, ticExpInt] = DenStepperBHAB1Pf( ...
    Lop, rhoVec_FT, GammaExVec_FT, dt / 2, dt );
  % Save prev Gamma if need be
  GammaExVec_FTprev = GammaExVec_FT;
elseif( flags.StepMeth == 6 ) % phiV
  [rhoVec_FTnext, ticExpInt] = DenStepperPhiV( ...
    Lop, rhoVec_FT, GammaExVec_FT, dt );
else
  error('No stepping method selected');
end
% Initialize some things and start time loop
tic
steadyState  = 0;
maxDrho   = 0; % Initialize this so things don't get messed up
whatBroke2 = [];
whatBroke3 = [];
fprintf(lfid,'Starting master time loop\n');
if shitIsFucked == 0 
  for t = 1:timeObj.N_time-1
    %Need to update rho!!!
    rhoVec_FT = rhoVec_FTnext;
    rhoPrev = rho;
    % Calculate rho if there is driving or interactions
    rho_FT = reshape(rhoVec_FT,n1,n2,n3);
    rho = real(ifftn(ifftshift(rho_FT)));
    if flags.dRhoCalc
      [GammaEx_FT, shitIsFuckedTemp1, whatBroke1] = ...
        dRhoMaster( rho, rho_FT, ...
        interObj, systemObj, diffObj, polarDrive, noise, dRhoFlux );
      GammaExVec_FT  = reshape( GammaEx_FT, N3,1);
    end
    % Take step
    if( flags.StepMeth == 0 )
      [rhoVec_FTnext, ticExptemp] = DenStepperAB1Pf( ...
        Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt  );
    elseif( flags.StepMeth == 1 )
      [rhoVec_FTnext, ticExptemp] = DenStepperAB2Pf( ...
        Lop, rhoVec_FT, GammaExVec_FT, GammaExVec_FTprev, NlPf, NlPrevPf, dt  );
      GammaExVec_FTprev = GammaExVec_FT;
    elseif( flags.StepMeth == 2 )
      [rhoVec_FTnext, ticExptemp] = DenStepperHAB1Pf( ...
        Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt );
    elseif( flags.StepMeth == 3 )
      [rhoVec_FTnext, ticExptemp] = DenStepperHAB2Pf( ...
        Lop, rhoVec_FT, GammaExVec_FT, GammaExVec_FTprev, NlPf, NlPrevPf, dt  );
      GammaExVec_FTprev = GammaExVec_FT;
    elseif( flags.StepMeth == 4 )
      [rhoVec_FTnext, ticExptemp] = DenStepperBHAB1Pf( ...
        Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt  );
    elseif( flags.StepMeth == 5 )
      [rhoVec_FTnext, ticExptemp] = DenStepperBHAB2Pf( ...
        Lop, rhoVec_FT, GammaExVec_FT,GammaExVec_FTprev, NlPf, NlExpPf, NlPrevPf, dt );
      GammaExVec_FTprev = GammaExVec_FT;
    elseif( flags.StepMeth == 6 )
      [rhoVec_FTnext, ticExptemp] = DenStepperPhiV( ...
        Lop, rhoVec_FT, GammaExVec_FT, dt );
    end
    %Save everything
    if ( mod(t,timeObj.N_dtRec) == 0 )
      % Turn it to a cube if it hasn't been yet
      if flags.dRhoCalc == 0
        rho_FT = reshape(rhoVec_FT,n1,n2,n3);
        rho    = real(ifftn(ifftshift(rho_FT)));
      end
      rho_FTnext = reshape(rhoVec_FTnext,n1,n2,n3);
      [steadyState,shitIsFuckedTemp2, whatBroke2, maxDrho] = ...
        BrokenSteadyDenTracker(rho, rhoPrev, rho_FT, constConc, timeObj, systemObj);
      if flags.SaveMe
        fprintf(lfid,'%f percent done\n',t./timeObj.N_time*100);
        DensityFT_rec(:,:,:,jrectemp)   = rho_FT;
        Density_rec(:,:,:,jrectemp)     = rho;
      end
      % Make sure neither dRho or Broken den broke
      shitIsFucked = shitIsFuckedTemp1 + shitIsFuckedTemp2;
      %Make sure things are taking too long. This is a sign density---> inf
      [TooLong] = ExpTooLongChecker(...
        ticExptemp,ticExpInt,rhoVec_FT,n1,n2,n3,jrec);
      if TooLong; shitIsFucked = 1; whatBroke3='Too long exp'; end
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
            jrec = jrecEnd + 1;
          end
        end
        break
      end
    end %end recording
  end %end time loop
end % dont start if already broken
% Say you're done
fprintf(lfid,'Finished master time loop\n');
% Update last rho
t =  t + 1;
rhoVec_FT  = rhoVec_FTnext;
rho_FT     = reshape(rhoVec_FT,n1,n2,n3);
rho        = real(ifftn(ifftshift(rho_FT)));
trun = toc;
%Save everything
if flags.SaveMe
  if ( mod(t,timeObj.N_dtRec)== 0 )
    fprintf(lfid,'%f percent done\n',t./timeObj.N_time*100);
    % Save it is there is unsaved data
    if jrectemp > 0
      if flags.dRhoCalc
        rho_FT = reshape(rhoVec_FT,n1,n2,n3);
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
  fprintf('I have done %i steps out of %i with %d jrecs\n',t, timeObj.N_time,jrec);
elseif steadyState
  fprintf('Things are going steady if you know what I mean.\n');
  fprintf('I have done %i steps out of %i with %d jrecs\n',t, timeObj.N_time,jrec);
else
  fprintf('Ran for the full time %.1f with %d jrecs\n',timeObj.t_tot,jrec);
end
% Save useful info
denRecObj.didIrun      = 1;
denRecObj.DidIBreak    = shitIsFucked;
denRecObj.whatBroke    = [whatBroke1 whatBroke2 whatBroke3];
denRecObj.SteadyState  = steadyState;
denRecObj.maxDrho      = maxDrho;
denRecObj.TimeRecVec   = TimeRecVec;
denRecObj.runTime      = trun;
denRecObj.simTime      = TimeRecVec(end);
end %function
