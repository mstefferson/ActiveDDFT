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
function [denRecObj]  = ...
  denEvolverFTDiagOp(rho,systemObj,particleObj,...
  timeObj,gridObj,diffObj,interObj,flags,lfid)
global runSave
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
  constConc = rho_FT( n1/2 + 1, n1/2 + 1);
else
  constConc = rho_FT( n1/2 + 1, n1/2 + 1, n3/2 + 1);
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
jrec     = 2; % Actual index for runSave
jchunk   = 1; % Write chunk index
%Set up Diffusion operator, discrete k-space propagator, and interaction
[Lop] = DiffOpBuilderIsoDiffCube(diffObj,gridObj,n1,n2,n3);
Prop = exp(Lop .* dt);   % Exponentiate the elements
%Driven Term
if flags.Drive
  % Build the sin and cos phi once
  phi = zeros( 1, 1, n3 );
  phi(1,1,:) = gridObj.x3;
  cosPhi3 = cos( repmat( phi, [n1, n2, 1] ) );
  sinPhi3 = sin( repmat( phi, [n1, n2, 1] ) );
else
  cosPhi3 = 0;
  sinPhi3 = 0;
end
% Interactions and driving
if interObj.anyInter || flags.Drive
  rho    = real(ifftn(ifftshift(rho_FT)));
  % Calculate dRho from interactions and driving
  [GammaCube_FT, shitIsFucked, whatBroke1] = dRhoMaster( rho, rho_FT, flags, ...
    interObj, systemObj, diffObj, particleObj, cosPhi3, sinPhi3 );
else
  GammaCube_FT = zeros( n1, n2, n3);
end
% Take the first step- Euler. Element by element mulitplication
if( flags.StepMeth == 0 ) % AB 1
  NlPf =  dt;
  [rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
elseif( flags.StepMeth == 1 ) % AB 2
  NlPf = 3 * dt / 2;
  NlPrevPf = dt / 2;
  [rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT,dt);
  % Save prev Gamma if need be
  GammaCube_FTprev = GammaCube_FT;
elseif( flags.StepMeth == 2 ) % HAB1
  NlPf = dt .* Prop;
  [rho_FTnext] = DenStepperHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
elseif( flags.StepMeth == 3 ) % HAB2
  NlPf = 3 * dt / 2 .* Prop;
  NlPrevPf = dt / 2 .* Prop .* Prop;
  [rho_FTnext] = DenStepperHAB1cPf( Prop, rho_FT, GammaCube_FT, dt .* Prop );
  % Save prev Gamma if need be
  GammaCube_FTprev = GammaCube_FT;
elseif( flags.StepMeth == 4 ) % BHAB1
  NlPf = dt / 2 .* ( 1 + Prop);
  [rho_FTnext] = DenStepperBHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf);
elseif( flags.StepMeth == 5 ) % BHAB2
  NlPf = dt / 2 * ( 2 + Prop );
  NlPrevPf = dt / 2;
  [rho_FTnext] = DenStepperBHAB1cPf( ...
    Prop, rho_FT, GammaCube_FT, dt / 2 .* ( 1 + Prop) );
  % Save prev Gamma if need be
  GammaCube_FTprev = GammaCube_FT;
elseif( flags.StepMeth == 6 ) % Exponential Euler
  GamProp = ( Prop - 1 ) ./ Lop;
  GamProp( gridObj.k1ind0, gridObj.k2ind0, gridObj.k3ind0 ) = 0;
  [rho_FTnext] = DenStepperEEM1c( Prop, GamProp, rho_FT,GammaCube_FT);
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
    % Calculate rho if there is driving or interactions
    if interObj.anyInter || flags.Drive
      rho    = real(ifftn(ifftshift(rho_FT)));
      % Calculate dRho from interactions and driving
  %     keyboard
      [GammaCube_FT,shitIsFuckedTemp1, whatBroke1] = dRhoMaster( rho, rho_FT, flags,...
        interObj, systemObj, diffObj, particleObj, cosPhi3, sinPhi3 );
    end
    % Take a step
    if( flags.StepMeth == 0 )
      [rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
    elseif( flags.StepMeth == 1 )
      [rho_FTnext] = DenStepperAB2cPf( ...
        Prop, rho_FT,GammaCube_FT,GammaCube_FTprev,NlPf, NlPrevPf );
      GammaCube_FTprev = GammaCube_FT;
    elseif( flags.StepMeth == 2 )
      [rho_FTnext] = DenStepperHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf);
    elseif( flags.StepMeth == 3 )
      [rho_FTnext] = DenStepperHAB2cPf( ...
        Prop, rho_FT, GammaCube_FT,GammaCube_FTprev, NlPf, NlPrevPf );
      GammaCube_FTprev = GammaCube_FT;
    elseif( flags.StepMeth == 4 )
      [rho_FTnext] = DenStepperBHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
    elseif( flags.StepMeth == 5 )
      [rho_FTnext] = DenStepperBHAB2cPf( ...
        Prop, rho_FT, GammaCube_FT,GammaCube_FTprev, NlPf, NlPrevPf );
      GammaCube_FTprev = GammaCube_FT;
    elseif( flags.StepMeth == 6 )
      [rho_FTnext] = DenStepperEEM1c( Prop, GamProp, rho_FT,GammaCube_FT);
    end
    %Save everything
    if ( mod(t,timeObj.N_dtRec) == 0 )
      if flags.Interactions == 0 && flags.Drive == 0
        rho    = real(ifftn(ifftshift(rho_FT)));
      end
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
          RecIndTemp = ...
            (jchunk-1) *  timeObj.N_recChunk + 1 : jchunk * timeObj.N_recChunk;
          % Shift by one because we include zero
          RecIndTemp = RecIndTemp + 1;
          runSave.Den_rec(:,:,:,RecIndTemp) = Density_rec;
          runSave.DenFT_rec(:,:,:,RecIndTemp) = DensityFT_rec;
          jrectemp = 0;
          jchunk = jchunk + 1;
        end
      end
      jrectemp = jrectemp + 1;
      jrec = jrec + 1;
      % Break out if shit is fucked or done. Write first though
      if shitIsFucked > 0 || steadyState == 1
        if flags.SaveMe
          jrectemp = jrectemp - 1;
          StartInd = (jchunk-1) *  timeObj.N_recChunk + 1;
          RecIndTemp = StartInd:StartInd + (jrectemp) - 1;
          % Shift by one because we include zero
          RecIndTemp = RecIndTemp + 1;
          % Save what remains
          if ~isempty(RecIndTemp)
            runSave.Den_rec(:,:,:,RecIndTemp) = Density_rec(:,:,:,1:jrectemp);
            runSave.DenFT_rec(:,:,:,RecIndTemp) = DensityFT_rec(:,:,:,1:jrectemp);
          end
        end
        break
      end
    end %end recording
  end %end time loop
end % broken from start
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
    % Turn it to a cube if it hasn't been yet
    if flags.Interactions == 0 && flags.Drive == 0
      rho    = real(ifftn(ifftshift(rho_FT)));
    end
    DensityFT_rec(:,:,:,jrectemp)   = rho_FT;
    Density_rec(:,:,:,jrectemp)     = rho;
    % Record Density_recs to file
    if ( mod(t, timeObj.N_dtChunk ) == 0 )
      RecIndTemp = ...
        (jchunk-1) *  timeObj.N_recChunk + 1 : jchunk * timeObj.N_recChunk;
      % Shift by one because we include zero
      RecIndTemp = RecIndTemp + 1;
      runSave.Den_rec(:,:,:,RecIndTemp) = Density_rec;
      runSave.DenFT_rec(:,:,:,RecIndTemp) = DensityFT_rec;
    end
    jrec = jrec + 1; % Still +1. Programs assumes this always happens
  end
end %end recording
%If something broke, return zeros. Else, return the goods
if shitIsFucked
  fprintf('Density is either negative or not conserved.\n');
  fprintf('I have done %i steps out of %i.\n',t, timeObj.N_time);
elseif steadyState
  fprintf('Things are going steady if you know what I mean.\n');
  fprintf('I have done %i steps out of %i.\n',t, timeObj.N_time);
else
  fprintf('Ran for the full time %.1f\n',timeObj.t_tot);
end
% Create vector of recorded times
jrec = jrec - 1;
TimeRecVec    = (0:jrec-1) * timeObj.t_rec;
% Save useful info
denRecObj.didIrun      = 1;
denRecObj.DidIBreak    = shitIsFucked;
denRecObj.whatBroke    = [whatBroke1 whatBroke2];
denRecObj.SteadyState  = steadyState;
denRecObj.maxDrho      = maxDrho;
denRecObj.TimeRecVec   = TimeRecVec;
denRecObj.rhoFinal     = rho;
denRecObj.runTime      = trun;
denRecObj.simTime      = TimeRecVec(end);
end %function
