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

function [DenRecObj]  = ...
    HR2DrotDenEvolverFTBodyIdC(rho,ParamObj,...
    TimeObj,GridObj,DiffMobObj,Flags, lfid)

global RunSave

fprintf(lfid,'In body of code\n');

%Set N since it used so frequently
Nx  = ParamObj.Nx;
Ny  = ParamObj.Ny;
Nm  = ParamObj.Nm;
N3 = Nx*Ny*Nm;

% Declare dt since it's used so much
dt = TimeObj.dt;

% FT initial density and max density
TotalDensity = sum(sum(sum(rho)));
rho_FT = fftshift(fftn(rho));

%Initialize matrices that change size the +1 is to include initial density
if Flags.SaveMe == 1
  Density_rec       = zeros( Nx, Ny, Nm, TimeObj.N_recChunk );    % Store density amplitudes
  DensityFT_rec      = zeros( Nx, Ny, Nm, TimeObj.N_recChunk );   % Store k-space amplitdues
else
  Density_rec = 0;
  DensityFT_rec = 0;
end

% Recording indecies
jrectemp = 1; % Temporary holder for Density_rec
jrec     = 2; % Actual index for RunSave
jchunk   = 1; % Write chunk index

%Set up Diffusion operator, discrete k-space propagator, and interaction
%Set up Diffusion operator in cube form
[Lop] = DiffOpBuilderIsoDiffCube(DiffMobObj,GridObj,Nx,Ny,Nm);
Prop = exp(Lop .* dt);   % Exponentiate the elements

%%%%%%%%%%%%%%%%%%%Mayer function stuff%%%%%%%%%%%%%%%%%%%%%%%%%%
Fm_FT = fftshift(fftn( mayerFncHr(...
    Nx, Ny, Nm, ParamObj.Lx, ParamObj.Ly, ParamObj.L_rod) ) );

%Hard rod interactions
if Flags.Interactions
    GammaExCube_FT = dRhoIntCalcVcFtId(rho,rho_FT,Fm_FT,ParamObj,...
        GridObj,DiffMobObj);
else
    GammaExCube_FT = zeros(Nx,Ny,Nm);
end

%Driven Term
if Flags.Drive
  % Build the sin and cos phi once
  phi = zeros( 1, 1, Nm );
  phi(1,1,:) = GridObj.phi;
  cosPhi3 = cos( repmat( phi, [Nx, Ny, 1] ) );
  sinPhi3 = sin( repmat( phi, [Nx, Ny, 1] ) );
 
  GammaDrCube_FT  = ...
      dRhoDriveCalcFtId(rho,ParamObj.vD,...
      cosPhi3, sinPhi3,DiffMobObj.ikx3,DiffMobObj.iky3);
else
    GammaDrCube_FT = zeros(Nx,Ny,Nm);
end

%Total
GammaCube_FT = GammaDrCube_FT + GammaExCube_FT ;

% Take the first step- Euler. Element by element mulitplication
if( Flags.StepMeth == 0 ) % AB 1
  NlPf =  dt;
 [rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
elseif( Flags.StepMeth == 1 ) % AB 2 
  NlPf = 3 * dt / 2;
  NlPrevPf = dt / 2;
  [rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT,dt);
elseif( Flags.StepMeth == 2 ) % HAB1
  NlPf = dt .* Prop;
  [rho_FTnext] = DenStepperHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
elseif( Flags.StepMeth == 3 ) % HAB2
  NlPf = 3 * dt / 2 .* Prop;
  NlPrevPf = dt / 2 .* Prop .* Prop;
  [rho_FTnext] = DenStepperHAB1cPf( Prop, rho_FT, GammaCube_FT, dt .* Prop );
elseif( Flags.StepMeth == 4 ) % BHAB1
  NlPf = dt / 2 .* ( 1 + Prop);
  [rho_FTnext] = DenStepperBHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf);
elseif( Flags.StepMeth == 5 ) % BHAB2
  NlPf = dt / 2 * ( 2 + Prop );
  NlPrevPf = dt / 2;
  [rho_FTnext] = DenStepperBHAB1cPf( ...
  Prop, rho_FT, GammaCube_FT, dt / 2 .* ( 1 + Prop) );
elseif( Flags.StepMeth == 6 ) % Exponential Euler
  GamProp = ( Prop - 1 ) ./ Lop;
  [rho_FTnext] = DenStepperEEM1c( Prop, GamProp, rho_FT,GammaCube_FT);
else
  error('No stepping method selected');
end


tic
ShitIsFucked = 0;
SteadyState  = 0;
MaxReldRho   = 0; % Initialize this so things don't get messed up

fprintf(lfid,'Starting master time loop\n');
for t = 1:TimeObj.N_time-1
    %Save the previous and take one step forward.
    % Save the old drho
    GammaCube_FTprev = GammaCube_FT;
    rho_prev       = rho;
    
    %Need to update rho!!!
    rho_FT      = rho_FTnext;
    
    % Calculate rho if there is driving or interactions
    if Flags.Interactions || Flags.Drive
        rho    = real(ifftn(ifftshift(rho_FT)));
    end
    
    %Hard rod interactions
    if Flags.Interactions
        GammaExCube_FT = dRhoIntCalcVcFtId(rho,rho_FT,Fm_FT,ParamObj,...
            GridObj,DiffMobObj);
    end
    
    %Driven Term
    if Flags.Drive
      GammaDrCube_FT  = ...
        dRhoDriveCalcFtId(rho,ParamObj.vD,...
        cosPhi3, sinPhi3,DiffMobObj.ikx3,DiffMobObj.iky3);
    end
    
    GammaCube_FT = GammaDrCube_FT + GammaExCube_FT ;
    % Take a step
    if( Flags.StepMeth == 0 ) 
        [rho_FTnext] = DenStepperAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
    elseif( Flags.StepMeth == 1 )
        [rho_FTnext] = DenStepperAB2cPf( ...
             Prop, rho_FT,GammaCube_FT,GammaCube_FTprev,NlPf, NlPrevPf );
    elseif( Flags.StepMeth == 2 ) 
      [rho_FTnext] = DenStepperHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf);
    elseif( Flags.StepMeth == 3 )
        [rho_FTnext] = DenStepperHAB2cPf( ...
            Prop, rho_FT, GammaCube_FT,GammaCube_FTprev, NlPf, NlPrevPf );
    elseif( Flags.StepMeth == 4 )
      [rho_FTnext] = DenStepperBHAB1cPf( Prop, rho_FT, GammaCube_FT, NlPf );
    elseif( Flags.StepMeth == 5 )
      [rho_FTnext] = DenStepperBHAB2cPf( ...
        Prop, rho_FT, GammaCube_FT,GammaCube_FTprev, NlPf, NlPrevPf );
    elseif( Flags.StepMeth == 6 ) 
      [rho_FTnext] = DenStepperEEM1c( Prop, GamProp, rho_FT,GammaCube_FT);
    end

  %Save everything
  if ( mod(t,TimeObj.N_dtRec) == 0 )
    
    if Flags.SaveMe
      fprintf(lfid,'%f percent done\n',t./TimeObj.N_time*100);
      [SteadyState,ShitIsFucked,MaxReldRho] = ...
        BrokenSteadyDenTracker(rho,rho_prev, TotalDensity ,TimeObj);
      DensityFT_rec(:,:,:,jrectemp)   = rho_FT;
      Density_rec(:,:,:,jrectemp)     = rho;
    else
      [SteadyState,ShitIsFucked,MaxReldRho] = ...
        BrokenSteadyDenTracker(rho,rho_prev,TotalDensity ,TimeObj);
    end

      % Write a chunk to disk
    if ( mod(t, TimeObj.N_dtChunk ) == 0 )
      % Record Density_recs to file
      RecIndTemp = ...
        (jchunk-1) *  TimeObj.N_recChunk + 1 : jchunk * TimeObj.N_recChunk;
      % Shift by one because we include zero
      RecIndTemp = RecIndTemp + 1;
      RunSave.Den_rec(:,:,:,RecIndTemp) = Density_rec;
      RunSave.DenFT_rec(:,:,:,RecIndTemp) = DensityFT_rec;
      jrectemp = 0;
      jchunk = jchunk + 1;
    end
    jrectemp = jrectemp + 1;
    jrec = jrec + 1;

    % Break out if shit is fucked or done. Write first though
    if ShitIsFucked == 1 || SteadyState == 1;
      jrectemp = jrectemp - 1;
      StartInd = (jchunk-1) *  TimeObj.N_recChunk + 1;
      RecIndTemp = StartInd:StartInd + (jrectemp) - 1;
      % Shift by one because we include zero
      RecIndTemp = RecIndTemp + 1;
      RunSave.Den_rec(:,:,:,RecIndTemp) = Density_rec(:,:,:,1:jrectemp);
      RunSave.DenFT_rec(:,:,:,RecIndTemp) = DensityFT_rec(:,:,:,1:jrectemp);
      break 
    end
  end %end recording
  
end %end time loop
fprintf(lfid,'Finished master time loop\n');

% Update last rho
t =  t + 1;
rho_prev  = rho;
rho_FT    = rho_FTnext;
rho       = real(ifftn(ifftshift(rho_FT)));

%Save everything
if ( mod(t,TimeObj.N_dtRec)== 0 )

  if Flags.SaveMe
    fprintf(lfid,'%f percent done\n',t./TimeObj.N_time*100);
    [SteadyState,ShitIsFucked,MaxReldRho] = ...
      BrokenSteadyDenTracker(rho,rho_prev, TotalDensity ,TimeObj);
    DensityFT_rec(:,:,:,jrectemp)   = rho_FT;
    Density_rec(:,:,:,jrectemp)     = rho;
  else
    [SteadyState,ShitIsFucked,MaxReldRho] = ...
      BrokenSteadyDenTracker(rho,rho_prev,TotalDensity ,TimeObj);
  end
  
  if ( mod(t, TimeObj.N_dtChunk ) == 0 )
      % Record Density_recs to file
      RecIndTemp = ...
        (jchunk-1) *  TimeObj.N_recChunk + 1 : jchunk * TimeObj.N_recChunk;
      % Shift by one because we include zero
      RecIndTemp = RecIndTemp + 1;
      RunSave.Den_rec(:,:,:,RecIndTemp) = Density_rec;
      RunSave.DenFT_rec(:,:,:,RecIndTemp) = DensityFT_rec;
  end
  jrec = jrec + 1; % Still +1. Programs assumes this always happens
end %end recording

%If something broke, return zeros. Else, return the goods
if ShitIsFucked
    fprintf('Density is either negative or not conserved.\n');
    fprintf('I have done %i steps out of %i.\n',t, TimeObj.N_time);
    
elseif SteadyState
    fprintf('Things are going steady if you know what I mean.\n');
    fprintf('I have done %i steps out of %i.\n',t, TimeObj.N_time);
end

% Create vector of recorded times 
jrec = jrec - 1;
TimeRecVec    = (0:jrec-1) * TimeObj.t_rec;

trun = toc;

% Save useful info
DenRecObj.DidIBreak    = ShitIsFucked;
DenRecObj.SteadyState  = SteadyState;
DenRecObj.MaxReldRho   = MaxReldRho;
DenRecObj.TimeRecVec   = TimeRecVec;
DenRecObj.rhoFinal     = rho;
DenRecObj.runTime      = trun;

end %function
