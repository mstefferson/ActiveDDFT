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
% into: a (Nx*Ny*Nm) x (Nx*Ny*Nm) linear operator and N^3 density vector
%
% Everything is sparsified
%
% Program never actually calculates the Lopagator, but uses expv from
% ExpoKit. Way Faster.
%
% Interactions handled using Mayer function.


function [DenRecObj] = HR2DrotDenEvolverFTBody(...
    wfid,lfid,rho,ParamObj, TimeObj,GridObj,DiffMobObj,feq)

fprintf(lfid,'In body of code\n');
% Create a text file that tells user what percent of the program has
% finished

PrgmTrackNm = sprintf('PrgmTracker_%i.txt',ParamObj.trial);
tfid        = fopen(PrgmTrackNm,'w');

%Set N since it used so frequently
Nx  = ParamObj.Nx;
Ny  = ParamObj.Ny;
Nm  = ParamObj.Nm;
N3 = Nx*Ny*Nm;
N2 = Nx*Ny;

% Declare dt since it's used so much
dt = TimeObj.delta_t;

% FT initial density and max density
TotalDensity = sum(sum(sum(rho)));
rho_FT = fftshift(fftn(rho));
rhoVec_FT = reshape(rho_FT,N3,1);

global Density_rec
global DensityFT_rec
%Initialize matrices that change size the +1 is to include initial density
if ParamObj.SaveMe == 1
    Density_rec       = zeros( Nx, Ny, Nm, TimeObj.N_record + 1 );      % Store density amplitudes
    DensityFT_rec      = zeros( Nx, Ny, Nm, TimeObj.N_record + 1 );      % Store k-space amplitdues
    
    %Initialize records
    DensityFT_rec(:,:,:,1)   = rho_FT;
    Density_rec(:,:,:,1)     = rho;
else
    Density_rec = 0;
    DensityFT_rec = 0;
end

j_record = 2;     %Record holder

%Set up Diffusion operator, discrete k-space Lopagator, and interaction
[Lop] = DiffOpBuilderDr(DiffMobObj,GridObj,Nm,N2,N3);

%%%%%%%%%%%%%%%%%%%Mayer function stuff%%%%%%%%%%%%%%%%%%%%%%%%%%
Fm_FT = fftshift(fftn( MayerFncDiffBtwPntsCalc(...
    Nx, Ny, Nm, ParamObj.Lx, ParamObj.Ly, ParamObj.L_rod) ));

%Hard rod interactions
if ParamObj.Interactions
    GammaExVec_FT  = reshape( ...
        dRhoIntCalcVcFt( rho,rho_FT,Fm_FT,ParamObj,GridObj,DiffMobObj), ...
        N3,1);
else
    GammaExVec_FT = zeros(N3,1);
    fprintf(wfid,'Interacts are off my Lord\n');
end

% Take first step- Euler
if( ParamObj.StepMeth == 0 ) % AB 1

  NlPf =  dt;
 [rhoVec_FTnext, ticExpInt] = DenStepperAB1Pf(...
  Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt );

elseif( ParamObj.StepMeth == 1 ) % AB 2

  NlPf = 3 * dt / 2;
  NlPrevPf = dt / 2;
  [rhoVec_FTnext, ticExpInt] = DenStepperAB1Pf( ...
   Lop, rhoVec_FT, GammaExVec_FT, dt, dt  );

elseif( ParamObj.StepMeth == 2 )  % HAB 1
  
  NlPf = dt;
  [rhoVec_FTnext, ticExpInt] = DenStepperHAB1Pf( ...
   Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt );

elseif( ParamObj.StepMeth == 3 ) % HAB 2

  NlPf = 3 * dt / 2 ;
  NlPrevPf = dt / 2 ;
  [rhoVec_FTnext, ticExpInt] = DenStepperHAB1Pf( ...
   Lop, rhoVec_FT, GammaExVec_FT, dt, dt  );

elseif( ParamObj.StepMeth == 4 ) % BHAB 1
  
  NlPf = dt / 2;
  [rhoVec_FTnext, ticExpInt] = DenStepperBHAB1Pf( ...
   Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt  );

elseif( ParamObj.StepMeth == 5 ) % BHAB 2
  
  NlPf = dt ;
  NlPrevPf = dt / 2;
  NlExpPf =  dt / 2;
  
  [rhoVec_FTnext, ticExpInt] = DenStepperBHAB1Pf( ...
   Lop, rhoVec_FT, GammaExVec_FT, dt / 2, dt );

elseif( ParamObj.StepMeth == 6 ) % phiV

  [rhoVec_FTnext, ticExpInt] = DenStepperPhiV( ...
   Lop, rhoVec_FT, GammaExVec_FT, dt );

else
  fprintf('No stepping method selected');
end

tic
ShitIsFucked = 0;
SteadyState  = 0;
MaxReldRho   = 0; % Initialize this so things don't get messed up

% keyboard
fprintf(lfid,'Starting master time loop\n');
for t = 1:TimeObj.N_time-1
    %Save the previous and take one step forward.
    % Save the old drho
    GammaExVec_FTprev = GammaExVec_FT;
    rhoVec_FTprev  = rhoVec_FT;
    
    %Need to update rho!!!
    rhoVec_FT      = rhoVec_FTnext;
    
    % Calculate rho if there is driving or interactions
    if ParamObj.Interactions || ParamObj.Drive
        rho_FT = reshape(rhoVec_FT,Nx,Ny,Nm);
        rho    = real(ifftn(ifftshift(rho_FT)));
    end
    
    %Hard rod interactions
    if ParamObj.Interactions == true
        GammaExVec_FT  = reshape( ...
            dRhoIntCalcVcFt( ...
            rho,rho_FT,Fm_FT,ParamObj,GridObj,DiffMobObj),...
            N3,1);
    end
       
    % Take step
    if( ParamObj.StepMeth == 0 ) 
%      [rhoVec_FTnext, ticExptemp] = DenStepperAB1( ...
%      Lop, rhoVec_FT, GammaExVec_FT, dt );
     [rhoVec_FTnext, ticExptemp] = DenStepperAB1Pf( ...
     Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt  );
    elseif( ParamObj.StepMeth == 1 )
       %[rhoVec_FTnext, ticExptemp] = DenStepperAB2(... 
          %Lop, rhoVec_FT, GammaExVec_FT, GammaExVec_FTprev, dt );
      [rhoVec_FTnext, ticExptemp] = DenStepperAB2Pf( ...
      Lop, rhoVec_FT, GammaExVec_FT, GammaExVec_FTprev, NlPf, NlPrevPf, dt  );
    elseif( ParamObj.StepMeth == 2 ) 
       %[rhoVec_FTnext, ticExptemp] = DenStepperHAB1( ...
       %Lop, rhoVec_FT, GammaExVec_FT,dt );
      [rhoVec_FTnext, ticExptemp] = DenStepperHAB1Pf( ...
      Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt );
    elseif( ParamObj.StepMeth == 3 )
       %[rhoVec_FTnext, ticExptemp] = DenStepperHAB2( ...
          %Lop, rhoVec_FT, GammaExVec_FT,GammaExVec_FTprev,dt );
      [rhoVec_FTnext, ticExptemp] = DenStepperHAB2Pf( ...
         Lop, rhoVec_FT, GammaExVec_FT, GammaExVec_FT, NlPf, NlPrevPf, dt  );
    elseif( ParamObj.StepMeth == 4 )
       %[rhoVec_FTnext, ticExptemp] = DenStepperBHAB1( ...
       %Lop, rhoVec_FT, GammaExVec_FT,dt );
      [rhoVec_FTnext, ticExptemp] = DenStepperBHAB1Pf( ...
      Lop, rhoVec_FT, GammaExVec_FT, NlPf, dt  );
    elseif( ParamObj.StepMeth == 5 )
      %[rhoVec_FTnext, ticExptemp] = DenStepperBHAB2( ...
      %Lop, rhoVec_FT, GammaExVec_FT,GammaExVec_FTprev,dt );
      [rhoVec_FTnext, ticExptemp] = DenStepperBHAB2Pf( ...
      Lop, rhoVec_FT, GammaExVec_FT,GammaExVec_FTprev, NlPf, NlExpPf, NlPrevPf, dt );
    elseif( ParamObj.StepMeth == 6 )
      [rhoVec_FTnext, ticExptemp] = DenStepperPhiV( ...
         Lop, rhoVec_FT, GammaExVec_FT, dt );
    end

        %Make sure things are taking too long. This is a sign density---> inf
    [ShitIsFucked] = ExpTooLongChecker(...
        wfid,ticExptemp,ticExpInt,rhoVec_FT,Nx,Ny,Nm,j_record);
    if ShitIsFucked
        j_record = j_record - 1;
        break
    end
    
    %Save everything (this includes the initial state)
    if (mod(t,TimeObj.N_count)== 0)
        if ParamObj.SaveMe
            [SteadyState,ShitIsFucked,MaxReldRho] = ...
                VarRecorderTracker(wfid,tfid,TimeObj,t,...
                Nx,Ny,Nm,rhoVec_FT,rhoVec_FTprev,TotalDensity ,j_record);
            
        else
            [SteadyState,ShitIsFucked,MaxReldRho] = ...
                VarRecorderTrackerNoSave(wfid,tfid,TimeObj,t,Nx,Ny,Nm,...
                rhoVec_FT,rhoVec_FTprev,TotalDensity);
        end
        if ShitIsFucked == 1 || SteadyState == 1
            break
        end
        j_record = j_record+1;
        %         keyboard
    end %end recording
    
end %end time loop
fprintf(lfid,'Finished master time loop\n');
%  keyboard

% Update last rho
if ParamObj.SaveMe 
    if ShitIsFucked == 0 && SteadyState == 0
        t =  t + 1;
        rhoVec_FTprev  = rhoVec_FT;
        rhoVec_FT      = rhoVec_FTnext;
        rho_FT         = reshape(rhoVec_FT,Nx,Ny,Nm);
        if ParamObj.SaveMe
            [SteadyState,ShitIsFucked,MaxReldRho] = ...
                VarRecorderTracker(wfid,tfid,TimeObj,t,...
                Nx,Ny,Nm,rhoVec_FT,rhoVec_FTprev,TotalDensity ,j_record);
            
        else
            [SteadyState,ShitIsFucked,MaxReldRho] = ...
                VarRecorderTrackerNoSave(wfid,tfid,TimeObj,t,Nx,Ny,Nm,...
                rhoVec_FT,rhoVec_FTprev,TotalDensity);
        end
    end
end %end if save

%If something broke, return zeros. Else, return the goods
if ShitIsFucked 
    fprintf(wfid,'Density is either negative or not conserved.\n');
    fprintf(wfid,'I have done %i steps out of %i.\n',t, TimeObj.N_time);
    
elseif SteadyState 
    fprintf(wfid,'Things are going steady if you know what I mean.\n');
    fprintf(wfid,'I have done %i steps out of %i.\n',t, TimeObj.N_time);
end

% Get rid of zeros in record matrices
Record_hold   = 1:j_record;
TimeRecVec    = (0:j_record-1) * TimeObj.t_record;
if ParamObj.SaveMe 
    Density_rec   = Density_rec(:,:,:,Record_hold);
    DensityFT_rec = DensityFT_rec(:,:,:,Record_hold);
else
    Density_rec   = rho;
    DensityFT_rec = rho_FT;
end %end if save

trun = toc;

% See how much memory this used
% [uVbody, sVbody] = memory;

%Save the structure
% keyboard

DenRecObj = struct('DidIBreak', ShitIsFucked,'SteadyState', SteadyState,...
    'TimeRecVec',TimeRecVec,...
    'RunTime', trun, ...
    'bc',ParamObj.bc,...
    'Density_rec',Density_rec,'DensityFT_rec', DensityFT_rec,...
    'MaxReldRho',MaxReldRho,'feq',feq);


fclose(wfid); %Close warning statement file
fclose(tfid); %Close program tracker file
% keyboard
end %functi
