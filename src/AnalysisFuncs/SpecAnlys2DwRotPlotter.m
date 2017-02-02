
function [] = ...
    SpecAnlys2DwRotPlotter(gridObj,ParamObj,DensityFT_record, D_rot, D_pos, n1,n2,n3, bc, ...
                           kx,ky,km,TimeRecVec,dt,xDispersion,yDispersion,angDispersion )
%Subroutine analyzes the decaying modes of the system and plots the
%dispersion relation for the system.

% Used for 2-d system with an orientational degree of freedom. Diffusion
% tensor is fram independent (Not rods)


%Plot the decay of a certain k-space amplitude
% Let's plot the different ky amplitudes with kx = 0 for different times.
%Find these amplitudes in the FT record


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% k_x Dispersion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if xDispersion == 1
    min_amp = 0.001;         % Minimum the amplitude value  needs to be to track
    NumModesMax = 5;       % Find which ks you want to track.
    kyholder = n2/2 + 2;    % Place holder of kx  = 0;
    kmholder = n3/2 + 5;    % Place holder of ky = 0;
    [ampl_record_x,fitParam_x] = SpecAnaly2DwRotSubRoutX(NumModesMax,DensityFT_record,...
                                 TimeRecVec, dt,kyholder,kmholder,...
                                 min_amp,kx,ky,km,D_rot, D_pos,n1,n2,n3,bc);

end %end if kx dispersion


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% k_y Dispersion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if yDispersion == 1
    min_amp = 0.001;         % Minimum the amplitude value  needs to be to track
    NumModesMax = 5;       % Find which ks you want to track.
    kxholder = n1/2 + 1;    % Place holder of kx  = 0;
    kmholder = n3/2 + 1;    % Place holder of ky = 0;
    [ampl_record_y,fitParam_y] = SpecAnaly2DwRotSubRoutY(NumModesMax,DensityFT_record, ...
                                 TimeRecVec, dt,kxholder,kmholder,...
                                 min_amp,kx,ky,km,D_rot, D_pos,n1,n2,n3,bc);

end %end if ky dispersion


%%%%%%%%%%%%%%%%%%%%%%%%% k_m dispersion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if angDispersion == 1
    min_amp = 0.0001;         % Minimum the amplitude value  needs to be to track
    NumModesMax = 15;       % Find which ks you want to track. Must be odd
    kxholder = n1/2 + 1 ;    % Place holder of kx  = 0;
    kyholder = n2/2 + 1;    % Place holder of ky = 0;
    [ampl_record_m] = SpecAnaly2DwRotSubRoutPhi(gridObj,ParamObj,NumModesMax,DensityFT_record,...
                                 TimeRecVec,dt ,kxholder,kyholder, min_amp,kx,ky,km,...
                                 D_rot, D_pos,n1,n2,n3,bc);
%     keyboard
end %end if ang dispersion
% keyboard
end %end function