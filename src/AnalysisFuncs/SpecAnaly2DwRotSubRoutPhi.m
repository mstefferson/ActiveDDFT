function [ampl_record] = ...
    SpecAnaly2DwRotSubRoutPhi(GridObj,ParamObj,NumModesMax,DensityFT_record,TimeRecVec, dt,kxholder,kyholder,...
    min_amp,kx,ky,km,D_rot, D_pos,Nx,Ny,Nm,bc)

DecayDisp   = 0;
AllKsVsTime = 1;

%Use squeeze to make a matrix (length(ky), length(j_record) ) of the
%amplitudes we want to look at
ampl_record = squeeze(DensityFT_record(kxholder,kyholder,:,:));

% keyboard
[k2plotInd,Nmodes] = Ks2TrackFinderSpec(NumModesMax,ampl_record,TimeRecVec,min_amp);

% keyboard
%Now plot all these amplitudes throughout the record
if DecayDisp
    if bc < 1.5 %Isotropic
        IsoDispPlotterBody(k2plotInd,ampl_record,Nmodes,TimeRecVec, dt,kxholder,kyholder,...
            kx,ky,km,D_rot, D_pos,Nx,Ny,Nm,bc)
    else % Nematic
        NemDispPlotterBody(GridObj,ParamObj,k2plotInd,ampl_record,Nmodes,TimeRecVec, kxholder,kyholder,...
    kx,ky,km,D_rot, D_pos,Nx,Ny,Nm,bc)

    end
    
end % End DecayDisp

if AllKsVsTime
    ModesVsTimePlotterKm(ampl_record,k2plotInd,TimeRecVec, Nmodes,kxholder,kyholder,...
        Nx,Ny,Nm)
end % End AllKsVsTime

end %end function