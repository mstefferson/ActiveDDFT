function [SteadyState,ShitIsFucked,MaxReldRho] = ...
VarRecorderTracker(wfid,tfid,TimeObj,t,Nx,Ny,Nm,rhoVec_FT,rhoVec_FT_prev,TotalDensity ,j_record)
% Track how mucht the wieghted density has changed.
%Check to see if steady state has been reached. If so, break the
%loop'
%         keyboard
global Density_rec 
global DensityFT_rec

% keyboard
fprintf(tfid,'%f percent done\n',t./TimeObj.N_time*100);
% fclose(tfid);
rho         = real(ifftn(ifftshift(reshape( rhoVec_FT,Nx,Ny,Nm ))));
rho_prev    = real(ifftn(ifftshift(reshape( rhoVec_FT_prev,Nx,Ny,Nm ))));

% See if things are broken
[SteadyState,ShitIsFucked,MaxReldRho] = BrokenSteadyDenTracker(wfid,rho,rho_prev,TotalDensity ,TimeObj);

rho_cube_FT = reshape( rhoVec_FT,Nx,Ny,Nm );
DensityFT_rec(:,:,:,j_record)   = rho_cube_FT;
Density_rec(:,:,:,j_record)     = rho;
%         keyboard
end

function [SteadyState,ShitIsFucked,MaxReldRho] = ...
    BrokenSteadyDenTracker(wfid,rho,rho_prev,TotalDensity ,TimeObj)
SteadyState = 0;
ShitIsFucked = 0;

AbsDensityChange = abs( rho - rho_prev );
WeightDensityChange = AbsDensityChange ./ rho;
MaxReldRho = max(max(max(WeightDensityChange)));
if MaxReldRho < TimeObj.ss_epsilon
    SteadyState = 1;
end

%See if something broke
%Negative Density check
if min(min(min(rho))) < 0
    fprintf(wfid,'Forgive me, your grace. Density has become negative\n');
    %         keyboard
    ShitIsFucked  = 1;
end
%Not conserving density check.
if abs( sum(sum(sum(rho)))- TotalDensity ) > TotalDensity / 1000;
    fprintf(wfid,'Forgive me, your grace. Density is not being conserved\n');
    ShitIsFucked  = 1;
end


% Nan or infinity
% keyboard
if find(isinf(rho)) ~= 0
    fprintf(wfid,'Forgive me, your grace. Density has gone infinite. ');
    fprintf(wfid,'Does that make sense? No. No it does not\n');
    ShitIsFucked  = 1;
end

if find(isnan(rho)) ~= 0
    fprintf(wfid,'Forgive me, your grace. Density elements are no longer numbers. ');
    fprintf(wfid,'Does that make sense? No. No it does not\n');
    ShitIsFucked  = 1;
end


end