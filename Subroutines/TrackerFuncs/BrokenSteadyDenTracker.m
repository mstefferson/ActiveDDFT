% [SteadyState,ShitIsFucked,MaxReldRho] = ...
%    BrokenSteadyDenTracker(rho,rho_prev,TotalDensity ,TimeObj)
%
% Description: Tracks broken densities and steady state

function [SteadyState,ShitIsFucked,MaxReldRho] = ...
    BrokenSteadyDenTracker(rho,rho_prev,TotalDensity ,TimeObj)
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
    fprintf('Forgive me, your grace. Density has become negative\n');
    ShitIsFucked  = 1;
end

%Not conserving density check.
if abs( sum(sum(sum(rho)))- TotalDensity ) > TotalDensity / 1000;
    fprintf('Forgive me, your grace. Density is not being conserved\n');
    ShitIsFucked  = 1;
end


% Nan or infinity
if find(isinf(rho)) ~= 0
    fprintf('Forgive me, your grace. Density has gone infinite. ');
    fprintf('Does that make sense? No. No it does not\n');
    ShitIsFucked  = 1;
end

if find(isnan(rho)) ~= 0
    fprintf('Forgive me, your grace. Density elements are no longer numbers. ');
    fprintf('Does that make sense? No. No it does not\n');
    ShitIsFucked  = 1;
end


end
