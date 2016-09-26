% [SteadyState,ShitIsFucked,MaxReldRho] = ...
%    BrokenSteadyDenTracker(rho,rho_prev,TotalDensity ,timeObj)
%
% Description: Tracks broken densities and steady state

function [SteadyState,ShitIsFucked,MaxReldRho] = ...
  BrokenSteadyDenTracker(rho, rhoPrev, rho_FT, constConc,timeObj,systemObj)

SteadyState = 0;
ShitIsFucked = 0;

WeightDensityChange =  abs( 1 - rhoPrev ./ rho  );

MaxReldRho = max( max( max(WeightDensityChange) ) );

if MaxReldRho < timeObj.ss_epsilon
  SteadyState = 1;
end

%See if something broke

%Negative Density check
if rho < 0
  fprintf('Forgive me, your grace. Density has become negative\n');
  ShitIsFucked  = 1;
end

%Not conserving density check.
if systemObj.Nm ~= 1
  constConcNow = rho_FT(systemObj.Nx/2+1, systemObj.Ny/2+1,systemObj.Nm/2+1 );
  if ( abs( constConcNow - constConc ) / constConc ) >  1e-10
    fprintf('Forgive me, your grace. Density is not being conserved\n');
    ShitIsFucked  = 1;
  end
end

% Nan or infinity
if isinf(rho)
  fprintf('Forgive me, your grace. Density has gone infinite. ');
  fprintf('Does that make sense? No. No it does not\n');
  ShitIsFucked  = 1;
end

if isnan(rho)
  fprintf('Forgive me, your grace. Density elements are no longer numbers. ');
  fprintf('Does that make sense? No. No it does not\n');
  ShitIsFucked  = 1;
end


end
