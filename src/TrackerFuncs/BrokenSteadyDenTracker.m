% [steadyState,shitIsFucked,maxDrho] = ...
%    BrokenSteadyDenTracker(rho,rho_prev,TotalDensity ,timeObj)
%
% Description: Tracks broken densities and steady state

function [steadyState,shitIsFucked, whatBroke, maxDrho] = ...
  BrokenSteadyDenTracker(rho, rhoPrev, rho_FT, constConc,timeObj,systemObj)
% Initialize
steadyState = 0;
shitIsFucked = 0;
whatBroke = [];
% Max density change
WeightDensityChange =  abs( 1 - rhoPrev ./ rho  );
maxDrho = max( max( max(WeightDensityChange) ) );
if maxDrho < timeObj.ss_epsilon_dt
  steadyState = 1;
  keyboard
end
%See if something broke
%Negative Density check
if any( rho(:) < 0 )
  fprintf('Forgive me, your grace. Density has become negative\n');
  shitIsFucked  = 1;
  whatBroke = 'negative density';
end
%Not conserving density check.
if systemObj.n3 ~= 1
  constConcNow = rho_FT(systemObj.n1/2+1, systemObj.n2/2+1,systemObj.n3/2+1 );
  if ( abs( constConcNow - constConc ) / constConc ) >  1e-10
    fprintf('Forgive me, your grace. Density is not being conserved\n');
    shitIsFucked  = 1;
    whatBroke = 'density not conserved';
  end
end
% Nan or infinity
if any( isinf( rho(:) ) )
  fprintf('Forgive me, your grace. Density has gone infinite. ');
  shitIsFucked  = 1;
  whatBroke = 'density inf';
end

if any( isnan( rho(:) ) )
  fprintf('Forgive me, your grace. Density elements are no longer numbers. ');
  shitIsFucked  = 1;
  whatBroke = 'density nan';
end
