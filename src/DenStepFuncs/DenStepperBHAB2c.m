% Better Hybrid AB2.
% Integrate using trapezoid method, say 
% Gamma_(n+1) = 2 * Gamma_n - Gamma_(n-1)

function [rho_FT_next] = ...
    DenStepperBHAB2c( Prop, rho_FT, GammaEx_FT, GammaEx_FTprev, dt)


rho_FT_next = Prop .* (rho_FT) + dt / 2 * ...
    (  (2 + Prop) .* GammaEx_FT - GammaEx_FTprev ) ;

end