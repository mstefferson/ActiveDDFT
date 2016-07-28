% Hybrid AB2 with prefactor

function [rho_FT_next] = DenStepperHAB2cPf(...
    Prop,rho_FT, GammaEx_FT,GammaEx_FT_prev,NlPf,NlPfprev)

% Step using the hybrid AB. Everything is still in cube form.
rho_FT_next = Prop .*  rho_FT + NlPf .* GammaEx_FT - NlPfprev .* GammaEx_FT_prev ;

end