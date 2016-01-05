
% Adams Bashforth 2
% Ignores integrating factor.

function [rho_FT_next] = ...
    DenStepperAB2cPf( Prop, rho_FT, GammaEx_FT,GammaEx_FTprev, NlPf, NlPfprev)

%Take the first step in k-space using Euler. Everything is a cube.
%All this is being done element by element.

rho_FT_next = Prop .* rho_FT + NlPf * GammaEx_FT - NlPfprev .* GammaEx_FTprev ;

end
