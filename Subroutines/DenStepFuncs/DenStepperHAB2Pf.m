% Hybrid AB 2
% Exponentiate actual propagator
% Uses a prefactor

function [rhoVec_FT_next,ticExpInt] = DenStepperHAB2Pf(...
    Lop ,rhoVec_FT, GammaEx_FT, GammaEx_FT_prev, NlPf, NlPfprev,dt)

%Take the first step in k-space using Euler. Save the time it takes
ticExpIntID = tic;

%Only Propagate if it is not zero
if GammaEx_FT_prev == 0
    GammaPrevPrpgtd = zeros(length(GammaEx_FT_prev),1);
else
[GammaPrevPrpgtd, err] = expv( dt, Lop, GammaEx_FT_prev) ;    
end

% [rhoVec_FT_next, err] = ...
%             expv( dt, Lop, rhoVec_FT + dt / 2 * (3 * GammaEx_FT - GammaEx_FT_prev) );
[rhoVec_FT_next, err] = ...
            expv( dt, Lop, rhoVec_FT + NlPf * GammaEx_FT - NlPfprev * GammaPrevPrpgtd );
ticExpInt = toc(ticExpIntID);



end
