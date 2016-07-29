% Hybrid AB 2
% Exponentiate actual propagator


function [rhoVec_FT_next,ticExpInt] = DenStepperHAB2(Lop,rhoVec_FT, GammaEx_FT,GammaEx_FT_prev,dt)

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
            expv( dt, Lop, rhoVec_FT + dt / 2 * (3 * GammaEx_FT - GammaPrevPrpgtd) );
ticExpInt = toc(ticExpIntID);



end