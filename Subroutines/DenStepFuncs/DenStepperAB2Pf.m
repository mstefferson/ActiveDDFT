% AB 2
% Exponentiate actual propagator
% With pre factor
function [rhoVec_FT_next,ticExpInt] = DenStepperAB2Pf(...
    Lop,rhoVec_FT, GammaEx_FT,GammaEx_FT_prev,NlPf, NlPfprev, dt)

%Take the first step in k-space using Euler. Save the time it takes
ticExpIntID = tic;

[rhoVec_FT_next, err] = ...
            expv( dt, Lop, rhoVec_FT) ;
ticExpInt = toc(ticExpIntID);

rhoVec_FT_next = rhoVec_FT_next  + NlPf * GammaEx_FT - NlPfprev * GammaEx_FT_prev ;

end
