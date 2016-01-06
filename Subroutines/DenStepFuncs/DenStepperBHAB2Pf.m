%  Better Hybrid AB 2
% Exponentiate actual propagator
% Prefactor

function [rhoVec_FT_next,ticExpInt] = DenStepperBHAB2Pf(...
    Lop, rhoVec_FT, GammaEx_FT,GammaEx_FTprev, NlPf, NlExpPf, NlPfprev,dt)

%Take the first step in k-space using Euler. Save the time it takes
ticExpIntID = tic;
[rhoVec_FT_next, err] = ...
                  expv( dt,Lop, rhoVec_FT + NlExpPf * GammaEx_FT );
rhoVec_FT_next = rhoVec_FT_next + ...
    NlPf *  GammaEx_FT - NlPfprev * GammaEx_FTprev ;
ticExpInt = toc(ticExpIntID);

end
