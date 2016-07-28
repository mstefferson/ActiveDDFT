% Hybrid AB 1
% Exponentiate actual propagator.
% Use a Pf

function [rhoVec_FT_next,ticExpInt] = DenStepperHAB1Pf( ...
    Lop, rhoVec_FT, GammaEx_FT, NlPf,dt)

%Take the first step in k-space using Euler. Save the time it takes
ticExpIntID = tic;
[rhoVec_FT_next, err] = ...
                  expv( dt, Lop, rhoVec_FT + NlPf * GammaEx_FT );
ticExpInt = toc(ticExpIntID);


end
