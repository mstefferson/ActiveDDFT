%  Better Hybrid AB 2
% Exponentiate actual propagator

function [rhoVec_FT_next,ticExpInt] = DenStepperBHAB2(...
    Lop, rhoVec_FT, GammaEx_FT,GammaEx_FTprev,dt)

%Take the first step in k-space using Euler. Save the time it takes
ticExpIntID = tic;
[rhoVec_FT_next, err] = ...
                  expv( dt, Lop, rhoVec_FT + dt/2 * GammaEx_FT );
rhoVec_FT_next = rhoVec_FT_next + ...
    dt * ( GammaEx_FT - 1/2 * GammaEx_FTprev );
ticExpInt = toc(ticExpIntID);

end