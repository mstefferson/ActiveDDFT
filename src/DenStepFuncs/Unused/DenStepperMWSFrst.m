function [rhoVec_FT_next,ticExpInt] = DenStepperMWSFrst( Lop, rhoVec_FT, GammaEx_FT,dt)

%Take the first step in k-space using Euler. Save the time it takes
ticExpIntID = tic;
[rhoVec_FT_next, err] = ...
    expv( dt, Lop, rhoVec_FT + dt/2 * GammaEx_FT );
rhoVec_FT_next = rhoVec_FT_next + dt / 2 * GammaEx_FT;
ticExpInt = toc(ticExpIntID);


end