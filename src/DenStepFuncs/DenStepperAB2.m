% AB 2
% Exponentiate actual propagator

function [rhoVec_FT_next,ticExpInt] = DenStepperAB2(Lop,rhoVec_FT, GammaEx_FT,GammaEx_FT_prev,dt)

%Take the first step in k-space using Euler. Save the time it takes
ticExpIntID = tic;

[rhoVec_FT_next, err] = ...
            expv( dt, Lop, rhoVec_FT) ;
ticExpInt = toc(ticExpIntID);

rhoVec_FT_next = rhoVec_FT_next  + dt / 2 * ( 3 * GammaEx_FT - GammaEx_FT_prev );

end