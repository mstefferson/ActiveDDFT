% AB 1
% Exponentiate actual propagator

function [rhoVec_FT_next,ticExpInt] = DenStepperAB1( Lop, rhoVec_FT, GammaEx_FT,dt)

%Take the first step in k-space using Euler. Save the time it takes
ticExpIntID = tic;
[rhoVec_FT_next, err] = expv( dt, Lop, rhoVec_FT );
             
rhoVec_FT_next = rhoVec_FT_next + dt * GammaEx_FT;

ticExpInt = toc(ticExpIntID);


end