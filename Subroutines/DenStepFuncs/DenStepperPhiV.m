% Use expokits phiv.

function [rhoVec_FT_next,ticExpInt] = DenStepperPhiV( Lop, rhoVec_FT, GammaEx_FT,dt)

% Use phiV
ticExpIntID = tic;
[rhoVec_FT_next, err] = phiv( dt, Lop, GammaEx_FT, rhoVec_FT );
ticExpInt = toc(ticExpIntID);


end
