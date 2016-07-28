function [rhoVec_FT_next,ticExpInt] = DenStepperMWSGnrl(Lop,rhoVec_FT, GammaEx_FT,GammaEx_FT_prev,dt)

ticExpIntID = tic;
% keyboard

[rhoVec_FT_next, err] = ...
    expv( dt, Lop, (rhoVec_FT + dt / 2 * (2 * GammaEx_FT - GammaEx_FT_prev) ) );

rhoVec_FT_next = rhoVec_FT_next + dt / 4 * ( 3 * GammaEx_FT - GammaEx_FT_prev );
% keyboard

ticExpInt = toc(ticExpIntID);


%
% ticExpIntID = tic;
% [rhoVec_FT_next, err] = ...
%             expv( dt, Lop, rhoVec_FT + dt / 2 * (3 * GammaEx_FT - GammaEx_FT_prev) );
% ticExpInt = toc(ticExpIntID);

end