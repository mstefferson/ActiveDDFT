function [rho_FT_next] = DenStepperMWSGnrlCube(Prop,rho_FT, GammaEx_FT,GammaEx_FT_prev,dt)

ticExpInt = toc(ticExpIntID);

rho_FT_next = Prop .* ( rho_FT + dt/2 * ( 2 * GammaEx_FT - GammaEx_FT_prev ) ) ...
              + dt / 4 * ( 3 * GammaEx_FT - GammaEx_FT_prev ) ;


%
% ticExpIntID = tic;
% [rhoVec_FT_next, err] = ...
%             expv( dt, Lop, rhoVec_FT + dt / 2 * (3 * GammaEx_FT - GammaEx_FT_prev) );
% ticExpInt = toc(ticExpIntID);

end