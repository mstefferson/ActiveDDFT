function [rhoVec_FT_next,ticExpInt] = DenStepperMDBGnrl(Lop,rhoVec_FT, GammaEx_FT,GammaEx_FT_prev,dt)

ticExpIntID = tic;

%Only Propagate if it is not zero
if GammaEx_FT_prev == 0
    GammaPrevPrpgtd = zeros(length(GammaEx_FT_prev),1);
else
[GammaPrevPrpgtd, err] = expv( dt, Lop, GammaEx_FT_prev) ;    
end
% keyboard
% [rhoVec_FT_next, err] = ...
%             expv( dt, Lop, rhoVec_FT + dt / 2 * (3 * GammaEx_FT - GammaEx_FT_prev) );
[rhoVec_FT_next, err] = ...
    expv( dt, Lop, (rhoVec_FT + dt / 2 * (3 * GammaEx_FT - GammaPrevPrpgtd) ) );

% keyboard

ticExpInt = toc(ticExpIntID);


%
% ticExpIntID = tic;
% [rhoVec_FT_next, err] = ...
%             expv( dt, Lop, rhoVec_FT + dt / 2 * (3 * GammaEx_FT - GammaEx_FT_prev) );
% ticExpInt = toc(ticExpIntID);

end