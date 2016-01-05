% Currently causes things to blow up
% Exponential Euler Methodf

function [rho_FT_next] = ...
    DenStepperEEMc1( Prop, GamProp, rho_FT, GammaEx_FT)

%Take the first step in k-space using Euler. Everything is a cube.
%All this is being done element by element.

rho_FT_next = Prop .* rho_FT + GamProp .* GammaEx_FT;

end