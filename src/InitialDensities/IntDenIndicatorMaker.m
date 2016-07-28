function [IntDenType, IntDenIndicator] = ...
    IntDenIndicatorMaker(IntGauss,IntPw,IntSepPw,IntEqPw,IntEqSepPw,IntLoad,IntNemPw)

if (IntGauss + IntPw + IntSepPw + IntEqPw +  IntEqSepPw + IntLoad + IntNemPw) ~= 1
    error('Choose one and only one intial conidtion your holiness')
end

if IntGauss
    IntDenType      = sprintf('Gaussian');
    IntDenIndicator = sprintf('G');
elseif IntPw
    IntDenType      = sprintf('PlaneWave');
    IntDenIndicator = sprintf('P');
elseif IntSepPw
    IntDenType      = sprintf('SepPlaneWave');
    IntDenIndicator = sprintf('SP');
elseif IntEqPw
    IntDenType      = sprintf('EquilibriumPW');
    IntDenIndicator = sprintf('EP');
elseif IntEqSepPw
    IntDenType      = sprintf('EquilibriumSPW');
    IntDenIndicator = sprintf('ESP');
elseif IntLoad
    IntDenType      = sprintf('Loaded');
    IntDenIndicator = sprintf('L');
elseif IntNemPw
    IntDenType      = sprintf('NematicPW');
    IntDenIndicator = sprintf('NP');
    
end

end