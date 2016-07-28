% Uses the initial density indicator to choose the correct intial density
% subroutine.

function [rho] = MakeConcFromInd(GridObj,ParamObj,IntDenType)

 if strcmp(IntDenType,'PlaneWave')
        [rho] = IntDenCalcPwModes2Drot(GridObj,ParamObj);
    elseif strcmp(IntDenType,'SepPlaneWave')
        [rho] = IntDenCalcSepPwModes2Drot(GridObj,ParamObj);
    elseif strcmp(IntDenType,'Gaussian')
        [rho] = IntDenCalcGauss2Drot(GridObj,ParamObj);
    elseif strcmp(IntDenType,'EquilibriumPW')
        [rho] = IntDenCalcPerturbEqPw2Drot(GridObj,ParamObj);
    elseif strcmp(IntDenType,'EquilibriumSPW')
        [rho] = IntDenCalcPerturbEqSepPw2Drot(GridObj,ParamObj);
    elseif strcmp(IntDenType,'Loaded')
        [rho] = IntDenCalcLoaded2Drot(DataTemp.textdata{4});
    elseif strcmp(IntDenType,'NematicPW')
        [rho] = IntDenCalcPerturbNemPw2Drot(GridObj,ParamObj);
 end
    
end