% Uses the initial density indicator to choose the correct intial density
% subroutine.

function [rho] = MakeConcFromInd(gridObj,ParamObj,IntDenType)

 if strcmp(IntDenType,'PlaneWave')
        [rho] = IntDenCalcPwModes2Drot(gridObj,ParamObj);
    elseif strcmp(IntDenType,'SepPlaneWave')
        [rho] = IntDenCalcSepPwModes2Drot(gridObj,ParamObj);
    elseif strcmp(IntDenType,'Gaussian')
        [rho] = IntDenCalcGauss2Drot(gridObj,ParamObj);
    elseif strcmp(IntDenType,'EquilibriumPW')
        [rho] = IntDenCalcPerturbEqPw2Drot(gridObj,ParamObj);
    elseif strcmp(IntDenType,'EquilibriumSPW')
        [rho] = IntDenCalcPerturbEqSepPw2Drot(gridObj,ParamObj);
    elseif strcmp(IntDenType,'Loaded')
        [rho] = IntDenCalcLoaded2Drot(DataTemp.textdata{4});
    elseif strcmp(IntDenType,'NematicPW')
        [rho] = IntDenCalcPerturbNemPw2Drot(gridObj,ParamObj);
 end
    
end