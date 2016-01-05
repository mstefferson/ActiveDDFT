% Uses the initial density indicator to choose the correct intial density
% subroutine.

function [rho] = MakeConc(GridObj,ParamObj,LoadName)

  if ParamObj.IntCond == 0
        [rho] = IntDenCalcPwModes2Drot(GridObj,ParamObj);
  elseif ParamObj.IntCond == 1
        [rho] = IntDenCalcPerturbEqPw2Drot(GridObj,ParamObj);
  elseif ParamObj.IntCond == 2
        [rho] = IntDenCalcPerturbNemPw2Drot(GridObj,ParamObj);
  elseif ParamObj.IntCond == 3
        [rho] = IntDenCalcLoaded2Drot(LoadName);
  elseif ParamObj.IntCond == 4
        [rho] = IntDenCalcSepPwModes2Drot(GridObj,ParamObj);
  elseif ParamObj.IntCond == 5
        [rho] = IntDenCalcPerturbEqSepPw2Drot(GridObj,ParamObj);
  elseif ParamObj.IntCond == 6
        [rho] = IntDenCalcGauss2Drot(GridObj,ParamObj);
  end
    
end
