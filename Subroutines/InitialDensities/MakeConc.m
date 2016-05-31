% Uses the initial density indicator to choose the correct intial density
% subroutine.

function [rho] = MakeConc(GridObj,ParamObj,RhoInit)

if RhoInit.IntCond == 0
  [rho] = IntDenCalcIsoPw2Drot(GridObj,ParamObj,RhoInit);
elseif RhoInit.IntCond == 1
  % if bc is too close to 1.5, errors arise. Fix this here.
  if 1.499 < ParamObj.bc && ParamObj.bc < 1.501
    ParamObjTemp = ParamObj;
    ParamObjTemp.bc = 1.502;
    [rho] = IntDenCalcEqPw2Drot(GridObj,ParamObjTemp,RhoInit);
  else
    [rho] = IntDenCalcEqPw2Drot(GridObj,ParamObj,RhoInit);
  end
elseif RhoInit.IntCond == 2
  [rho] = IntDenCalcNemPw2rot(GridObj,ParamObj,RhoInit);
elseif RhoInit.IntCond == 3
  [rho] = IntDenCalcLoaded2Drot(RhoInit.LoadName);
elseif RhoInit.IntCond == 4
  [rho] = IntDenCalcIsoSepPw2Drot(GridObj,ParamObj,RhoInit);
elseif RhoInit.IntCond == 5
  % if bc is too close to 1.5, errors arise. Fix this here.
  if 1.499 < ParamObj.bc && ParamObj.bc < 1.501
    ParamObjTemp = ParamObj;
    ParamObjTemp.bc = 1.502;
    [rho] = IntDenCalcEqSepPw2Drot(GridObj,ParamObjTemp,RhoInit);
  else
    [rho] = IntDenCalcEqSepPw2Drot(GridObj,ParamObj,RhoInit);
  end
elseif RhoInit.IntCond == 6
  [rho] = IntDenCalcGauss2Drot(GridObj,ParamObj,RhoInit);
end

end
