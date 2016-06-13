% Uses the initial density indicator to choose the correct intial density
% subroutine.

function [rho] = MakeConc(ParamObj,RhoInit,x,y,phi)

if RhoInit.IntCond == 0
  [rho] = IntDenCalcIsoPw2Drot(ParamObj,RhoInit);
elseif RhoInit.IntCond == 1
  % if bc is too close to 1.5, errors arise. Fix this here.
  if 1.499 < ParamObj.bc && ParamObj.bc < 1.501
    ParamObjTemp = ParamObj;
    ParamObjTemp.bc = 1.502;
    [rho] = IntDenCalcEqPw2Drot(ParamObjTemp,RhoInit,x,y,phi);
  else
    [rho] = IntDenCalcEqPw2Drot(ParamObj,RhoInit,x,y,phi);
  end
elseif RhoInit.IntCond == 2
  [rho] = IntDenCalcNemPw2rot(ParamObj,RhoInit,x,y,phi);
elseif RhoInit.IntCond == 3
  [rho] = IntDenCalcLoaded2Drot(RhoInit.LoadName);
elseif RhoInit.IntCond == 4
  %[rho] = IntDenCalcIsoSepPw2Drot(GridObj,ParamObj,RhoInit);
  fprintf('Error not written');
  error('Error not written');
  rho = 0;
elseif RhoInit.IntCond == 5
  % if bc is too close to 1.5, errors arise. Fix this here.
  if 1.499 < ParamObj.bc && ParamObj.bc < 1.501
    ParamObjTemp = ParamObj;
    ParamObjTemp.bc = 1.502;
  else
  %[rho] = IntDenCalcEqSepPw2Drot(GridObj,ParamObjTemp,RhoInit);
  fprintf('Error not written');
  error('Error not written');
  rho = 0;
  end
elseif RhoInit.IntCond == 6
  [rho] = IntDenCalcGauss2Drot(GridObj,ParamObj,RhoInit);
end

end
