% Uses the initial density indicator to choose the correct intial density
% subroutine.

function [rho] = MakeConc(systemObj,rhoInit,x,y,phi)

if rhoInit.IntCond == 0
  [rho] = IntDenCalcIsoPw2Drot(systemObj,rhoInit);
elseif rhoInit.IntCond == 1
  % if bc is too close to 1.5, errors arise. Fix this here.
  if 1.499 < systemObj.bc && systemObj.bc < 1.501
    systemObjTemp = systemObj;
    systemObjTemp.bc = 1.502;
    % Initial distribution
    [rho] = IntDenCalcEqPw2Drot(systemObjTemp,rhoInit,x,y,phi);
    % Perturb it
    [rho] = PwPerturbFT(rho,systemObjTemp,rhoInit);
  else
    % Initial distribution
    [rho] = IntDenCalcEqPw2Drot(systemObj,rhoInit,x,y,phi);
    % Perturb it
    [rho] = PwPerturbFT(rho,systemObj,rhoInit);
  end
elseif rhoInit.IntCond == 2
  % Initial distribution
  [rho] = IntDenCalcNemPw2rot(systemObj,rhoInit,x,y,phi);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif rhoInit.IntCond == 3
  [rho] = IntDenCalcLoaded2Drot(rhoInit.LoadName);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif rhoInit.IntCond == 4
  %[rho] = IntDenCalcIsoSepPw2Drot(gridObj,ParamObj,rhoInit);
  fprintf('Error not written');
  error('Error not written');
  rho = 0;
elseif rhoInit.IntCond == 5
  %[rho] = IntDenCalcEqSepPw2Drot(gridObj,ParamObjTemp,rhoInit);
  fprintf('Error not written');
  error('Error not written');
  rho = 0;
elseif rhoInit.IntCond == 6
  % Initial distribution
  [rho] = IntDenCalcGauss2Drot(gridObj,ParamObj,rhoInit);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
end

end
