% Uses the initial density indicator to choose the correct intial density
% subroutine.

function [rho] = MakeConc(systemObj,rhoInit,gridObj)

if rhoInit.IntCond == 0
  [rho] = IntDenCalcIsoPw2Drot(systemObj,rhoInit);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif rhoInit.IntCond == 1
  % if bc is too close to 1.5, errors arise. Fix this here.
  if 1.499 < systemObj.bc && systemObj.bc < 1.501
    systemObjTemp = systemObj;
    systemObjTemp.bc = 1.502;
    % Initial distribution
    [rho] = IntDenCalcEqPw2Drot(...
      systemObjTemp,rhoInit,gridObj.x,gridObj.y,gridObj.phi);
    % Perturb it
    [rho] = PwPerturbFT(rho,systemObjTemp,rhoInit);
  else
    % Initial distribution
    [rho] = IntDenCalcEqPw2Drot(systemObj,rhoInit,gridObj.x,gridObj.y,gridObj.phi);
    % Perturb it
    [rho] = PwPerturbFT(rho,systemObj,rhoInit);
  end
elseif rhoInit.IntCond == 2
  % Initial distribution
  [rho] = IntDenCalcNemPw2rot(systemObj,rhoInit,gridObj.x,gridObj.y,gridObj.phi);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif rhoInit.IntCond == 3
  [rho] = IntDenCalcLoaded2Drot(rhoInit.LoadName);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif rhoInit.IntCond == 4 || rhoInit.IntCond == 5
  % Build eq first
  % if bc is too close to 1.5, errors arise. Fix this here.
  if 1.499 < systemObj.bc && systemObj.bc < 1.501
    systemObjTemp = systemObj;
    systemObjTemp.bc = 1.502;
    % Initial distribution
    [rho] = IntDenCalcEqPw2Drot(systemObjTemp,rhoInit,gridObj.x,gridObj.y,gridObj.phi);
  else
    % Initial distribution
    [rho] = IntDenCalcEqPw2Drot(systemObj,rhoInit,gridObj.x,gridObj.y,gridObj.phi);
  end
  % Gaussian perturb
  [rho] = polarPerturbGauss( rho, systemObj, rhoInit, gridObj );
else
  fprintf('Error not written');
  error('Error not written');
  rho = 0;
end

end
