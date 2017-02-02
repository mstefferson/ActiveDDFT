% Uses the initial density indicator to choose the correct intial density
% subroutine.

function [rho] = MakeConc(systemObj,rhoInit,gridObj)

if rhoInit.IntCond == 0
  [rho] = IntDenCalcIsoPw2Drot(systemObj);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif rhoInit.IntCond == 1
  % if bc is too close to 1.5, errors arise. Fix this here.
  if 1.499 < systemObj.bc && systemObj.bc < 1.501
    systemObjTemp = systemObj;
    systemObjTemp.bc = 1.502;
    % Initial distribution
    [rho] = IntDenCalcEqPw2Drot(systemObjTemp,rhoInit);
    % Perturb it
    [rho] = PwPerturbFT(rho,systemObjTemp,rhoInit);
  else
    % Initial distribution
    [rho] = IntDenCalcEqPw2Drot(systemObj,rhoInit);
    % Perturb it
    [rho] = PwPerturbFT(rho,systemObj,rhoInit);
  end
elseif rhoInit.IntCond == 2
  % Initial distribution
  [rho] = IntDenCalcNemPw2rot(systemObj,gridObj.phi,rhoInit.shiftAngle);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif rhoInit.IntCond == 3
  [rho] = IntDenCalcLoaded2Drot(rhoInit.LoadName,systemObj);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif rhoInit.IntCond == 4 || rhoInit.IntCond == 5
  % Build eq first
  % if bc is too close to 1.5, errors arise. Fix this here.
  if 1.499 < systemObj.bc && systemObj.bc < 1.501
    systemObjTemp = systemObj;
    systemObjTemp.bc = 1.502;
    % Initial distribution
    [rho] = IntDenCalcEqPw2Drot(systemObjTemp,rhoInit);
  else
    % Initial distribution
    [rho] = IntDenCalcEqPw2Drot(systemObj,rhoInit);
  end
  % Gaussian perturb
  [rho] = polarPerturbGauss( rho, systemObj, rhoInit, gridObj );
elseif rhoInit.IntCond == 6
  % delta function in polar order
  [rho] = deltaPolarIc(systemObj);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
else
  fprintf('Error not written');
  error('Error not written');
  rho = 0;
end
% Renormalize out here
% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
if systemObj.Nm > 1
  int1 = trapz_periodic(gridObj.phi,rho,3);
else
  int1 = rho;
end
CurrentNorm = trapz_periodic(gridObj.x,...
  trapz_periodic(gridObj.y,int1,2),1);
rho = rho .* systemObj.numPart ./ CurrentNorm;

end
