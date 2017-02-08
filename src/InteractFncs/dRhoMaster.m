% Handles all the dRho contributions that are not in Lop
function [GammaCube_FT] = dRhoMaster( rho, rho_FT, ...
  flags, interObj,  systemObj, diffObj, particleObj, cosPhi3, sinPhi3 )
% Initialize
GammaCube_FT = 0;
% Interactions
% Hard rod
if interObj.hard == 'mayer'
  GammaExCube_FT = dRhoIntCalcVcFt( rho, rho_FT, interObj.FmFt, ...
    systemObj, diffObj, interObj.muExScale );
  GammaCube_FT = GammaCube_FT + GammaExCube_FT;
end
% Driving
if flags.Drive && flags.DiagLop
  GammaDrCube_FT  = ...
    dRhoDriveCalcFtId(rho,particleObj.vD,...
    cosPhi3, sinPhi3,diffObj.ik1rep3,diffObj.ik2rep3);
  GammaCube_FT = GammaCube_FT + GammaDrCube_FT;
end
% Hard rod
if interObj.extFlag
  GammaPotCube_FT = 0;
  GammaCube_FT = GammaCube_FT + GammaPotCube_FT;
end
