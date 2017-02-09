% Handles all the dRho contributions that are not in Lop
function [GammaCube_FT, shitIsFucked] = dRhoMaster( rho, rho_FT, ...
  flags, interObj,  systemObj, diffObj, particleObj, cosPhi3, sinPhi3,t )
% Initialize
GammaCube_FT = 0;
shitIsFucked = 0;
% Interactions
% mayers

if interObj.hardId == 1 % mayers
  muExFt = muExCalcVc2Ft(rho_FT, interObj.FmFt,systemObj,interObj.muExScale);
  GammaExCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
  GammaCube_FT = GammaCube_FT + GammaExCube_FT;
end
% spt
if interObj.hardId == 2 % spt
  if interObj.typeId == 2 % disks
    if systemObj.n3 == 1
      nu = interObj.b .* rho;
    else
      nu = interObj.b .* ifftn( ifftshift( rho_FT(:,:,interObj.k30) ) );
    end
    if nu <  1
      [muExFt] = muExDisksSPT(nu);
    else
      error('Density is too high!')
      fprintf('Density is too high!');
      shitIsFucked = 1;
      muExFt = 0;
    end
    GammaExCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
    GammaCube_FT = GammaCube_FT + GammaExCube_FT;
  end
end
% keyboard
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
