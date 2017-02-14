% Handles all the dRho contributions that are not in Lop
function [GammaCube_FT, shitIsFucked] = dRhoMaster( rho, rho_FT, ...
  flags, interObj,  systemObj, diffObj, particleObj, cosPhi3, sinPhi3,t )
% Initialize
GammaCube_FT = 0;
shitIsFucked = 0;
% Interactions
% mayers
if interObj.hardId == 1 % mayers
  if interObj.typeId == 1; % rods
    muExFt = muExCalcVc2Ft(rho_FT, interObj.FmFt,systemObj,interObj.muMayerScale);
    GammaExCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
    GammaCube_FT = GammaCube_FT + GammaExCube_FT;
  end
  if interObj.typeId == 2; % disks
    muExFt = muExCalcVc2Ft(rho_FT(:,:,interObj.k3ind0 ), interObj.FmFt,systemObj,interObj.muMayerScale);
    GammaExCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
    GammaCube_FT = GammaCube_FT + GammaExCube_FT;
  end
end
% spt
if interObj.hardId == 2 % spt
  if interObj.typeId == 2 % disks
    if systemObj.n3 == 1
      nu = interObj.sptScale .* rho;
    else
      nu = interObj.sptScale .* ifftn( ifftshift( rho_FT(:,:,interObj.k30) ) );
    end
    if any(nu <  1)
      error('Density is too high!')
      fprintf('Density is too high!');
      shitIsFucked = 1;
      muExFt = 0;
    else
      [muExFt] = muExDisksSPT(nu);
    end
    GammaExCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
    GammaCube_FT = GammaCube_FT + GammaExCube_FT;
  end
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
