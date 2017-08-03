% Handles all the dRho contributions that are not in Lop
function [GammaCube_FT, shitIsFucked, whatBroke] = dRhoMaster( rho, rho_FT, ...
  flags, interObj,  systemObj, diffObj, particleObj, cosPhi3, sinPhi3 )
% Initialize
GammaCube_FT = 0;
shitIsFucked = 0;
whatBroke = [];
% Interactions: short range
% mayers
if interObj.hardId == 1 % mayers
  if interObj.typeId == 1 % rods
    % calculate excess chemical potential
    muExFt = muExCalcMayerLF(rho_FT, interObj.FmFt, systemObj,...
      interObj.muMayerScale, interObj.muMayerInds, interObj.muMayerMinusInds);
    GammaCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
  end
  if interObj.typeId == 2 % disks
    muExFt = muExCalcVc2Ft( rho_FT(:,:,interObj.k3ind0), interObj.FmFt, systemObj, interObj.muMayerScale );
    GammaCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
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
    % check if too high or too low
    if any(nu(:) >  1)
      fprintf('Density is too high!');
      shitIsFucked = 1;
      muExFt = 0;
      whatBroke = 'density too high from spt';
    elseif any(nu(:) <  0)
      fprintf('Density is negative!');
      shitIsFucked = 1;
      whatBroke = 'density negative';
      muExFt = 0;
    else
      [muExFt] = muExDisksSPT(nu);
    end
    GammaCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
  end
end
% Interactions: long range
if interObj.longFlag  % mean field
  [muExFt] =  muExCalcMfFt( rho_FT(interObj.ind1,interObj.ind2,interObj.ind3),...
  interObj.vFt, systemObj, interObj.muMfScale );
  GammaExCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj, interObj);
  GammaCube_FT = GammaCube_FT + GammaExCube_FT;
end
% driving
if flags.Drive && flags.DiagLop
  GammaDrCube_FT  = ...
    dRhoDriveCalcFtId(rho,particleObj.vD,...
    cosPhi3, sinPhi3,diffObj.ik1rep3,diffObj.ik2rep3);
  GammaCube_FT = GammaCube_FT + GammaDrCube_FT;
end
% external potential
if interObj.extFlag
  GammaPotCube_FT = 0;
  GammaCube_FT = GammaCube_FT + GammaPotCube_FT;
end
