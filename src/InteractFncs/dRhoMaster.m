% Handles all the dRho contributions that are not in Lop
function [GammaCube_FT, shitIsFucked] = dRhoMaster( rho, rho_FT, ...
  flags, interObj,  systemObj, diffObj, particleObj, cosPhi3, sinPhi3 )
% Initialize
debug  = 1;
GammaCube_FT = 0;
shitIsFucked = 0;
% Interactions
% short range
% mayers
if interObj.hardId == 1 % mayers
  if interObj.typeId == 1 % rods
    
    if debug
      muExFt = muExCalcVc2Ft(rho_FT, interObj.FmFt,systemObj,interObj.muMayerScale*systemObj.n3);
      muExFt2 = muExCalcMayerLF(rho_FT, interObj.FmFt2, systemObj,...
        interObj.muMayerScale, interObj.muMayerInds, interObj.muMayerMinusInds);
      muEx = real( ifftn( ifftshift( muExFt ) ) );
      muEx2 = real( ifftn( ifftshift( muExFt2 ) ) );
      figure()
      if 1
        ii = 1;
        subplot(3,1,1); imagesc( muEx(:,:,ii) ); colorbar;
        subplot(3,1,2); imagesc( muEx2(:,:,ii) ); colorbar;
        subplot(3,1,3); imagesc( muEx(:,:,ii) - muEx2(:,:,ii) ); colorbar;
      end
    else
%       muExFt = muExCalcVc2Ft(rho_FT, interObj.FmFt,systemObj,interObj.muMayerScale*systemObj.n3);
      muExFt = muExCalcMayerLF(rho_FT, interObj.FmFt2, systemObj,...
        interObj.muMayerScale, interObj.muMayerInds, interObj.muMayerMinusInds);
%       muEx = real( ifftn( ifftshift( muExFt ) ) );
%       muEx2 = real( ifftn( ifftshift( muExFt2 ) ) );
    end
    GammaExCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
    GammaCube_FT = GammaCube_FT + GammaExCube_FT;
    keyboard
  end
  if interObj.typeId == 2 % disks
    muExFt = muExCalcVc2Ft( rho_FT(:,:,interObj.k3ind0), interObj.FmFt, systemObj, interObj.muMayerScale );
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
    % check if too high or too low
    if any(nu >  1)
      error('Density is too high!')
      fprintf('Density is too high!');
      shitIsFucked = 1;
      muExFt = 0;
    elseif any(nu <  0)
      fprintf('Density is negative!');
      error('Density is negative!');
      shitIsFucked = 1;
      muExFt = 0;
    else
      [muExFt] = muExDisksSPT(nu);
    end
    GammaExCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
    GammaCube_FT = GammaCube_FT + GammaExCube_FT;
  end
end
% long range
if interObj.longId == 1 % mean field
  [muExFt] =  muExCalcMfFt( rho_FT(:,:,interObj.k3ind0), interObj.vFt, systemObj, interObj.muMfScale );
  GammaExCube_FT = dRhoIntCalcMu( rho, muExFt, systemObj, diffObj );
  GammaCube_FT = GammaCube_FT + GammaExCube_FT;
end
% driving
if flags.Drive && flags.DiagLop
  GammaDrCube_FT  = ...
    dRhoDriveCalcFtId(rho,particleObj.vD,...
    cosPhi3, sinPhi3,diffObj.ik1rep3,diffObj.ik2rep3);
  GammaCube_FT = GammaCube_FT + GammaDrCube_FT;
end
% external
if interObj.extFlag
  GammaPotCube_FT = 0;
  GammaCube_FT = GammaCube_FT + GammaPotCube_FT;
end
