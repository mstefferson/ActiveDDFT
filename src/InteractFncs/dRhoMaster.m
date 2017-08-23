% Handles all the dRho contributions that are not in Lop
function [GammaCube_FT, shitIsFucked, whatBroke] = dRhoMaster( rho, rho_FT, ...
  flags, interObj,  systemObj, diffObj, particleObj, cosPhi3, sinPhi3 )
% Initialize
GammaCube_FT = 0;
shitIsFucked = 0;
calcGamma = 0;
whatBroke = [];
% Calculate dV master
if interObj.extFlag
  calcGamma = 1;
  dVmaster.dx1 = interObj.dVExt.dx1;
  dVmaster.dx2 = interObj.dVExt.dx2;
  dVmaster.dx3 = interObj.dVExt.dx3;
else
  dVmaster.dx1 = 0;
  dVmaster.dx2 = 0;
  dVmaster.dx3 = 0;
end
% Interactions: short range
if interObj.hardFlag 
  calcGamma = 1;
  if interObj.hardId == 1 % mayers
    % calculate excess chemical potential
    muExFt = muExCalcMayerLF(rho_FT(interObj.srInd1,interObj.srInd2,interObj.srInd3),...
    interObj.FmFt, systemObj,...
      interObj.muMayerScale, interObj.muMayerInds, interObj.muMayerMinusInds);
  elseif interObj.hardId == 2 % spt
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
    end
  else 
    fprintf('Cannot find hardID\n');
    error('Cannot find hardID');
  end
  dMu = dVCalc(muExFt, systemObj, diffObj, ...
    interObj.srInd1, interObj.srInd2, interObj.srInd3);
  dVmaster.dx1 = dMu.dx1  + dVmaster.dx1;
  dVmaster.dx2 = dMu.dx2  + dVmaster.dx2;
  dVmaster.dx3 = dMu.dx2  + dVmaster.dx3;
end % hard
% Interactions: long range
if interObj.longFlag  % mean field
  calcGamma = 1;
  [muExFt] =  muExCalcPDirCorrFt( rho_FT(interObj.lrInd1,interObj.lrInd2,interObj.lrInd3),...
  interObj.c2Ft, systemObj, interObj.muMfScale );
  dMu = dVCalc(muExFt, systemObj, diffObj, ...
    interObj.lrInd1, interObj.lrInd2, interObj.lrInd3);
  dVmaster.dx1 = dMu.dx1  + dVmaster.dx1;
  dVmaster.dx2 = dMu.dx2  + dVmaster.dx2;
  dVmaster.dx3 = dMu.dx3  + dVmaster.dx3;
end
% calculate gamma once
if calcGamma
  GammaCube_FT = dRhoIntCalcMu( rho, dVmaster, systemObj, diffObj, interObj);
end
% driving
if flags.Drive && flags.DiagLop
  GammaDrCube_FT  = ...
    dRhoDriveCalcFtId(rho,particleObj.vD,...
    cosPhi3, sinPhi3,diffObj.ik1rep3,diffObj.ik2rep3);
  GammaCube_FT = GammaCube_FT + GammaDrCube_FT;
end
