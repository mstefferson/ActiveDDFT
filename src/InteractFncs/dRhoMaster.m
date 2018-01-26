% Handles all the dRho contributions that are not in Lop
function [gammaCubeFt, shitIsFucked, whatBroke] = dRhoMaster( rho, rho_FT, ...
  interObj,  systemObj, diffObj, polarDrive, noise, densityDepDiff )
% Initialize
gammaCubeFt = 0;
shitIsFucked = 0;
calcGamma = 0;
whatBroke = [];
muExFt = 0;
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
end % hard
% Interactions: long range
if interObj.longFlag  % mean field
  calcGamma = 1;
  [muExFtTemp] =  muExCalcPDirCorrFt( rho_FT(interObj.lrInd1,interObj.lrInd2,interObj.lrInd3),...
    interObj.c2Ft, systemObj, interObj.muMfScale );
  muExFt = muExFt + muExFtTemp;
end
% calculate gamma once
if calcGamma
  dMu = dVCalc(muExFt, diffObj, ...
    interObj.dv1IntFlag, interObj.dv2IntFlag, interObj.dv3IntFlag, ...
    interObj.intInd1, interObj.intInd2, interObj.intInd3);
  dVmaster.dx1 = dMu.dx1  + dVmaster.dx1;
  dVmaster.dx2 = dMu.dx2  + dVmaster.dx2;
  dVmaster.dx3 = dMu.dx3  + dVmaster.dx3;
  [gammaCubeFt, iotaEx1, iotaEx2, iotaEx3] = dRhoIntCalcMu( ...
    rho, dVmaster, systemObj, diffObj, interObj);
else
  iotaEx1 = 0; 
  iotaEx2 = 0; 
  iotaEx3 = 0;
end
% driving
if polarDrive.Flag
  gammaDrCubeFt = polarDrive.calcDrho( rho );
  gammaCubeFt = gammaCubeFt + gammaDrCubeFt;
  iotaEx1 = polarDrive.Iota1 + iotaEx1;
  iotaEx2 = polarDrive.Iota2 + iotaEx2;
end
% noise
if noise.Flag
  gammaNoiseFt = noise.calcDrho( rho );
  gammaCubeFt = gammaCubeFt + gammaNoiseFt;
end
% density dep diffusion
if densityDepDiff.Flag
  gammaDiffFt = densityDepDiff.calcDrho( rho_FT,...
    {iotaEx1, iotaEx2, iotaEx3} );
  gammaCubeFt = gammaCubeFt + gammaDiffFt;
end
