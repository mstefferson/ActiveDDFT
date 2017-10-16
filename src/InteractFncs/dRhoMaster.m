% Handles all the dRho contributions that are not in Lop
function [gammaCubeFt, shitIsFucked, whatBroke] = dRhoMaster( rho, rho_FT, ...
  interObj,  systemObj, diffObj, polarDrive, noise, densityDepDr, prop )
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
    interObj.dv1Flag, interObj.dv2Flag, interObj.dv3Flag, ...
    interObj.intInd1, interObj.intInd2, interObj.intInd3);
  dVmaster.dx1 = dMu.dx1  + dVmaster.dx1;
  dVmaster.dx2 = dMu.dx2  + dVmaster.dx2;
  dVmaster.dx3 = dMu.dx3  + dVmaster.dx3;
  [gammaCubeFt, j3Ex] = dRhoIntCalcMu( rho, dVmaster, systemObj, diffObj, interObj);
end
%% Debug
debug = 1;
if debug
  gammaIntFt = gammaCubeFt;
end
%%
% driving
if polarDrive.Flag
  gammaDrCubeFt = polarDrive.calcDrho( rho );
  gammaCubeFt = gammaCubeFt + gammaDrCubeFt;
end
% noise
if noise.Flag
  gammaNoiseFt = noise.calcDrho( rho );
  gammaCubeFt = gammaCubeFt + gammaNoiseFt;
end
% density dep diffusion
if densityDepDr.Flag
  gammaRotDiffFt = densityDepDr.calcDrho( rho, rho_FT, j3Ex );
  gammaCubeFt = gammaCubeFt + gammaRotDiffFt;
end

if debug
  gammaDiffFt = prop .* rho_FT;
  gammaDiff = real( ifftn(ifftshift( gammaDiffFt ) ) );
  gammaInt = real( ifftn(ifftshift( gammaIntFt ) ) );
  gammaRotDiff = real( ifftn(ifftshift( gammaRotDiffFt ) ) );
  gammaDrive = real( ifftn(ifftshift( gammaDrCubeFt  ) ) );
  %%
  %   gammaRotDiff2Plot = sum( gammaRotDiff, 3);
  %   gammaInt2Plot = sum( gammaInt, 3);
  %   gammaDrive2Plot = sum( gammaDrive, 3);
  indWant = 17;
  gammaDiff2Plot = gammaDiff(:,:,indWant);
  gammaRotDiff2Plot = gammaRotDiff(:,:,indWant);
  gammaInt2Plot = gammaInt(:,:,indWant);
  gammaDrive2Plot = gammaDrive(:,:,indWant);
  rho2plot = rho(:,:,indWant);
  subplot(2,2,1)
  pcolor( gammaRotDiff2Plot ); colorbar;
  title('Rot Diffusion')
  subplot(2,2,2)
  pcolor( gammaInt2Plot ); colorbar;
  title('Interaction')
  subplot(2,2,3)
  pcolor( gammaDiff2Plot ); colorbar;
  title('Diffusion')
  subplot(2,2,4)
%   pcolor( gammaDrive2Plot ); colorbar;
%   title('Drive')
  pcolor( rho2plot ); colorbar;
  title('rho')
  
  % check conservation
  if 0
    normRotDiff = sum( gammaRotDiff(:) )
    normInt = sum( gammaInt(:) )
    normDrive = sum( gammaDrCubeFt(:) )
    
    ind1 = systemObj.n1/2 + 1;
    ind2 = systemObj.n2/2 + 1;
    ind3 = systemObj.n3/2 + 1;
    
    gammaIntFt( ind1, ind2, ind3 )
    gammaRotDiffFt( ind1, ind2, ind3 )
    gammaDrCubeFt( ind1, ind2, ind3 )
  end
  %%
  keyboard
end
