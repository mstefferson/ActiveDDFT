% Handles all the dRho contributions that are not in Lop
function [gammaCubeFt, shitIsFucked, whatBroke] = dRhoMaster( rho, rho_FT, ...
  interObj,  systemObj, diffObj, polarDrive, noise, densityDepDr, lop, prop, dt, gamProp)
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
  [gammaRotDiffFt, gammaRotDiffFtDiff, gammaRotDiffFtInt] = densityDepDr.calcDrho( rho, rho_FT, j3Ex );
  gammaCubeFt = gammaCubeFt + gammaRotDiffFt;
end

if debug
  %gammaDiffFt = ( prop - 1 ) .* rho_FT;
  gammaDiffFt = ( lop ) .* rho_FT;
  
  % calc in real space
  gammaDiff = real( ifftn(ifftshift( gammaDiffFt ) ) );
  gammaInt = real( ifftn(ifftshift( gammaIntFt ) ) );
  gammaRotDiff = real( ifftn(ifftshift( gammaRotDiffFt ) ) );
  gammaRotDiffDiff = real( ifftn(ifftshift( gammaRotDiffFtDiff ) ) );
  gammaRotDiffInt = real( ifftn(ifftshift( gammaRotDiffFtInt ) ) );

  % scale is
  nlProp = 1;
  gammaInt = nlProp .* gammaInt;
  gammaRotDiff = nlProp .* gammaRotDiff;
  
  if 0
    rho2plot = sum( rho, 3);
    gammaDiff2Plot = sum( gammaDiff, 3);
    gammaRotDiff2Plot = sum( gammaRotDiff, 3);
    gammaRotDiffDiff2Plot = sum( gammaRotDiff, 3);
    gammaRotDiffInt2Plot = sum( gammaRotDiff, 3);
    gammaInt2Plot = sum( gammaInt, 3);
  end
  
  if 1
    indWant = 17;
    gammaDiff2Plot = gammaDiff(:,:,indWant);
    gammaRotDiff2Plot = gammaRotDiff(:,:,indWant);
    gammaRotDiffDiff2Plot = gammaRotDiffDiff(:,:,indWant);
    gammaRotDiffInt2Plot = gammaRotDiffInt(:,:,indWant);
    gammaInt2Plot = gammaInt(:,:,indWant);
    rho2plot = rho(:,:,indWant);
  end
   
  % plot different contributions
  figure
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
  pcolor( rho2plot ); colorbar;
  title('rho')
  
  % compate rotational diffusion parts
  figure
  subplot(2,2,1)
  pcolor( gammaRotDiff2Plot ); colorbar;
  title('Rot Diff Total')
  subplot(2,2,2)
  pcolor( gammaRotDiffDiff2Plot ); colorbar;
  title('Rot Diff Diff ')
  subplot(2,2,3)
  pcolor( gammaRotDiffInt2Plot ); colorbar;
  title('Rot Diff Int ')
  subplot(2,2,4)
  pcolor( gammaDiff2Plot ); colorbar;
  title('Diffusion')
  
  figure()
  hold
  n = systemObj.n3;
  phi = 1:n;
  ind1 = 1;
  ind2 = 1;

  plot( phi, reshape( gammaDiff(ind1,ind2,:), [1 n] ) );
  plot( phi, reshape( gammaInt(ind1,ind2,:), [1 n] ) );
  plot( phi, reshape( gammaRotDiffDiff(ind1,ind2,:), [1 n] ) );
  plot( phi, reshape( gammaRotDiffInt(ind1,ind2,:), [1 n] ) );
  legend( 'pure diff', 'gamma Int', 'rot diff diff',...
    'rot diff int' )
  
  figure()
  hold
   plot( phi, reshape( gammaDiff(ind1,ind2,:) + gammaRotDiffDiff(ind1,ind2,:),...
     [1 n] ) );
  plot( phi, reshape( gammaInt(ind1,ind2,:) + gammaRotDiffInt(ind1,ind2,:), [1 n] ) );
  legend( 'diff', 'Int')
  
  % check flux
  jDiff = -densityDepDr.Dr0 .* ...
    real( ifftn( ifftshift( densityDepDr.Ik3 .* rho_FT ) ) );
  jDiffRot = -real( ifftn( ifftshift( ...
    densityDepDr.DrNl .* densityDepDr.Ik3 .* rho_FT ) ) );
  figure()
  hold
  plot( phi, reshape( jDiff(ind1,ind2,:), [1 n] ) );
  plot( phi, reshape( jDiffRot(ind1,ind2,:), [1 n] ) );
  
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
