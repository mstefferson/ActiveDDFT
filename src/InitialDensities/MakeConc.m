% Uses the initial density indicator to choose the correct intial density
% subroutine.
%
function [rho] = MakeConc(systemObj,particleObj,rhoInit,gridObj)
% go through conditions
if strcmp( rhoInit.intCond, 'iso' ) % iso
  fprintf('IC: isotropic \n' );
  % Initial distribution
  [rho] = IntDenCalcIso(systemObj);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif strcmp( rhoInit.intCond, 'eq' ) % eq
  fprintf('IC: equilbrium \n' );
  % Initial distribution
  [rho] = IntDenCalcEq(systemObj, particleObj, rhoInit);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif strcmp( rhoInit.intCond, 'nem' ) % nematic
  if systemObj.n3 > 1
    fprintf('IC: nematic \n' );
    % Initial distribution
    [rho] = IntDenCalcNem(systemObj,gridObj.x3,rhoInit.shiftAngle);
  else
    fprintf('IC: nematic requested, but no angular grid! Going iso\n' );
    [rho] = IntDenCalcIso(systemObj);
  end
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif strcmp( rhoInit.intCond, 'load' ) % load
  fprintf('IC: loading %s \n', rhoInit.LoadName );
  % Initial distribution
  [rho] = IntDenCalcLoaded2Drot(rhoInit.LoadName,systemObj);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif strcmp( rhoInit.intCond, 'delP' ) % delta function in polar
   fprintf('IC: delta function in polar\n' );
  % delta function in polar order
  [rho] = deltaPolarIc(systemObj);
  % Perturb it
  [rho] = PwPerturbFT(rho,systemObj,rhoInit);
elseif strcmp( rhoInit.intCond, 'crys' ) % hexagon lattice
  fprintf('IC: hexagonal lattice with spacing %.2f\n', ...
    rhoInit.crystalLattice(1) );
  if systemObj.n1 == systemObj.n2 && systemObj.l1 == systemObj.l2 
    rho = hexCrystal(systemObj.l1, systemObj.n1, systemObj.numPart, ...
      rhoInit.crystalLattice(1), rhoInit.crystalLattice(2) );
    % rep it
    rho =  1 / systemObj.l3 * repmat( rho, [1,1,systemObj.n3] );
  else
    fprintf('Error: box must be symmetric\n');
    error('Box must be symmetric');
  end
else
  fprintf('Error not written');
  error('Error not written');
end
% perturb it
elseif strcmp( rhoInit.intCond, 'lorenz' ) % lorenzian lattice
  rho = lorenzianIc(systemObj,gridObj);
elseif rhoInit.IntCond == 4 || rhoInit.IntCond == 5
  if rhoInit.IntCond == 4
    fprintf('IC: gaussian perp with homogenous concentration \n' );
  else
    fprintf('IC: gaussian perp with inhomogenous concentration \n' );
  end
  % Build eq first
  [rho] = IntDenCalcEq(systemObj, particleObj, rhoInit);
  % Gaussian perturb
  [rho] = polarPerturbGauss( rho, systemObj, rhoInit, gridObj );
end
% Renormalize out here
% Integrate first along the depth of matrix w.r.t theta, then across the
% columns w.r.t x, then down the rows w.r.t. y
if systemObj.n3 > 1
  int1 = trapz_periodic(gridObj.x3,rho,3);
else
  int1 = rho;
end
CurrentNorm = trapz_periodic(gridObj.x1,...
  trapz_periodic(gridObj.x2,int1,2),1);
rho = rho .* systemObj.numPart ./ CurrentNorm;
end
