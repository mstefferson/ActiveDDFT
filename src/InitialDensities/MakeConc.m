% Uses the initial density indicator to choose the correct intial density
% subroutine.
%
function [rho] = MakeConc(systemObj,particleObj,rhoInit,gridObj)
% go through conditions
if strcmp( rhoInit.intCond, 'iso' ) % iso
  fprintf('IC: isotropic \n' );
  % Initial distribution
  [rho] = IntDenCalcIso(systemObj);
elseif strcmp( rhoInit.intCond, 'eq' ) % eq
  fprintf('IC: equilbrium \n' );
  % Initial distribution
  [rho] = IntDenCalcEq(systemObj, particleObj, rhoInit);
elseif strcmp( rhoInit.intCond, 'nem' ) % nematic
  if systemObj.n3 > 1
    fprintf('IC: nematic \n' );
    % Initial distribution
    [rho] = IntDenCalcNem(systemObj,gridObj.x3,rhoInit.shiftAngle);
  else
    fprintf('IC: nematic requested, but no angular grid! Going iso\n' );
    [rho] = IntDenCalcIso(systemObj);
  end
elseif strcmp( rhoInit.intCond, 'load' ) % load
  fprintf('IC: loading %s \n', rhoInit.loadName );
  % Initial distribution
  [rho] = IntDenCalcLoaded2Drot( [rhoInit.pathName rhoInit.loadName],systemObj);
elseif strcmp( rhoInit.intCond, 'delP' ) % delta function in polar
   fprintf('IC: delta function in polar\n' );
  % delta function in polar order
  [rho] = deltaPolarIc(systemObj);
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
elseif strcmp( rhoInit.intCond, 'gauss' )
  rho = gaussCalc( systemObj, rhoInit, gridObj );
elseif strcmp( rhoInit.intCond, 'lorenz' )
  rho = lorenzianCalc( systemObj, rhoInit, gridObj );
else
  fprintf('Cannot find desired IC. Setting to iso');
  % Initial distribution
  [rho] = IntDenCalcIso(systemObj);
end

% Perturb it
for ii = 1:rhoInit.numPerturb
  perturbTemp = rhoInit.perturb{1};
  end
  % Plane wave perturb
  if strcmp( perturbTemp, 'pw' )
    rho = PwPerturbFT( rho, systemObj, perturbTemp );
  % Lorenzian perturb
  elseif strcmp( perturbTemp, 'lorenz' )
    perturb = lorenzianCalc( systemObj, perturbTemp, gridObj );
    rho = rho + perturb;
  % Gaussian perturb
  elseif strcmp( perturbTemp, 'lorenz' )
    perturb = gaussCalc( systemObj, perturbTemp, gridObj );
    rho = rho + perturb;
  else
    printf('Cannot find perturbation. Not perturbing\n')
  end
end
% Fix negative rho if that happened.
[rho] = FixNegDenFnc(rho);
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
