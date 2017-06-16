% IntDenCalcEq(systemObj, particleObj rhoInit)
function [rho] = IntDenCalcEq(systemObj, particleObj, rhoInit)
% make equilbrium concentration based on interaction
if strcmp( particleObj.interHb, 'mayer')
  fprintf('Choosing hardrod IN equilbrium distribution\n');
  % if bc is too close to 1.5, errors arise. Fix this here.
  if 1.499 < systemObj.bc && systemObj.bc < 1.501
    systemObjTemp = systemObj;
    systemObjTemp.bc = 1.502;
    % Initial distribution
    [rho] = IntDenCalcEqHr(systemObjTemp,rhoInit);
  else
    % Initial distribution
    [rho] = IntDenCalcEqHr(systemObj,rhoInit);
  end
elseif strcmp( particleObj.interLr, 'softshoulder')
  % find equilbrium crystal
  paramVec = [systemObj.n1 systemObj.n2 systemObj.l1 systemObj.l2...
    particleObj.lrEs1 particleObj.lrEs2 particleObj.lrLs1 particleObj.lrLs2...
    systemObj.c];
  disp = dispersionSoftShoulder( paramVec );
  if systemObj.n1 == systemObj.n2 && systemObj.l1 == systemObj.l2 
    if strcmp( disp.phase, 'crystal A' )
      a = 2.26;
      rho = hexCrystal(systemObj.l1, systemObj.n1, systemObj.numPart, ...
        a, rhoInit.crystalLattice(2) );
    elseif strcmp( disp.phase, 'crystal B' )
      a = 1.21;
      rho = hexCrystal(systemObj.l1, systemObj.n1, systemObj.numPart, ...
        a, rhoInit.crystalLattice(2) );
    else
      a = Inf;
      [rho] = IntDenCalcIso(systemObj);
    end
    fprintf('Choosing soft shoulder %s with %.2f lattice spacing\n', ...
      disp.phase, a );
    % rep it
    rho =  1 / systemObj.l3 * repmat( rho, [1,1,systemObj.n3] );
  else
    fprintf('Not making crystal, box not symmetric\n');
    error('Box must be symmetric');
    rho = 0;
  end
else
  fprintf('Cannot find equilbrium distribution based on interaction, making iso\n');
  [rho] = IntDenCalcIso(systemObj);
end
