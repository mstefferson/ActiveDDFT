% [interObj] = interObjMaker( particleObj, systemObj )
% Creates an interaction object that is used by interaction functions.
% Includes interaction types and calculates necessary functions.
%
% Id key:
% typeId = 1, 2, 3 rods, disks, spheres
% hardId = 1, 2, 3 mayer, spt, fmt
% softId = 1 soft shoulder
function [interObj] = interObjMaker( particleObj, systemObj, gridObj )
% Store interactions in interobj
interObj.anyInter = 0;
% save k3 index
interObj.k3ind0 = floor( systemObj.n3 / 2 ) + 1;
% flags
interObj.dv1Flag = 0;
interObj.dv2Flag = 0;
interObj.dv3Flag = 0;
% Short range interactions
if isempty(particleObj.interHb)
  interObj.hardFlag = 0;
  interObj.hardId = 0;
  interObj.hard = 'none';
  fprintf('No short interactions\n');
else
  interObj.hardFlag = 1;
  interObj.hard = particleObj.interHb;
  interObj.anyInter = 1;
  % rods
  if  strcmp( particleObj.type, 'rods' )
    interObj.typeId = 1;
    % flags/inds
    interObj.dv1Flag = 1;
    interObj.dv2Flag = 1;
    interObj.dv3Flag = 1;
    interObj.srInd1 = 1:systemObj.n1;
    interObj.srInd2 = 1:systemObj.n2;
    interObj.srInd3 = 1:systemObj.n3;
    % mayer
    if strcmp( interObj.hard, 'mayer' )
      fprintf('%s hard %s\n', interObj.hard, particleObj.type);
      interObj.hardSpec = 'rodsMayer';
      interObj.hardId = 1;
      % grab lap frame mayer function
      [~,interObj.FmFt] = mayerFncHrLabFrame(...
        systemObj.n1, systemObj.n2, systemObj.n3, ...
        systemObj.l1, systemObj.l2, particleObj.lMaj);
      interObj.muMayerScale = ( systemObj.tmp * systemObj.l3 * systemObj.l1 * systemObj.l2) ./ ...
        (systemObj.n1 * systemObj.n2 * systemObj.n3 ^ 2);
      interObj.muMayerInds = 1:systemObj.n3;
      interObj.muMayerMinusInds = [1 systemObj.n3:-1:2];
    else
      fprintf('Cannot find hard rod interactions\n')
    end
    % disks
  elseif strcmp( particleObj.type, 'disks' )
    % flags/inds
    interObj.dv1Flag = 1;
    interObj.dv2Flag = 1;
    interObj.typeId = 2;
    interObj.srInd1 = 1:systemObj.n1;
    interObj.srInd2 = 1:systemObj.n2;
    interObj.srInd3 = floor(systemObj.n3 / 2) + 1;
    % mayer
    if strcmp( interObj.hard, 'mayer' )
      interObj.hardSpec = 'disksMayer';
      interObj.hardId = 1;
      % mayer only a function of x1, x2 for disks
      [~,interObj.FmFt] = mayerFncHd(...
        systemObj.n1, systemObj.n2, ...
        systemObj.l1, systemObj.l2, particleObj.lMaj) ;
      interObj.muMayerScale = (systemObj.l3 * systemObj.l1 * systemObj.l2) ./ ...
        (systemObj.n1 * systemObj.n2 * systemObj.n3);
      fprintf('%s hard %s\n', interObj.hard, particleObj.type);
      % scaled particle theory
    elseif strcmp( interObj.hard, 'spt' )
      interObj.hardSpec = 'disksSpt';
      interObj.hardId = 2;
      % packing fraction factor mulitplied by FT factors (for concentration)
      interObj.sptScale = systemObj.l3 ./ systemObj.n3 .* particleObj.b;
      interObj.k30 = gridObj.k3ind0;
      fprintf('%s hard %s\n', interObj.hard, particleObj.type);
    else
      fprintf('Cannot find hard disk interactions\n')
    end
    % spheres
  elseif strcmp( particleObj.type, 'spheres' )
    interObj.typeId = 3;
    % flags/inds
    interObj.dv1Flag = 1;
    interObj.dv2Flag = 1;
    interObj.dv3Flag = 1;
    interObj.srInd1 = 1:systemObj.n1;
    interObj.srInd2 = 1:systemObj.n2;
    interObj.srInd3 = 1:systemObj.n3;
    % mayer
    if strcmp( interObj.hard, 'mayer' )
      interObj.hardSpec = 'spheresMayer';
      interObj.hardId = 1;
      [~,interObj.FmFt] = mayerFncHs(...
        systemObj.n1, systemObj.n2, systemObj.n3,...
        systemObj.l1, systemObj.l2, systemObj.l3, particleObj.lMaj) ;
      interObj.muExScale = (systemObj.l3 * systemObj.l1 * systemObj.l2) ./ ...
        (systemObj.n1 * systemObj.n2 * systemObj.n3);
      fprintf('%s hard %s\n', interObj.hard, particleObj.type);
    else
      fprintf('Cannot find hard sphere interactions\n')
    end
  else
    fprintf('Cannot find particle type\n');
    error('Cannot find particle type\n');
  end % particle type
end % short range interactions
% Long range interactions
if isempty(particleObj.interactLrV)
  interObj.longFlag = 0;
  interObj.longId = 0;
  interObj.long = 'none';
  fprintf('No long range interactions\n');
else
  % store some things
  interObj.longFlag = 1;
  numPotentials = length( particleObj.interactLrV );
  interObj.anyInter = 1;
  fprintf('Long interactions\n');
  % Scale
  interObj.muMfScale = (systemObj.l3 * systemObj.l1 * systemObj.l2) ./ ...
    (systemObj.n1 * systemObj.n2 * systemObj.n3);
  % add the potenials, v
  v = 0;
  for ii = 1:numPotentials
    currPot = particleObj.interactLrV{ii};
    fprintf('Interaction %s along dim \n', currPot{1});
    % soft shoulder 2D
    if strcmp( currPot{1}, 'ss2d' )
      % build potential
      vTemp = SoftShoulderClass( currPot{1}, ...
        currPot{2}(1), currPot{2}(2), currPot{2}(3),currPot{2}(4),...
        systemObj.n1, systemObj.l1, systemObj.n2, systemObj.l2 );
      interObj.interactLrV{ii} = vTemp;
      v = v + reshape( vTemp.V, vTemp.ReshapeInds );
    end
    % polar align 2D
    if strcmp( currPot{1}, 'pa2d' )
      % build potential
      vTemp = PolarAlignClass( currPot{1}, ...
        currPot{2}(1), systemObj.n3, systemObj.l3);
      interObj.interactLrV{ii} = vTemp;
      v = v + reshape( vTemp.V, vTemp.ReshapeInds );
    end
    % polar align 2D
    if strcmp( currPot{1}, 'pag2d' )
      % build potential
      vTemp = PolarAlignGaussClass( currPot{1}, ...
        currPot{2}(1), currPot{2}(2),...
        systemObj.n1, systemObj.n2, systemObj.n3, ...
        systemObj.l1, systemObj.l2, systemObj.l3);
      interObj.interactLrV{ii} = vTemp;
      v = v + vTemp.V;
    end
    % decaying exponential 2D
    if strcmp( currPot{1}, 'de2d' )
      % build potential
      vTemp = DecayExpClass( currPot{1}, ...
        currPot{2}(1), currPot{2}(2), ...
        systemObj.n1, systemObj.n2, ...
        systemObj.l1, systemObj.l2 );
      interObj.interactLrV{ii} = vTemp;
      v = v + reshape( vTemp.V, vTemp.ReshapeInds );
    end
  end % loop over potentials
  % get vFt and find shape
  [vn1, vn2, vn3] = size( v );
  if vn1 == 1
    interObj.lrInd1 = floor(systemObj.n1/2) + 1;
  else
    interObj.lrInd1 = 1:systemObj.n1;
    interObj.dv1Flag = 1;
  end
  if vn2 == 1
    interObj.lrInd2 = floor(systemObj.n2/2) + 1;
  else
    interObj.lrInd2 = 1:systemObj.n2;
    interObj.dv2Flag = 1;
  end
  if vn3 == 1
    interObj.lrInd3 = floor(systemObj.n3/2) + 1;
  else
    interObj.lrInd3 = 1:systemObj.n3;
    interObj.dv3Flag = 1;
  end
  % store it
  interObj.vInt = v;
  interObj.vIntFt = fftshift( fftn( v ) );
end % long range interaction
% External potentional
if isempty(particleObj.externalV)
  interObj.extFlag = 0;
  interObj.ext = 'none';
  fprintf('No external potential\n');
else
  interObj.extFlag = 1;
  interObj.anyInter = 1;
  numexternalVs = length( particleObj.externalV );
  interObj.externalV = cell(numexternalVs,1);
  nVec = [systemObj.n1 systemObj.n2 systemObj.n3];
  % reinitialize
  v = 0;
  dV.dx1 = 0;
  dV.dx2 = 0;
  dV.dx3 = 0;
  for ii = 1:numexternalVs
    currPot = particleObj.externalV{ii};
    fprintf('External %s along dim %d\n', currPot{1}, currPot{2}(1) );
    if strcmp( currPot{1}, 'linV' )
      % run a check just in case
      if nVec( currPot{2}(1) ) > 1
        vTemp = LinearVClass( currPot{1}, currPot{2}(1), currPot{2}(2), gridObj.x1 );
        interObj.externalV{ii} = vTemp;
        v = v + vTemp.VvReshape;
        dV.dx1 = dV.dx1 + vTemp.DvDx1;
        dV.dx2 = dV.dx2 + vTemp.DvDx2;
        dV.dx3 = dV.dx3 + vTemp.DvDx3;
      end
    elseif strcmp( currPot{1}, 'quadV' )
      % run a check just in case
      if nVec( currPot{2}(1) ) > 1
        vTemp = QuadVClass( currPot{1}, currPot{2}(1), currPot{2}(2), gridObj.x1 );
        interObj.externalV{ii} = vTemp;
        v = v + vTemp.VvReshape;
        dV.dx1 = dV.dx1 + vTemp.DvDx1;
        dV.dx2 = dV.dx2 + vTemp.DvDx2;
        dV.dx3 = dV.dx3 + vTemp.DvDx3;
      end
    else
      fprintf('Not building potential, dimension not available\n' );
    end
  end % for loop
  interObj.dVExt = dV;
  if currPot{2}(1)  == 1
    interObj.dv1Flag = 1;
  end
  if currPot{2}(1)  == 2
    interObj.dv2Flag = 1;
  end
  if currPot{2}(1)  == 3
    interObj.dv3Flag = 1;
  end
end % external