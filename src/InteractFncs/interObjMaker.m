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
interObj.hardFlag = 0;
interObj.longFlag = 0;
interObj.extFlag = 0;
% save k3 index
interObj.k3ind0 = floor( systemObj.n3 / 2 ) + 1;
% flags, any int
interObj.dv1Flag = 0;
interObj.dv2Flag = 0;
interObj.dv3Flag = 0;
% flags, any int
interObj.dv1IntFlag = 0;
interObj.dv2IntFlag = 0;
interObj.dv3IntFlag = 0;
% initialize ints
interObj.srInd1 = [];
interObj.srInd2 = [];
interObj.srInd3 = [];
interObj.lrInd1 = [];
interObj.lrInd2 = [];
interObj.lrInd3 = [];
% Short range interactions
if isempty(particleObj.interHb)
  interObj.hardFlag = 0;
  interObj.hardId = 0;
  interObj.hard = 'none';
else
  interObj.hardFlag = 1;
  interObj.hard = particleObj.interHb;
  interObj.anyInter = 1;
  % rods
  if  strcmp( particleObj.type, 'rods' )
    interObj.typeId = 1;
    % flags/inds
    interObj.dv1IntFlag = 1;
    interObj.dv2IntFlag = 1;
    interObj.dv3IntFlag = 1;
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
    interObj.dv1IntFlag = 1;
    interObj.dv2IntFlag = 1;
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
    interObj.dv1IntFlag = 1;
    interObj.dv3IntFlag = 1;
    interObj.dv3IntFlag = 1;
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
      interObj.muExScale = (systemObj.tmp * ...
        systemObj.l3 * systemObj.l1 * systemObj.l2) ./ ...
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
else
  strInd = 1;
  corrInd = 2;
  paramInd = 3;
  % store some things
  numPotentials = length( particleObj.interactLrV );
  % Scale
  interObj.muMfScale = (systemObj.tmp * ...
    systemObj.l3 * systemObj.l1 * systemObj.l2) ./ ...
    (systemObj.n1 * systemObj.n2 * systemObj.n3);
  % add the potenials, c2
  c2 = 0;
  v = 0;
  for ii = 1:numPotentials
    foundV = 0;
    currPot = particleObj.interactLrV{ii};
    % soft shoulder 2D
    if strcmp( currPot{strInd}, 'ss' )
      foundV = 1;
      % build potential
      vTemp = SoftShoulderClass( currPot{strInd}, currPot{corrInd},...
        currPot{paramInd}(1), currPot{paramInd}(2), currPot{paramInd}(3),currPot{paramInd}(4),...
        systemObj.tmp, systemObj.n1, systemObj.l1, systemObj.n2, systemObj.l2 );
    end
    % polar align
    if strcmp( currPot{strInd}, 'pa' )
      foundV = 1;
      % build potential
      vTemp = PolarAlignClass( currPot{strInd}, currPot{corrInd},...
        currPot{paramInd}(1),  systemObj.tmp, systemObj.n3, systemObj.l3);
    end
    % polar align gaussian drop off
    if strcmp( currPot{strInd}, 'pag' )
      foundV = 1;
      % build potential
      vTemp = PolarAlignGaussClass( currPot{strInd}, currPot{corrInd},...
        currPot{paramInd}(1), currPot{paramInd}(2),...
        systemObj.tmp, systemObj.n1, systemObj.n2, systemObj.n3, ...
        systemObj.l1, systemObj.l2, systemObj.l3);
    end
    % decaying exponential 2D
    if strcmp( currPot{strInd}, 'de' )
      foundV = 1;
      % build potential
      vTemp = DecayExpClass( currPot{strInd}, currPot{corrInd},...
        currPot{paramInd}(1), currPot{paramInd}(2), ...
        systemObj.tmp, systemObj.n1, systemObj.n2, ...
        systemObj.l1, systemObj.l2 );
    end
    % gaussian
    if strcmp( currPot{strInd}, 'gauss' )
      foundV = 1;
      % build potential
      vTemp = PolarAlignGaussClass( currPot{strInd}, currPot{corrInd},...
        currPot{paramInd}(1), currPot{paramInd}(2),...
        systemObj.tmp, systemObj.n1, systemObj.n2, systemObj.n3, ...
        systemObj.l1, systemObj.l2, systemObj.l3);
    end
    if foundV
      interObj.anyInter = 1;
      interObj.longFlag = 1;
      fprintf('Long inter %s using %s correlations\n', ...
        currPot{strInd}, currPot{corrInd});
      interObj.interactLrV{ii} = vTemp;
      v = v + reshape( vTemp.V, vTemp.ReshapeInds );
      c2 = c2 + reshape( vTemp.C2, vTemp.ReshapeInds );
    end
  end % loop over potentials
  % get vFt and find shape
  [vn1, vn2, vn3] = size( c2 );
  if vn1 == 1
    interObj.lrInd1 = floor(systemObj.n1/2) + 1;
  else
    interObj.lrInd1 = 1:systemObj.n1;
    interObj.dv1IntFlag = 1;
  end
  if vn2 == 1
    interObj.lrInd2 = floor(systemObj.n2/2) + 1;
  else
    interObj.lrInd2 = 1:systemObj.n2;
    interObj.dv2IntFlag = 1;
  end
  if vn3 == 1
    interObj.lrInd3 = floor(systemObj.n3/2) + 1;
  else
    interObj.lrInd3 = 1:systemObj.n3;
    interObj.dv3IntFlag = 1;
  end
  % store it
  interObj.v = v;
  interObj.c2 = c2;
  interObj.c2Ft = fftshift( fftn( c2 ) );
end % long range interaction
% if there is long or short range interaction, set interaction flag
if interObj.dv1IntFlag; interObj.dv1Flag = 1; end
if interObj.dv2IntFlag; interObj.dv2Flag = 1; end
if interObj.dv3IntFlag; interObj.dv3Flag = 1; end
% External potentional
if isempty(particleObj.externalV)
  interObj.ext = 'none';
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
    fprintf('External %s\n', currPot{1} );
    addPot = 0;
    posVecs = cell(3,1);
    posVecs{1} = gridObj.x1;
    posVecs{2} = gridObj.x2;
    posVecs{3} = gridObj.x3;
    if strcmp( currPot{1}, 'linV' ) % linear
      % run a check just in case
      if nVec( currPot{2}(1) ) > 1
        fprintf('Along dim %d\n', currPot{2}(1) );
        vTemp = LinearVClass( currPot{1}, currPot{2}(1), ...
          currPot{2}(2), posVecs{ currPot{2}(1) } );
        interObj.externalV{ii} = vTemp;
        addPot = 1;
      end
    elseif strcmp( currPot{1}, 'quadV' )
      % run a dim check just in case
      if nVec( currPot{2}(1) ) > 1
        fprintf('Along dim %d\n', currPot{2}(1) );
        vTemp = QuadVClass( currPot{1}, currPot{2}(1), ...
          currPot{2}(2), posVecs{ currPot{2}(1) } );
        interObj.externalV{ii} = vTemp;
        addPot = 1;
      end
    elseif strcmp( currPot{1}, 'nemV' ) % nematic
      % run a dim check just in case
      if nVec( 3 ) > 1
        vTemp = NematicExternalVClass( currPot{1}, currPot{2}(1), ...
          currPot{2}(2), gridObj.x3 );
        interObj.externalV{ii} = vTemp;
        addPot = 1;
        %interObj.dv3Flag = 1;
      end
    elseif strcmp( currPot{1}, 'polV' ) % polar
      % run a dim check just in case
      if nVec( 3 ) > 1
        vTemp = PolarExternalVClass( currPot{1}, currPot{2}(1), ...
          currPot{2}(2), gridObj.x3 );
        interObj.externalV{ii} = vTemp;
        addPot = 1;
        %interObj.dv3Flag = 1;
      end
    else
      fprintf('Not building potential, potential not available\n' );
    end
    if addPot
      v = v + vTemp.VvReshape;
      dV.dx1 = dV.dx1 + vTemp.DvDx1;
      dV.dx2 = dV.dx2 + vTemp.DvDx2;
      dV.dx3 = dV.dx3 + vTemp.DvDx3;
    end
    if currPot{2}(1)  == 1
      interObj.dv1Flag = 1;
    end
    if currPot{2}(1)  == 2
      interObj.dv2Flag = 1;
    end
    if currPot{2}(1)  == 3
      interObj.dv3Flag = 1;
    end
  end % for loop
  interObj.dVExt = dV;
end % external
% print it
if interObj.hardFlag == 0; fprintf('No hard interactions\n'); end
if interObj.longFlag == 0; fprintf('No long interactions\n'); end
if interObj.extFlag == 0; fprintf('No external potentials\n'); end
% get ind
interObj.intInd1 = unique( [interObj.lrInd1, interObj.srInd1] );
interObj.intInd2 = unique( [interObj.lrInd2, interObj.srInd2] );
interObj.intInd3 = unique( [interObj.lrInd3, interObj.srInd3] );