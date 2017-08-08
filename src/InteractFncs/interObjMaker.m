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
    % mayer
    if strcmp( interObj.hard, 'mayer' )
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
      fprintf('%s hard %s\n', interObj.hard, particleObj.type);
      interObj.srInd1 = 1:systemObj.n1;
      interObj.srInd2 = 1:systemObj.n2;
      interObj.srInd3 = 1:systemObj.n3;
    else
      fprintf('Cannot find hard rod interactions\n')
    end
  % disks
  elseif strcmp( particleObj.type, 'disks' ) 
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
if isempty(particleObj.interLr) || strcmp( particleObj.interLr{1}, '' )
  interObj.longFlag = 0;
  interObj.longId = 0;
  interObj.long = 'none';
  fprintf('No long range interactions\n');
else
  % store some things
  interObj.longFlag = 1;
  interObj.long = particleObj.interLr;
  numPotentials = length( interObj.long );
  interObj.anyInter = 1;
  fprintf('Long interactions\n');
  % Scale
  interObj.muMfScale = (systemObj.l3 * systemObj.l1 * systemObj.l2) ./ ...
    (systemObj.n1 * systemObj.n2 * systemObj.n3);
  % factors
  interObj.lrLs1 = particleObj.lrLs1;
  interObj.lrLs2 = particleObj.lrLs2;
  interObj.lrEs1 = particleObj.lrEs1;
  interObj.lrEs2 = particleObj.lrEs1;
  % add the potenials, v
  v = 0;
  for ii = 1:numPotentials
    % soft shoulder 2D
    if strcmp( interObj.long{ii}, 'ss2d' )
      fprintf('Long interactions %s, soft shoulder 2d\n', interObj.long{ii});
      interObj.longId(ii) = 1;
      % build potential
      [vTemp] = softShoulder2d( interObj.lrEs1(ii), interObj.lrEs2(ii), ...
        interObj.lrLs1(ii), interObj.lrLs2(ii), ...
        systemObj.n1, systemObj.l1, systemObj.n2, systemObj.l2 );
      v = v + reshape( vTemp, [systemObj.n1 systemObj.n2 1] );
    end
    % polar align 2D
    if strcmp( interObj.long{ii}, 'pa2d' ) 
      fprintf('Long interactions %s, polar align 2d\n', interObj.long{ii});
      interObj.longId(ii) = 2;
      % build potential
      [vTemp] = polarAlign2d( interObj.lrEs1(ii),  systemObj.n3, systemObj.l3 );
      v = v + reshape( vTemp, [1, 1, systemObj.n3] );
    end
    % polar alig gauss 2D
    if strcmp( interObj.long{ii}, 'pag2d' )
      fprintf('Long interactions %s, polar align gauss 2d\n', interObj.long{ii});
      interObj.longId(ii) = 3;
      % build potential
      [vTemp] = polarAlignGaussian2d( interObj.lrEs1(ii), interObj.lrLs1(ii), systemObj.n3, systemObj.l3 );
      v = v + vTemp;
    end
    % decaying exponential 2D
    if strcmp( interObj.long{ii}, 'de2d' )
      fprintf('Long interactions %s, decaying exponential 2d\n', interObj.long{ii});
      interObj.longId(ii) = 4;
      % build potential
      [vTemp] = decayexp2d( interObj.lrEs1(ii), interObj.lrLs1(ii), ...
      systemObj.n1, systemObj.l1, systemObj.n2, systemObj.l2);
      v = v + vTemp;
    end
  end % loop over potentials
  % get vFt and find shape
  [vn1, vn2, vn3] = size( v );
  if vn1 == 1
    interObj.lrInd1 = floor(systemObj.n1/2) + 1;
  else
    interObj.lrInd1 = 1:systemObj.n1;
  end
  if vn2 == 1
    interObj.lrInd2 = floor(systemObj.n2/2) + 1;
  else
    interObj.lrInd2 = 1:systemObj.n2;
  end
  if vn3 == 1
    interObj.lrInd3 = floor(systemObj.n3/2) + 1;
  else
    interObj.lrInd3 = 1:systemObj.n3;
  end
  % store it
  interObj.v = v;
  interObj.vFt = fftshift( fftn( v ) );
end % long range interaction
% External potentional
if isempty(particleObj.externalPot)
  interObj.extFlag = 0;
  interObj.ext = 'none';
  fprintf('No external potential\n');
else
  interObj.extFlag = 1;
  interObj.ext = particleObj.externalPot;
  fprintf('External %s\n', interObj.ext);
  if strcmp( interObj.ext, 'linV1' )
    x = reshape( gridObj.x1, [ systemObj.n1, 1 ] );
    [interObj.v, interObj.vFt] = linearV( particleObj.exEs1, x );
    dV.f1 = dVCalc(interObj.vFt, systemObj, diffObj, interObj);
  end
  if strcmp( interObj.ext, 'linV1' )
    x = reshape( gridObj.x1, [ systemObj.n1, 1 ] );
    [interObj.v, interObj.vFt, dv] = linearV( particleObj.exEs1, x );
    dV.dx1 =  dv;
    dV.dx2 =  0;
    dV.dx3 =  0;
  end
  if strcmp( interObj.ext, 'quadV1' )
    x = reshape( gridObj.x1, [ systemObj.n1, 1 ] );
    [interObj.v, interObj.vFt, dv] = linearV( particleObj.exEs1, x );
    dV.dx1 =  dv;
    dV.dx2 =  0;
    dV.dx3 =  0;
  end
  interObj.dV = dV;
  keyboard
end
