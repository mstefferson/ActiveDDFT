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
      [~,interObj.FmFt] = mayerFncHr(...
        systemObj.n1, systemObj.n2, systemObj.n3, ...
        systemObj.l1, systemObj.l2, particleObj.lMaj) ;
      interObj.muMayerScale = - ( systemObj.tmp * systemObj.l3 * systemObj.l1 * systemObj.l2) ./ ...
        (systemObj.n1 * systemObj.n2 * systemObj.n3);
      fprintf('%s hard %s\n', interObj.hard, particleObj.type);
    else
      fprintf('Cannot find hard rod interactions\n')
    end
  % disks
  elseif strcmp( particleObj.type, 'disks' ) 
    interObj.typeId = 2;
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
if isempty(particleObj.interLr)
  interObj.longFlag = 0;
  interObj.longId = 0;
  interObj.long = 'none';
  fprintf('No long range interactions\n');
else 
  interObj.longFlag = 1;
  interObj.long = particleObj.interLr;
  interObj.anyInter = 1;
  fprintf('Long interactions %s\n', interObj.long);
  % Scale
  interObj.muMfScale = (systemObj.l3 * systemObj.l1 * systemObj.l2) ./ ...
    (systemObj.n1 * systemObj.n2 * systemObj.n3);
  % rods
  if  strcmp( particleObj.type, 'rods' ) 
    interObj.longFlag = 0;
    interObj.long = 'none';
    fprintf('No long range interactions written for rods\n')
  % disks
  elseif  strcmp( particleObj.type, 'disks' ) 
    if strcmp( interObj.long, 'softshoulder' )
      interObj.longSpec = 'disksSoft';
      interObj.longId = 1;
      % scales
      interObj.lrLs1 = particleObj.lrLs1;
      interObj.lrLs2 = particleObj.lrLs2;
      interObj.lrEs1 = particleObj.lrEs1;
      interObj.lrEs2 = particleObj.lrEs1;
      % build soft should potential
      [~, interObj.vFt] = softShoulder2D( interObj.lrEs1, interObj.lrEs2, ...
        interObj.lrLs1, interObj.lrLs2, ...
        systemObj.n1, systemObj.l1, systemObj.n2, systemObj.l2 );
    end
  % spheres
  elseif strcmp( particleObj.type, 'spheres' )
    interObj.longFlag = 0;
    interObj.long = 'none';
    fprintf('No long range interactions written for spheres\n')
    fprintf('No long range interactions\n');
  else
    fprintf('Cannot find particle type\n');
    error('Cannot find particle type\n');
  end % particle type
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
end
