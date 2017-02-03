function [interObj] = InterObjMaker( particleObj, systemObj)
% Store interactions in interobj
% Short range interactions
if isempty(particleObj.interHb)
  interObj.hardFlag = 0;
  interObj.hard = 'none';
  fprintf('No long range interactions\n');
else
  interObj.hardFlag = 1;
  interObj.hard = particleObj.interHbT;
  fprintf('Hard interactions %s\n', interObj.hard);
end
% Long range interactions
if isempty(particleObj.interLr)
  interObj.longFlag = 0;
  interObj.long = 'none';
  fprintf('No long range interactions\n');
else
  interObj.longFlag = 1;
  interObj.long = particleObj.interLrT;
  fprintf('Long interactions %s\n', interObj.long);
end
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
% Short range hard interactions
if  particleObj.type == 'rods' 
  if interObj.hard == 'mayer'
    interObj.hardSpec = 'rodsMayer';
    interObj.FmFt = fftshift(fftn( mayerFncHr(...
      systemObj.n1, systemObj.n2, systemObj.n3, systemObj.l1, ...
      systemObj.l2, particleObj.lMaj) ));
  else
    fprintf('Cannot find hard interactions\n')
  end
elseif particleObj.type == 'disks' 
  fprintf('Not written\n')
elseif particleObj.type == 'spheres' 
  fprintf('Not written\n')
else
  fprintf('Cannot find particle type\n');
  error('Cannot find particle type\n');
end
% Long range hard interactions not written
% External potential not written

