function particleMaster = particleMobInit( particleMaster, diagFlag )
% Rods
if strcmp( particleMaster.type, 'rods' )
  particleMaster.lMaj = 1;  % Length along the major axis
  particleMaster.lMin = 0;  % Length along the minor axis
  particleMaster.b  = particleMaster.lMaj^2/pi; % Average excluded volume per particle
if diagFlag 
    particleMaster.mobPar = particleMaster.mob; 
  else
    particleMaster.mobPar  = 2*particleMaster.mob; 
  end
  particleMaster.mobPerp   = particleMaster.mob;
  particleMaster.mobRot   = 6 * particleMaster.mob / particleMaster.lMaj^2;
% Disks 
elseif strcmp( particleMaster.type, 'disks' )
  particleMaster.lMaj = 1;  % Length along the major axis
  particleMaster.lMin = 1;  % Length along the minor axis
  particleMaster.b  = pi .* particleMaster.lMaj^2 / 4; % packing fraction
  particleMaster.mobPar = particleMaster.mob;
  particleMaster.mobPerp = particleMaster.mob;
  particleMaster.mobRot = particleMaster.mob / ( 3 * particleMaster.lMaj^2 );
% Spheres 
elseif strcmp( particleMaster.type, 'spheres' )
  fprintf('Not written\n')
else
  fprintf('Cannot find particle type\n');
  error('Cannot find particle type\n');
end

