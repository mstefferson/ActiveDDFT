% Calculates the mobility and diffusion tensor based on dimensions of rod
% and fluid

function [DiffMobObj] = DiffMobCoupCoeffCalcIsoDiff(...
  T,Mob_pos,Mob_rot, kx, ky,km)

% Use Einstein diffusion relations
DiffMobObj.Mob_pos = Mob_pos; 
DiffMobObj.Mob_rot = Mob_rot; 

DiffMobObj.D_pos = Mob_pos * T; 
DiffMobObj.D_rot = Mob_rot * T; 

% Multiply by derivatives in k space repeatedily,
% they are large, but shouild be worth allocating
[ DiffMobObj.iky3, DiffMobObj.ikx3, DiffMobObj.ikm3 ] = ...
  meshgrid(ky,kx,km);

DiffMobObj.ikx3 = sqrt(-1) .* DiffMobObj.ikx3;
DiffMobObj.iky3 = sqrt(-1) .* DiffMobObj.iky3;
DiffMobObj.ikm3 = sqrt(-1) .* DiffMobObj.ikm3;

end
