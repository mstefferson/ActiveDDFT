% Calculates the mobility and diffusion tensor based on dimensions of rod
% and fluid

function [diffObj] = DiffMobCoupCoeffCalcIsoDiff(...
  T,Mob_pos,Mob_rot, kx, ky,km)

% Use Einstein diffusion relations
diffObj.Mob_pos = Mob_pos; 
diffObj.Mob_rot = Mob_rot; 

diffObj.D_pos = Mob_pos * T; 
diffObj.D_rot = Mob_rot * T; 

% Multiply by derivatives in k space repeatedily,
% they are large, but shouild be worth allocating
[ diffObj.iky3, diffObj.ikx3, diffObj.ikm3 ] = ...
  meshgrid(ky,kx,km);

diffObj.ikx3 = sqrt(-1) .* diffObj.ikx3;
diffObj.iky3 = sqrt(-1) .* diffObj.iky3;
diffObj.ikm3 = sqrt(-1) .* diffObj.ikm3;

end
