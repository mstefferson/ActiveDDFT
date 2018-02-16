% function DiffMobCoupCoeffCalc
% 
% Description: Builds diffusion object which includes vectors that will go
% into operator
%
function [diffObj] = DiffMobCoupCoeffCalc(...
  T,Mob,Mob_par,Mob_perp,Mob_rot,...
  kx,ky,km,kx2D,ky2D,vd, phi) 
% Use Einstein diffusion relations
diffObj.Mob_pos  = Mob;
diffObj.Mob_par  = Mob_par;
diffObj.Mob_perp = Mob_perp;
diffObj.Mob_rot  = Mob_rot; 
diffObj.D_pos = Mob_perp * T; % Perpendicular coeff
diffObj.D_par  = Mob_par * T; % Parallel diffusion coeff
diffObj.D_perp = Mob_perp * T; % Perpendicular coeff
diffObj.D_rot  = Mob_rot * T; % Rotational diffusion
% Build the actual mobility matrix
phi = reshape( phi, [1 1 length(phi) ] );
cosphi = cos(phi);
sinphi = sin(phi);
if Mob_par == Mob_perp % iso
  diffObj.Ani == 0 
  diffObj.MobMat11 = diffObj.Mob_par;
  diffObj.MobMat12 = 0;
  diffObj.MobMat22 = diffObj.Mob_par;
else % aniso
  diffObj.Ani == 1 
  diffObj.MobMat11 = diffObj.Mob_par .* cosphi .^ 2 
    + diffObj.Mob_perp  .* sinphi .^ 2;
  diffObj.MobMat12 = (diffObj.Mob_par - diffObj.Mob_perp) / 2 * cosphi .* sinphi;
  diffObj.MobMat22 = diffObj.Mob_par .* sinphi .^ 2 
    + diffObj.Mob_perp  .* cosphi .^ 2;
end
diffObj.MobMat33 = diffObj.Mob_rot; 
% Aniso rods diffusion coupling
% Constant in front of cross terms
CrossTermFactor = (diffObj.Mob_par - diffObj.Mob_perp)/4; 
% Coupling coefficent
if CrossTermFactor == 0
  diffObj.Ani = 0;
  diffObj.CfMplus2 = 0; 
  diffObj.CfMminus2 = 0;
else
  diffObj.Ani = 1;
  diffObj.CfMplus2 = CrossTermFactor.*(ky2D - 1i.*kx2D).^2; 
  diffObj.CfMminus2 = CrossTermFactor.*(ky2D + 1i.*kx2D).^2; 
end
% k vector for gradient
diffObj.ik1 = sqrt(-1) .* reshape( kx, [ systemObj.n1 1 ] );
diffObj.ik2 = sqrt(-1) .* reshape( ky, [ 1 systemObj.n2 ] );
diffObj.ik3 = sqrt(-1) .* reshape( km, [ 1 1 systemObj.n3 ] );
end
