% function DiffMobCoupCoeffCalc
% 
% Description: Builds diffusion object which includes vectors that will go
% into operator
%
function [diffObj] = DiffMobCoupCoeffCalc(...
  T,Mob,Mob_par,Mob_perp,Mob_rot,...
  kx,ky,km,kx2D,ky2D,vd) 
% Use Einstein diffusion relations
diffObj.Mob_pos  = Mob;
diffObj.Mob_par  = Mob_par;
diffObj.Mob_perp = Mob_perp;
diffObj.Mob_rot  = Mob_rot; 
diffObj.D_pos = Mob_perp * T; % Perpendicular coeff
diffObj.D_par  = Mob_par * T; % Parallel diffusion coeff
diffObj.D_perp = Mob_perp * T; % Perpendicular coeff
diffObj.D_rot  = Mob_rot * T; % Rotational diffusion

%% DEBUG %%
diffObj.D_par = 0;
diffObj.D_perp = 0;
diffObj.D_pos = 0;

%Aniso rods diffusion coupling
% Constant in front of cross terms
CrossTermFactor = (diffObj.D_par - diffObj.D_perp)/4; 
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
% Driven part
if vd == 0
  diffObj.dr = 0;
  diffObj.CfMplus1  = 0;
  diffObj.CfMminus1 = 0;
else
  diffObj.dr = 1;
  diffObj.CfMplus1  = -vd/2 * ( 1i .* kx2D - ky2D  ) ;
  diffObj.CfMminus1 = -vd/2 * ( 1i .* kx2D + ky2D  ) ;
end
% Multiply by derivatives in k space repeatedily,
% they are large, but shouild be worth allocating
[ diffObj.ik2rep3, diffObj.ik1rep3, diffObj.ik3rep3 ] = ...
  meshgrid(ky,kx,km);
% 3d k vectors
diffObj.ik1rep3 = sqrt(-1) .* diffObj.ik1rep3;
diffObj.ik2rep3 = sqrt(-1) .* diffObj.ik2rep3;
diffObj.ik3rep3 = sqrt(-1) .* diffObj.ik3rep3;
% NL couplings
% key: jiMm2f = ji m -2 coupling factor
% jx
n3 = length(km);
j1f    = - ( diffObj.D_par + diffObj.D_perp ) .* ...
  ( sqrt(-1) * kx2D ) / 2;
diffObj.j1f_reps = repmat( j1f, [1,1,n3]);
if diffObj.Ani == 0 || diffObj.D_perp == diffObj.D_par
  diffObj.j1Mm2f_reps = 0;
  diffObj.j1Mp2f_reps = 0;
else
  jxMm2f = - ( diffObj.D_par - diffObj.D_perp ) .* ...
    ( sqrt(-1) * kx2D + ky2D ) / 4;
  jxMp2f = - ( diffObj.D_par - diffObj.D_perp ) .* ...
    ( sqrt(-1) * kx2D - ky2D ) / 4;
  diffObj.j1Mm2f_reps = repmat( jxMm2f, [1,1,n3]);
  diffObj.j1Mp2f_reps = repmat( jxMp2f, [1,1,n3]);
end
% jy
j2f    = - ( diffObj.D_perp + diffObj.D_par ) .* ...
  ( sqrt(-1) * ky2D ) / 2;
diffObj.j2f_reps = repmat( j2f, [1,1,n3]);
if diffObj.Ani == 0 || diffObj.D_perp == diffObj.D_par
  diffObj.j2Mm2f_reps = 0;
  diffObj.j2Mp2f_reps = 0;
else
  jyMm2f = - ( diffObj.D_perp - diffObj.D_par ) .* ...
    ( sqrt(-1) * ky2D - kx2D ) / 4;
  jyMp2f = - ( diffObj.D_perp - diffObj.D_par ) .* ...
    ( sqrt(-1) * ky2D + kx2D ) / 4;
  diffObj.j2Mm2f_reps = repmat( jyMm2f, [1,1,n3]);
  diffObj.j2Mp2f_reps = repmat( jyMp2f, [1,1,n3]);
end
% jm
diffObj.jx3f = - diffObj.ik3rep3 .* diffObj.D_rot;

end
