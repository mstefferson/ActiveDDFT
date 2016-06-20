function [DiffMobObj] = DiffMobCoupCoeffCalc(...
  T,Mob_par,Mob_perp,Mob_rot,dt,delta_x,delta_phi,...
  kx,ky,km,kx2D,ky2D,vd)
         
% Use Einstein diffusion relations
DiffMobObj.Mob_par  = Mob_par;
DiffMobObj.Mob_perp = Mob_perp;
DiffMobObj.Mob_rot  = Mob_rot; 

DiffMobObj.D_par  = Mob_par * T; % Parallel diffusion coeff
DiffMobObj.D_perp = Mob_perp * T; % Perpendicular coeff
DiffMobObj.D_rot  = Mob_rot * T; % Rotational diffusion

% Check stability condition
StabCoeffPar  = DiffMobObj.D_par  .* dt / (delta_x ^2 );
StabCoeffPerp = DiffMobObj.D_perp .* dt / (delta_x ^2 );
StabCoeffRot  = DiffMobObj.D_rot  .* dt / (delta_phi ^2 );

if StabCoeffPar > 1/2
    fprintf('StabCoeffPar = %f (should be less than 1/2) \n', StabCoeffPar);
end
if StabCoeffPerp > 1/2
    fprintf('StabCoeffPerp = %f (should be less than 1/2) \n', StabCoeffPerp);
end
if StabCoeffRot > 1/2
    fprintf('StabCoeffRot = %f (should be less than 1/2) \n', StabCoeffRot);
end

%Aniso diffusion coupling
% Constant in front of cross terms
CrossTermFactor = - (DiffMobObj.D_par - DiffMobObj.D_perp) / 4; 
% Coupling coefficent
DiffMobObj.CfMminus2  = CrossTermFactor .* ( kx2D - sqrt(-1) .* ky2D ) .^ 2; 
DiffMobObj.CfMplus2   = CrossTermFactor .* ( kx2D + sqrt(-1) .* ky2D ) .^ 2; 

% Driven part
DiffMobObj.CfMminus1 = -vd/2 * ( sqrt(-1) .* kx2D + ky2D  ) ;
DiffMobObj.CfMplus1  = -vd/2 * ( sqrt(-1) .* kx2D - ky2D  ) ;

% Multiply by derivatives in k space repeatedily,
% they are large, but shouild be worth allocating
[ DiffMobObj.iky3, DiffMobObj.ikx3, DiffMobObj.ikm3 ] = ...
  meshgrid(ky,kx,km);

DiffMobObj.ikx3 = sqrt(-1) .* DiffMobObj.ikx3;
DiffMobObj.iky3 = sqrt(-1) .* DiffMobObj.iky3;
DiffMobObj.ikm3 = sqrt(-1) .* DiffMobObj.ikm3;

% NL couplings
% key: jxMm2f = jx m -2 coupling factor

% jx
Nm = length(km);
DiffMobObj.jxf    = - ( DiffMobObj.D_par + DiffMobObj.D_perp) / 2 * ...
  ( sqrt(-1) * kx2D );
DiffMobObj.jxMm2f = - ( DiffMobObj.D_par - DiffMobObj.D_perp) / 4 * ...
  ( sqrt(-1) * kx2D + ky2D );
DiffMobObj.jxMp2f = - ( DiffMobObj.D_par - DiffMobObj.D_perp) / 4 * ...
  ( sqrt(-1) * kx2D - ky2D );

% jy
DiffMobObj.jyf    = - ( DiffMobObj.D_par + DiffMobObj.D_perp) / 2 * ...
  ( sqrt(-1) * ky2D );
DiffMobObj.jyMm2f = - ( DiffMobObj.D_perp - DiffMobObj.D_par) / 4 * ...
  ( sqrt(-1) * ky2D - kx2D );
DiffMobObj.jyMp2f = - ( DiffMobObj.D_perp - DiffMobObj.D_par) / 4 * ...
  ( sqrt(-1) * ky2D + kx2D );

% Reps
DiffMobObj.jxf_reps = repmat( DiffMobObj.jxf, [1,1,Nm]);
DiffMobObj.jyf_reps = repmat( DiffMobObj.jyf, [1,1,Nm]);
DiffMobObj.jxMm2f_reps = repmat( DiffMobObj.jxMm2f, [1,1,Nm]);
DiffMobObj.jxMp2f_reps = repmat( DiffMobObj.jxMp2f, [1,1,Nm]);
DiffMobObj.jyMm2f_reps = repmat( DiffMobObj.jyMm2f, [1,1,Nm]);
DiffMobObj.jyMp2f_reps = repmat( DiffMobObj.jyMp2f, [1,1,Nm]);

end
