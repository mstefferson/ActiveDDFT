function [diffObj] = DiffMobCoupCoeffCalc(...
  T,Mob_par,Mob_perp,Mob_rot,dt,delta_x,delta_phi,...
  kx,ky,km,kx2D,ky2D,vd)
         
% Use Einstein diffusion relations
diffObj.Mob_par  = Mob_par;
diffObj.Mob_perp = Mob_perp;
diffObj.Mob_rot  = Mob_rot; 

diffObj.D_par  = Mob_par * T; % Parallel diffusion coeff
diffObj.D_perp = Mob_perp * T; % Perpendicular coeff
diffObj.D_rot  = Mob_rot * T; % Rotational diffusion

% Check stability condition
StabCoeffPar  = diffObj.D_par  .* dt / (delta_x ^2 );
StabCoeffPerp = diffObj.D_perp .* dt / (delta_x ^2 );
StabCoeffRot  = diffObj.D_rot  .* dt / (delta_phi ^2 );

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
CrossTermFactor = (diffObj.D_par - diffObj.D_perp)/4; 
% Coupling coefficent
diffObj.CfMplus2 = CrossTermFactor.*(ky2D - 1i.*kx2D).^2; 
diffObj.CfMminus2  = CrossTermFactor.*(ky2D + 1i.*kx2D).^2; 

% Driven part
diffObj.CfMplus1  = -vd/2 * ( 1i .* kx2D - ky2D  ) ;
diffObj.CfMminus1 = -vd/2 * ( 1i .* kx2D + ky2D  ) ;

% Multiply by derivatives in k space repeatedily,
% they are large, but shouild be worth allocating
[ diffObj.iky3, diffObj.ikx3, diffObj.ikm3 ] = ...
  meshgrid(ky,kx,km);

diffObj.ikx3 = sqrt(-1) .* diffObj.ikx3;
diffObj.iky3 = sqrt(-1) .* diffObj.iky3;
diffObj.ikm3 = sqrt(-1) .* diffObj.ikm3;

% NL couplings
% key: jxMm2f = jx m -2 coupling factor

% jx
Nm = length(km);
diffObj.jxf    = - ( diffObj.D_par + diffObj.D_perp ) .* ...
  ( sqrt(-1) * kx2D ) / 2;
diffObj.jxMm2f = - ( diffObj.D_par - diffObj.D_perp ) .* ...
  ( sqrt(-1) * kx2D + ky2D ) / 4;
diffObj.jxMp2f = - ( diffObj.D_par - diffObj.D_perp ) .* ...
  ( sqrt(-1) * kx2D - ky2D ) / 4;

% jy
diffObj.jyf    = - ( diffObj.D_perp + diffObj.D_par ) .* ...
  ( sqrt(-1) * ky2D ) / 2;
diffObj.jyMm2f = - ( diffObj.D_perp - diffObj.D_par ) .* ...
  ( sqrt(-1) * ky2D - kx2D ) / 4;
diffObj.jyMp2f = - ( diffObj.D_perp - diffObj.D_par ) .* ...
  ( sqrt(-1) * ky2D + kx2D ) / 4;

% Reps
diffObj.jxf_reps = repmat( diffObj.jxf, [1,1,Nm]);
diffObj.jyf_reps = repmat( diffObj.jyf, [1,1,Nm]);
diffObj.jxMm2f_reps = repmat( diffObj.jxMm2f, [1,1,Nm]);
diffObj.jxMp2f_reps = repmat( diffObj.jxMp2f, [1,1,Nm]);
diffObj.jyMm2f_reps = repmat( diffObj.jyMm2f, [1,1,Nm]);
diffObj.jyMp2f_reps = repmat( diffObj.jyMp2f, [1,1,Nm]);

end
