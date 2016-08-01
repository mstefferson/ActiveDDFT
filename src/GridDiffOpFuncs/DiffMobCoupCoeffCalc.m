function [diffObj] = DiffMobCoupCoeffCalc(...
  T,Mob_par,Mob_perp,Mob_rot,...
  kx,ky,km,kx2D,ky2D,vd)
         
% Use Einstein diffusion relations
diffObj.Mob_par  = Mob_par;
diffObj.Mob_perp = Mob_perp;
diffObj.Mob_rot  = Mob_rot; 

diffObj.D_par  = Mob_par * T; % Parallel diffusion coeff
diffObj.D_perp = Mob_perp * T; % Perpendicular coeff
diffObj.D_rot  = Mob_rot * T; % Rotational diffusion

%Aniso diffusion coupling
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
  diffObj.CfMplus1  = -vd/2 * ( 1i .* kx2D - ky2D  ) ;
  diffObj.CfMminus1 = -vd/2 * ( 1i .* kx2D + ky2D  ) ;
end

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
jxf    = - ( diffObj.D_par + diffObj.D_perp ) .* ...
  ( sqrt(-1) * kx2D ) / 2;
if diffObj.Ani == 0
  jxMm2f = 0;
  jxMp2f = 0; 
else
  jxMm2f = - ( diffObj.D_par - diffObj.D_perp ) .* ...
    ( sqrt(-1) * kx2D + ky2D ) / 4;
  jxMp2f = - ( diffObj.D_par - diffObj.D_perp ) .* ...
    ( sqrt(-1) * kx2D - ky2D ) / 4;
end

% jy
jyf    = - ( diffObj.D_perp + diffObj.D_par ) .* ...
  ( sqrt(-1) * ky2D ) / 2;
if diffObj.Ani == 0
  jyMm2f = 0;
  jyMp2f = 0; 
else
  jyMm2f = - ( diffObj.D_perp - diffObj.D_par ) .* ...
    ( sqrt(-1) * ky2D - kx2D ) / 4;
  jyMp2f = - ( diffObj.D_perp - diffObj.D_par ) .* ...
    ( sqrt(-1) * ky2D + kx2D ) / 4;
end

% Rep
diffObj.jxf_reps = repmat( jxf, [1,1,Nm]);
diffObj.jyf_reps = repmat( jyf, [1,1,Nm]);
diffObj.jxMm2f_reps = repmat( jxMm2f, [1,1,Nm]);
diffObj.jxMp2f_reps = repmat( jxMp2f, [1,1,Nm]);
diffObj.jyMm2f_reps = repmat( jyMm2f, [1,1,Nm]);
diffObj.jyMp2f_reps = repmat( jyMp2f, [1,1,Nm]);

end
