function [DiffMobObj]...
             = DiffMobCoupCoeffCalc(wfid,T,Mob_par,Mob_perp,Mob_rot,delta_t,delta_x,delta_phi,kx2D,ky2D,vd)
         
% Use Einstein diffusion relations
D_par  = Mob_par * T;                                            % Parallel diffusion coeff
D_perp = Mob_perp * T;                                           % Perpendicular coeff
D_rot  = Mob_rot * T;                                            % Rotational diffusion

% Check stability condition
StabCoeffPar  = D_par  .* delta_t / (delta_x ^2 );
StabCoeffPerp = D_perp .* delta_t / (delta_x ^2 );
StabCoeffRot  = D_rot  .* delta_t / (delta_phi ^2 );

if StabCoeffPar > 1/2
    fprintf(wfid,'My lord! I foresee bad times ahead. Your gridspacing might not lead to convergence. The parallel diffusion coefficien is to blame\n');
    fprintf(wfid,'StabCoeffPar = %f (should be less than 1/2) \n', StabCoeffPar);
end
if StabCoeffPerp > 1/2
    fprintf(wfid,'My lord! I foresee bad times ahead. Your gridspacing might not lead to convergence. The perpendicular diffusion coefficien is to blame\n');
    fprintf('StabCoeffPerp = %f (should be less than 1/2) \n', StabCoeffPerp);
end
if StabCoeffRot > 1/2
    fprintf(wfid,'My lord! I foresee bad times ahead. Your gridspacing might not lead to convergence. The rotational diffusion coefficien is to blame\n')
    fprintf(wfid,'StabCoeffRot = %f (should be less than 1/2) \n', StabCoeffRot);
end

%Aniso diffusion coupling
CrossTermFactor = (D_par - D_perp)/4;                       % Constant in front of cross terms
CoupFacMplus2   = CrossTermFactor.*(ky2D - 1i.*kx2D).^2;      % Coupling coefficent
CoupFacMminus2  = CrossTermFactor.*(ky2D + 1i.*kx2D).^2;     % Coupling coefficent

% Driven part
CoupFacMplus1  = -vd/2 * ( 1i .* kx2D - ky2D  ) ;
CoupFacMminus1 = -vd/2 * ( 1i .* kx2D + ky2D  ) ;

% keyboard

DiffMobObj = struct('Mob_par', Mob_par, 'Mob_perp', Mob_perp, 'Mob_rot',Mob_rot,...
    'D_par',D_par, 'D_perp',D_perp, 'D_rot',D_rot, ...
    'CfMplus2',CoupFacMplus2,'CfMminus2',CoupFacMminus2, ...
    'CfMplus1',CoupFacMplus1,'CfMminus1',CoupFacMminus1);
