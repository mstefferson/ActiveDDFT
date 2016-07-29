% Calculates the mobility and diffusion tensor based on dimensions of rod
% and fluid

function [diffObj]...
             = DiffMobCoupCoeffCalcOld(wfid,Eta_visc,L_rod,Diam,T,dt,delta_x,delta_phi,kx2D,ky2D)
%Coefficients of the friction matrix
Fric_par  = pi * Eta_visc * L_rod / (log(L_rod/Diam));          % for motion parallel to rods
Fric_perp = 2*Fric_par;                                         % for motion perpendicular to rods
Fric_rot  = pi * Eta_visc * L_rod ^3 / (3 * log(L_rod/Diam));   % for rotational motion

Mob_par  = Fric_par^(-1);                                           % Parallel mobility coefficient
Mob_perp = Fric_perp^(-1);                                          % Perpendicular mobility coefficient
Mob_rot  = Fric_rot^(-1);                                           % Rotational mobility coefficient

% Use Einstein diffusion relations
D_par  = Mob_par * T;                                            % Parallel diffusion coeff
D_perp = Mob_perp * T;                                           % Perpendicular coeff
D_rot  = Mob_rot * T;                                            % Rotational diffusion

% Check stability condition
StabCoeffPar  = D_par  .* dt / (delta_x ^2 );
StabCoeffPerp = D_perp .* dt / (delta_x ^2 );
StabCoeffRot  = D_rot  .* dt / (delta_phi ^2 );

if StabCoeffPar > 1/2
    fprintf(wfid,'My lord! I foresee bad times ahead. Your gridspacing might not lead to convergence. The parallel diffusion coefficien is to blame\n');
    fprintf(wfid,'StabCoeffPar = %f (should be less than 1/2) \n', StabCoeffPar);
end
if StabCoeffPerp > 1/2
    fprintf(wfid,'My lord! I foresee bad times ahead. Your gridspacing might not lead to convergence. The perpendicular diffusion coefficien is to blame\n');
    fprintf('StabCoeffPerp = %f (should be less than 1/2) \n', StabCoeffPerp);
end
if StabCoeffRot > 1/2
    fprintf(wfid,'My lord! I foresee bad times ahead. Your gridspacing might not lead to convergence. The rotational diffusion coefficien is to blame\n');
    fprintf(wfid,'StabCoeffRot = %f (should be less than 1/2) \n', StabCoeffRot);
end

CrossTermFactor = (D_par - D_perp)/4;                       % Constant in front of cross terms
CoupFacMplus2   = CrossTermFactor.*(ky2D - 1i.*kx2D).^2;      % Coupling coefficent
CoupFacMminus2  = CrossTermFactor.*(ky2D + 1i.*kx2D).^2;     % Coupling coefficent

diffObj = struct('Mob_par', Mob_par, 'Mob_perp', Mob_perp, 'Mob_rot',Mob_rot,...
    'D_par',D_par, 'D_perp',D_perp, 'D_rot',D_rot, ...
    'CfMplus2',CoupFacMplus2,'CfMminus2',CoupFacMminus2);