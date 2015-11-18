
% Function returns the spatial concentration and Order Parameters of a 2-D system of particles
% with an orientation

function [C,POP,nx_POP,ny_POP,NOP,NADx,NADy] = OpCPNCalc(Nx,Ny,rho,phi,x,y,phi3d)

% Director is chosen to be in the +y direction for all gridpoints

%Output key
% C: Concentration
% POP: Scalar polar order parameter
% nx_POP: polar order parameter in the x direction. 1st moment of the
% distribution w.r.t. orientation.
% ny_POP:polar order parameter in the y direction. 1st moment of the
% distribution w.r.t. orientation.
% QeigMax_MB: Eigenvalue of the matrix nn - I/2
% NOP: The nematic order parameter. Defined as 3/2*  max(eigenvalue( nematic
% order parameter S_NOP)
% NADx: Nematic alignment director in x direction. Nematic alignment
% director is the eigenvector corresponding to max(eig(S_NOP))
% NADy: Nematic alignment director in y direction. Nematic alignment
% director is the eigenvector corresponding to max(eig(S_NOP))

%Concentration is the first moment of the distribution. Integrate over all
%angles
if Nx == 1 && Ny == 1
    C = trapz_periodic(phi,rho);
else
    C = trapz_periodic(phi,rho,3);
end
% C = trapz(phi,rho,3);
% keyboard
% Calculate the first moment of the distribution orientation. This gives
% the orientation field
%
if Nx == 1 && Ny == 1
    nx_POP = trapz_periodic(phi, cos(phi).*rho) ./ C; %Polar order parameter in x-direction
    ny_POP = trapz_periodic(phi, sin(phi).*rho) ./ C; %Polar order parameter in y-direction
else
    nx_POP = trapz_periodic(phi, cos(phi3d).*rho,3) ./ C; %Polar order parameter in x-direction
    ny_POP = trapz_periodic(phi, sin(phi3d).*rho,3) ./ C; %Polar order parameter in y-direction
end
% nx_POP = trapz(phi, cos(phi3d).*rho,3) ./ C; %Polar order parameter in x-direction
% ny_POP = trapz(phi, sin(phi3d).*rho,3) ./ C; %Polar order parameter in y-direction


POP = sqrt(nx_POP.^2 + ny_POP.^2);

%%%%%%%%%%%%%%Q matrix%%%%%%%%%%%%%%%%%
% Nematic Order parameter Q.

% keyboard

eigMaxQ_NOP = zeros(Nx,Ny);    % Eigenvalue of nemativ order parameter matrix
NADx = zeros(Nx,Ny);           % Nematic alignment x-direction
NADy = zeros(Nx,Ny);           % Nematic alignment y-direction
if Nx == 1 && Ny == 1
    Q_NOPxx_temp = trapz_periodic(phi,rho .* (cos(phi).*cos(phi) - 1/2)) ./ C;
    Q_NOPxy_temp = trapz_periodic(phi,rho .* (cos(phi).*sin(phi))) ./ C;
    Q_NOPyy_temp = trapz_periodic(phi,rho .* (sin(phi).*sin(phi) - 1/2)) ./ C;
else
    Q_NOPxx_temp = trapz_periodic(phi,rho .* (cos(phi3d).*cos(phi3d) - 1/2),3) ./ C;
    Q_NOPxy_temp = trapz_periodic(phi,rho .* (cos(phi3d).*sin(phi3d)),3) ./ C;
    Q_NOPyy_temp = trapz_periodic(phi,rho .* (sin(phi3d).*sin(phi3d) - 1/2),3) ./ C;
end


for i = 1:1:length(x)
    for j = 1:length(y)
        Q_temp = [Q_NOPxx_temp(i,j) Q_NOPxy_temp(i,j); Q_NOPxy_temp(i,j) Q_NOPyy_temp(i,j)];
        [EigVec,EigS] = eigs(Q_temp);
        eigMaxQ_NOP(i,j) = max(max(EigS));
        %Find the eigenvector corresponding to this eigenvalue
        NADtemp = EigVec( EigS == repmat(max(EigS,[],2),[1,2]) );
        NADx(i,j) = NADtemp(1);
        NADy(i,j) = NADtemp(2);
    end
end
% We need to build 2x2 matrices and diagonalize them. But, we know how to
% easily diagonalize a 2x2 matrix so can we find all the eigenvalues in 1
% step.
%  eigMaxQ_NOP = sqrt( ((Q_NOPxx_temp - Q_NOPyy_temp )./2).^2 + Q_NOPxy_temp.^2);
% keyboard
%Calculate Nematic Scalar parameter
NOP = 2*eigMaxQ_NOP;

%  keyboard
end
