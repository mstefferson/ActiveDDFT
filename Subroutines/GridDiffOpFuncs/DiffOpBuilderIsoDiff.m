% Build the diffusion operator
% Currently uncalled

function [Lop] = DiffOpBuilderIsoDiff(DiffMobObj,GridObj,Nx,Ny,Nm,N3)

% Build a strange km for repmat
km = zeros( 1, 1, Nm );
km(1,1,:) = GridObj.km;
%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( ( (DiffMobObj.D_par + DiffMobObj.D_perp)./ 2 ) .* ...
  ( repmat( GridObj.kx', [1, Ny, Nm] ) .^ 2 + ...
  repmat( GridObj.ky, [Nx, 1, Nm ]) .^ 2 ) + ...
  DiffMobObj.D_rot .* repmat( km, [Nx, Ny, 1] ) .^ 2  );

%Diagonal matrix part of the operator (no interactions)
Lop = spdiags( reshape( Lop_kcube, N3, 1 ), 0, N3, N3 );

end
