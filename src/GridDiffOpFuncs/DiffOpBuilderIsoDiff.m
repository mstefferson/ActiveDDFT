% Build the diffusion operator
% Currently uncalled

function [Lop] = DiffOpBuilderIsoDiff(diffObj,gridObj,Nx,Ny,Nm,N3)

% Build a strange km for repmat
km = zeros( 1, 1, Nm );
km(1,1,:) = gridObj.km;
%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( ( (diffObj.D_par + diffObj.D_perp)./ 2 ) .* ...
  ( repmat( gridObj.kx', [1, Ny, Nm] ) .^ 2 + ...
  repmat( gridObj.ky, [Nx, 1, Nm ]) .^ 2 ) + ...
  diffObj.D_rot .* repmat( km, [Nx, Ny, 1] ) .^ 2  );

%Diagonal matrix part of the operator (no interactions)
Lop = spdiags( reshape( Lop_kcube, N3, 1 ), 0, N3, N3 );

end
