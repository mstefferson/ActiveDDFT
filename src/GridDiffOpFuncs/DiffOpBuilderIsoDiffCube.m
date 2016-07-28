% Build the diffusion operator
% Currently called by HR2DrotDenEvolverFTBodyIdC

function [Lop_kcube] = DiffOpBuilderIsoDiffCube(DiffMobObj,GridObj,Nx,Ny,Nm)

% Build a strange km for repmat
km = zeros( 1, 1, Nm );
km(1,1,:) = GridObj.km;

%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( DiffMobObj.D_pos  .* ...
  ( repmat( GridObj.kx', [1, Ny, Nm] ) .^ 2 + ...
  repmat( GridObj.ky, [Nx, 1, Nm ]) .^ 2 ) + ...
  DiffMobObj.D_rot .* repmat( km, [Nx, Ny, 1] ) .^2  );

end
