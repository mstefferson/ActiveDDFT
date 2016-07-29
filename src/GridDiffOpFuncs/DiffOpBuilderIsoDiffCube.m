% Build the diffusion operator
% Currently called by HR2DrotDenEvolverFTBodyIdC

function [Lop_kcube] = DiffOpBuilderIsoDiffCube(diffObj,gridObj,Nx,Ny,Nm)

% Build a strange km for repmat
km = zeros( 1, 1, Nm );
km(1,1,:) = gridObj.km;

%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( diffObj.D_pos  .* ...
  ( repmat( gridObj.kx', [1, Ny, Nm] ) .^ 2 + ...
  repmat( gridObj.ky, [Nx, 1, Nm ]) .^ 2 ) + ...
  diffObj.D_rot .* repmat( km, [Nx, Ny, 1] ) .^2  );

end
