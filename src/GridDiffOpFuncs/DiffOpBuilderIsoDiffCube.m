% Build the diffusion operator
% Currently called by HR2DrotDenEvolverFTBodyIdC

function [Lop_kcube] = DiffOpBuilderIsoDiffCube(diffObj,gridObj,n1,n2,n3)

% Build a strange km for repmat
km = zeros( 1, 1, n3 );
km(1,1,:) = gridObj.km;

%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( diffObj.D_pos  .* ...
  ( repmat( gridObj.k1', [1, n2, n3] ) .^ 2 + ...
  repmat( gridObj.k2, [n1, 1, n3 ]) .^ 2 ) + ...
  diffObj.D_rot .* repmat( km, [n1, n2, 1] ) .^2  );

end
