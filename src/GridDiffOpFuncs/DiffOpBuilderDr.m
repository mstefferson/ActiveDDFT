% Build the diffusion operator
% Called by HR2DrotDenEvolverFTBody

function [Lop] = DiffOpBuilderDr(diffObj,gridObj,n1,n2,n3,N2,N3)

% Build a strange km for repmat
k3 = zeros( 1, 1, n3 );
k3(1,1,:) = gridObj.k3;
%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( ( (diffObj.D_par + diffObj.D_perp)./ 2 ) .* ...
  ( repmat( gridObj.k1', [1, n2, n3] ) .^ 2 + ...
  repmat( gridObj.k2, [n1, 1, n3 ]) .^ 2 ) + ...
  diffObj.D_rot .* repmat( k3, [n1, n2, 1] ) .^ 2  );

%Diagonal matrix part of the operator (no interactions)
Lop = spdiags( reshape( Lop_kcube, N3, 1 ), 0, N3, N3 );

%%%%%%%%%%%%%%% COMPLEX CONJUGATION!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%
%Handle the cross terms of the PDE for +/- 2 coupling
if diffObj.Ani == 1
  CpMplus2Vec  = reshape( repmat( diffObj.CfMplus2,  [1 1 n3]), 1, N3 );
  CpMminus2Vec = reshape( repmat( diffObj.CfMminus2, [1 1 n3]), 1, N3 );
  
  % Put them in the matrix
  Mplus2Mtx  = spdiags( [ zeros(1,2*N2) CpMplus2Vec( 1:(N2*(n3-2)) ) ].', 2*N2, N3, N3 );
  Mminus2Mtx = spdiags( [ CpMminus2Vec(1:(N2*(n3-2))) zeros(1,2*N2) ].' , -2*N2, N3, N3 );
  
  %Periodic terms
  %m = N-1, N coupling to m = 1,2
  Mplus2Mtx =  Mplus2Mtx  + ...
    spdiags( [ CpMplus2Vec( (N2*(n3-2)+1):N3 ) zeros( 1, N2*(n3-2) ) ].' , -(N2*(n3-2)), N3, N3 );
  %m = 1,2 coupling to m = N-1,N
  Mminus2Mtx = Mminus2Mtx + ...
    spdiags( [ zeros( 1, N2*(n3-2) ) CpMminus2Vec( (N2*(n3-2)+1):N3 ) ].', N2*(n3-2), N3, N3 ) ;
  Lop = Lop + Mplus2Mtx + Mminus2Mtx;
end

end
