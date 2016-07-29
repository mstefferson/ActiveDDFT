% Build the diffusion operator
% Called by HR2DrotDenEvolverFTBody

function [Lop] = DiffOpBuilderDr(diffObj,gridObj,Nx,Ny,Nm,N2,N3)

% Build a strange km for repmat
km = zeros( 1, 1, Nm );
km(1,1,:) = gridObj.km;
%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( ( (diffObj.D_par + diffObj.D_perp)./ 2 ) .* ...
  ( repmat( gridObj.kx', [1, Ny, Nm] ) .^ 2 + ...
  repmat( gridObj.ky, [Nx, 1, Nm ]) .^ 2 ) + ...
  diffObj.D_rot .* repmat( km, [Nx, Ny, 1] ) .^ 2  );

%Diagonal matrix part of the operator (no interactions)
Lop_DiagMtx = spdiags( reshape( Lop_kcube, N3, 1 ), 0, N3, N3 );

%Handle the cross terms of the PDE for +/- 1 coupling
CpMplus1Vec  = reshape( repmat( diffObj.CfMplus1,  [1 1 Nm]), 1, N3 );
CpMminus1Vec = reshape( repmat( diffObj.CfMminus1, [1 1 Nm]), 1, N3 );


%Handle the cross terms of the PDE for +/- 2 coupling
CpMplus2Vec  = reshape( repmat( diffObj.CfMplus2,  [1 1 Nm]), 1, N3 );
CpMminus2Vec = reshape( repmat( diffObj.CfMminus2, [1 1 Nm]), 1, N3 );

% m +/- 1 terms
% Put them in the matrix
Mplus1Mtx  = spdiags( [ zeros(1,N2) CpMplus1Vec( 1:(N2*(Nm-1)) ) ]', N2, N3, N3 );
Mminus1Mtx = spdiags( [ CpMminus1Vec(1:(N2*(Nm-1))) zeros(1,N2) ]' , -N2, N3, N3 );

%Periodic terms
Mplus1Mtx =  Mplus1Mtx  + ...
    spdiags( [ CpMplus1Vec( (N2*(Nm-1)+1) :N3 ) zeros( 1, N2*(Nm-2) ) ]' , -(N2*(Nm-1)), N3, N3 ); %m = N coupling to m = 2
Mminus1Mtx = Mminus1Mtx + ...
    spdiags( [ zeros( 1, N2*(Nm-1) ) CpMminus1Vec( 1:N2 ) ]', N2*(Nm-1), N3, N3 ) ;                  %m = 1 coupling to m = N-1

% keyboard
% m +/- 2 terms
% Put them in the matrix
Mplus2Mtx  = spdiags( [ zeros(1,2*N2) CpMplus2Vec( 1:(N2*(Nm-2)) ) ]', 2*N2, N3, N3 );
Mminus2Mtx = spdiags( [ CpMminus2Vec(1:(N2*(Nm-2))) zeros(1,2*N2) ]' , -2*N2, N3, N3 );

%Periodic terms
Mplus2Mtx =  Mplus2Mtx  + ...
    spdiags( [ CpMplus2Vec( (N2*(Nm-2)+1):N3 ) zeros( 1, N2*(Nm-2) ) ]' , -(N2*(Nm-2)), N3, N3 ); %m = N coupling to m = 2
Mminus2Mtx = Mminus2Mtx + ...
    spdiags( [ zeros( 1, N2*(Nm-2) ) CpMminus2Vec( 1:2*N2 ) ]', N2*(Nm-2), N3, N3 ) ;                  %m = 1 coupling to m = N-1

%%%%%%%%%%%%%%%%%%%%%%Non-interacting Propagator%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop = Lop_DiagMtx + Mplus1Mtx + Mminus1Mtx + Mplus2Mtx + Mminus2Mtx;

end
