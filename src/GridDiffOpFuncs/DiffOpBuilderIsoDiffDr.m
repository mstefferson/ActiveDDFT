% Build the diffusion operator

function [Lop] = DiffOpBuilderIsoDiffDr(diffObj,gridObj,N3)

%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( diffObj.D_pos .* (gridObj.kx3D.^2 + gridObj.ky3D.^2) ...
            + diffObj.D_rot .* gridObj.km3D.^2 );

%Diagonal matrix part of the operator (no interactions)
Lop = spdiags( reshape( Lop_kcube, N3, 1 ), 0, N3, N3 );

end