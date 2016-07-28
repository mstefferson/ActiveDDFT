% Build the diffusion operator

function [Lop] = DiffOpBuilderIsoDiffDr(DiffMobObj,GridObj,N3)

%%%%%%%%%%%%%%%%%%Diagonal operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lop_kcube = -( DiffMobObj.D_pos .* (GridObj.kx3D.^2 + GridObj.ky3D.^2) ...
            + DiffMobObj.D_rot .* GridObj.km3D.^2 );

%Diagonal matrix part of the operator (no interactions)
Lop = spdiags( reshape( Lop_kcube, N3, 1 ), 0, N3, N3 );

end