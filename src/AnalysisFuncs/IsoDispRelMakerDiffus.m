function [lambda] = IsoDispRelMakerDiffus(kx,ky,km,D_pos,D_rot)

lambda = - ( D_pos .* ( kx.^2 + ky.^2 ) + D_rot .* km.^2  ) ;

end