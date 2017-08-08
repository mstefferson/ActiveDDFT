function [v, vFt, dv] = linearV( a, x )
v = a * x;
vFt = fftshift( fftn( v, a ) ) ;
dv = a * ones( size(x) );
end
  
