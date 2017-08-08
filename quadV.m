function [v, vFt, dv] = quadV( a, x )
v = 1 / 2 * a * x .^ 2;
vFt = fftshift( fftn( v, a ) ) ;
dv = - a .* x;
end
  
