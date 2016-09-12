function g = gaussShift( x, Lx, gAmp, gVar , xCenter)

  % fix xCenter in case it's crazy
  xCenter   = mod( xCenter, Lx );
  [~, indCenter] = min( abs( x - xCenter ) );
  indMiddle = floor( length(x) / 2 ) + 1;

  % build gaussian
  g = gAmp .* exp( - ( x- Lx/2) .^ 2 / ( 2 .* gVar .^ 2 ) );
  g = circshift( g, [0 , indCenter - indMiddle] );
  