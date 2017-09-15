function g = gaussShift( x, l1, gAmp, gVar , xCenter)

  % fix xCenter in case it's crazy
  xCenter   = mod( xCenter, l1 );
  [~, indCenter] = min( abs( x - xCenter ) );
  indMiddle = floor( length(x) / 2 ) + 1;

  % build gaussian
  g = gAmp .* exp( - ( x- l1/2) .^ 2 / ( 2 .* gVar .^ 2 ) );
  g = circshift( g, [0 , indCenter - indMiddle] );
  