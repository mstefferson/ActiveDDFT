% Feed it a slice (1D vec)
function [maxF, minF, aveF, diffMaxMin, fwhd, maxInd] = ...
  bandStats( F, Lvar, NposVar )

% turn it into a row vector
[r, c] = size( F );
if r > c
  c = r;
  F = F';
end
% find max then shift it to center
[maxF, maxInd] = max( F );
[minF, ~] = min( F );
aveF = (maxF + minF) / 2;
diffMaxMin = ( maxF - minF );
lF2 = ceil( c / 2);
F = circshift( F, lF2 - maxInd , 2 );
% find half max and delta ind
[~, halfmax]= min( abs( F - (maxF + minF) / 2 ) );
deltaIndhd = abs( halfmax - lF2 );
% want full width
fwhd = 2 .* deltaIndhd  .* Lvar./ NposVar;
