% Feed it a slice (1D vec)
function [maxF, minF, aveF, diffMaxMin, fwhd, maxInd] = ...
  bandStats( F, Lvar, NposVar )

%%
[maxF, maxInd] = max( F );
[minF, ~] = min( F );
aveF = (maxF + minF) / 2;
diffMaxMin = ( maxF - minF );
lF2 = ceil( length(F) / 2);

% find left half max if max point is larger than N/2
if maxInd >= lF2
  % half delta
  [~, halfmaxIndhd]= min( abs( F( maxInd - lF2 +1 : maxInd ) - (maxF + minF) / 2 ) );
  deltaIndhd = lF2 - halfmaxIndhd;
else %right hand max /2
  % half delta
  [~, halfmaxIndhd]= min( abs( F( maxInd : maxInd + lF2 ) - (maxF + minF) / 2 ) );
  deltaIndhd = halfmaxIndhd - 1;
end


% want full width
fwhd = 2 .* deltaIndhd ./ NposVar .* Lvar;

if fwhd > Lvar / 2
  fwhd
  keyboard
end