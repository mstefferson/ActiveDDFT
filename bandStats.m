% Feed it a slice (1D vec)
function [maxF, minF, aveF, diffMaxMin, fwhm, fwhd, maxInd] = ...
  bandStats( F, Lvar, NposVar ) 
%%
[maxF, maxInd] = max( F );
[minF, ~] = min( F );
aveF = (maxF + minF) / 2;
diffMaxMin = ( maxF - minF );
lF2 = length(F) / 2;

% find left half max if max point is larger than N/2
if maxInd >= lF2
  % half max
  [~, halfmaxIndhm]= min( abs( F( maxInd - lF2 +1 : maxInd ) - (maxF) / 2 ) );
  deltaIndhm = lF2 - halfmaxIndhm;
  % half delta
  [~, halfmaxIndhd]= min( abs( F( maxInd - lF2 +1 : maxInd ) - (maxF - minF) / 2 ) );
  deltaIndhd = lF2 - halfmaxIndhd;
else %right hand max /2
   % half max
  [~, halfmaxIndhm]= min( abs( F( maxInd : maxInd + lF2 ) - (maxF) / 2 ) );
  deltaIndhm = halfmaxIndhm - 1;
    % half delta
     [~, halfmaxIndhd]= min( abs( F( maxInd : maxInd + lF2 ) - (maxF - minF) / 2 ) );
  deltaIndhd = halfmaxIndhd - 1;
end


% want full width
fwhm = 2 .* deltaIndhm ./ NposVar .* Lvar;
fwhd = 2 .* deltaIndhd ./ NposVar .* Lvar;

