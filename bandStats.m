% Feed it a slice (1D vec)
function [maxF, minF, aveF, diffSc, fwhm, maxInd] = bandAnalysis( F ) 

[maxF, maxFind] = max( F );
[minF, minFind] = max( F );
aveF = (maxF + minF) / 2;
diffSc = ( maxF - minF ) ./ aveF;

% find left half max if max point is larger than N/2
if maxFind >= length(F) / 2
  [~, halfmaxInd]= min( F(1:maxFind) - maxF );
  deltaInd = maxFind - halfmaxInd;
else %right hand max
  [~, halfmaxInd]= min( F(maxFind:end) - maxF );
  deltaInd = halfmaxInd- maxFind ;
end

fwhm = 2 .* deltaInd ./ NposVar .* Lvar;
