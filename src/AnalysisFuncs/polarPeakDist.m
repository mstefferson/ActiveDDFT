function [polarP2Pdist] = polarPeakDist( F, L )

% only do it if there is a Polar order above 0.1%
maxF = max(F);
if maxF > 1e-3
  % find one max
  [~, maxInd1] = max(F);
  % find the other
  F(maxInd1) = 0;
  [~,maxInd2] = max(F);
  polarP2Pdist = abs( maxInd1 - maxInd2 ) ./ length(F) .* L;
else
  polarP2Pdist = 0;
end

end

