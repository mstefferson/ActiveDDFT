function [cStats, pStats, nStats] = bandStatsCPNwrap( Cslice, Pslice, Nslice );

% C slice
[cStats.maxV, cStats.minV, cStats.aveV, cStats.vdiff, cStats.fwhm, maxInd] = ...
  bandAnalysis( Cslice ) ;

% P slice: two peaks so be careful
if maxInd > length(Pslice) / 2;
  pInd = 1:maxInd;
else
  pInd = maxInd:length(Pslice);
end
[pStats.maxV, pStats.minV, pStats.aveV, pStats.vdiff, pStats.fwhm, ~] = ...
  bandAnalysis( Pslice(pInd) ) ;

% N slice
[nStats.maxV, nStats.minV, nStats.aveV, nStats.vdiff, nStats.fwhm, ~] = ...
  bandAnalysis( Nslice ) ;

end

