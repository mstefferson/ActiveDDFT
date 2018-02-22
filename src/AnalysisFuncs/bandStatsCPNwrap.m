function [cStats, pStats, nStats] = ...
  bandStatsCPNwrap( Cslice, Pslice, Nslice, Lvar, NposVar )
% C slice
[cStats.maxV, cStats.minV, cStats.aveV, ...
  cStats.vdiff, cStats.fwhd, maxInd] = ...
  bandStats( Cslice, Lvar, NposVar ) ;
% P slice: two peaks so be careful
 deltaInd = round( cStats.fwhd / Lvar .* NposVar);
if maxInd > length(Pslice) / 2
  pInd = max(1, maxInd-deltaInd):maxInd;
else
  pInd = maxInd:min( maxInd+deltaInd, NposVar );
end
[pStats.maxV, pStats.minV, pStats.aveV, ...
  pStats.vdiff, pStats.fwhd,~] = ...
  bandStats( Pslice(pInd), Lvar, NposVar ) ;
pStats.p2p = polarPeakDist( Pslice, Lvar );

% N slice
[nStats.maxV, nStats.minV, nStats.aveV, ...
  nStats.vdiff, nStats.fwhd,~] = ...
  bandStats( Nslice, Lvar, NposVar ) ;
end

