function [cStats, pStats, nStats] = ...
  bandStatsCPNwrap( Cslice, Pslice, Nslice, Lvar, NposVar )
% C slice
[cStats.maxV, cStats.minV, cStats.aveV, ...
  cStats.vdiff, cStats.fwhm, cStats.fwhd, maxInd] = ...
  bandStats( Cslice, Lvar, NposVar ) ;
% Scale c's 
cStats.maxV = cStats.maxV / pi;
cStats.minV = cStats.minV / pi;
cStats.aveV = cStats.aveV / pi;
cStats.vdiff = cStats.vdiff / cStats.aveV;

% P slice: two peaks so be careful
if maxInd > length(Pslice) / 2;
  pInd = 1:maxInd;
else
  pInd = maxInd:length(Pslice);
end

[pStats.maxV, pStats.minV, pStats.aveV, ...
  pStats.vdiff, pStats.fwhm, pStats.fwhd,~] = ...
  bandStats( Pslice(pInd), Lvar, NposVar ) ;

% N slice
[nStats.maxV, nStats.minV, nStats.aveV, ...
  nStats.vdiff, nStats.fwhm, nStats.fwhd,~] = ...
  bandStats( Nslice, Lvar, NposVar ) ;

end

