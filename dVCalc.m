function [dV] = dVCalc(vFt,systemObj, diffObj, interObj)

% Allocate and calulate
n3 = systemObj.n3;
[n1v,n2v,n3v] = size(vFt);
%Takes its derivative in k-space, product in real, then back to k-space
ik1 = diffObj.ik1rep3( interObj.ind1, interObj.ind2, interObj.ind3 );
ik2 = diffObj.ik2rep3( interObj.ind1, interObj.ind2, interObj.ind3 );
ik3 = diffObj.ik3rep3( interObj.ind1, interObj.ind2, interObj.ind3 );
% coordinate 1
if n1v > 1
  dVdx1Ft   = ik1 .*  vFt; % derivative of v in k space
  dV.dVdx1   =  real(ifftn(ifftshift(dVdx1Ft))); % back to real
else
  dV.dVdx1 = 0;
end
% coordinate 2
if n2v > 1
  dVdx2Ft   = ik2 .*  vFt;
  dV.dVdx2   =  real(ifftn(ifftshift(dVdx2Ft)));
else
  dV.dVdx2 = 0;
end
% coordinate 3
if n3v > 1
  dVdx3Ft = ik3 .*  vFt;
  dV.dVdx3 =  real(ifftn(ifftshift(dVdx3Ft)));
else
  dV.dVdx3 = 0;
end 
% Store it
