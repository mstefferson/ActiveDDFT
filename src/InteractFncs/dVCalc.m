% dV = dVCalc(vFt,systemObj, diffObj, interObj)
% calculate the gradients in the potential (forces)
% outputs:
% dV.fx1 = dV/dx1
% dV.fx2 = dV/dx2
% dV.fx3 = dV/dx3

function [dV] = dVCalc(vFt,systemObj, diffObj, ind1, ind2, ind3)
% Allocate and calulate
n3 = systemObj.n3;
[n1v,n2v,n3v] = size(vFt);
%Takes its derivative in k-space, product in real, then back to k-space
ik1 = diffObj.ik1rep3( ind1, ind2, ind3 );
ik2 = diffObj.ik2rep3( ind1, ind2, ind3 );
ik3 = diffObj.ik3rep3( ind1, ind2, ind3 );
% coordinate 1
if n1v > 1
  dVdx1Ft   = ik1 .*  vFt; % derivative of v in k space
  dV.dx1   =  real(ifftn(ifftshift(dVdx1Ft))); % back to real
else
  dV.dx1 = 0;
end
% coordinate 2
if n2v > 1
  dVdx2Ft   = ik2 .*  vFt;
  dV.dx2   =  real(ifftn(ifftshift(dVdx2Ft)));
else
  dV.dx2 = 0;
end
% coordinate 3
if n3v > 1
  dVdx3Ft = ik3 .*  vFt;
  dV.dx3 =  real(ifftn(ifftshift(dVdx3Ft)));
else
  dV.dx3 = 0;
end 
% Store it
