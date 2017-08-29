% dV = dVCalc(vFt,systemObj, diffObj, interObj)
% calculate the gradients in the potential (forces)
% outputs:
% dV.fx1 = dV/dx1
% dV.fx2 = dV/dx2
% dV.fx3 = dV/dx3

function [dV] = dVCalc(vFt, diffObj, dv1Flag, dv2Flag, dv3Flag, ind1, ind2, ind3)
%Takes its derivative in k-space, product in real, then back to k-space
% coordinate 1
if dv1Flag == 1
  ik1 = diffObj.ik1rep3( ind1, ind2, ind3 );
  dVdx1Ft   = ik1 .*  vFt; % derivative of v in k space
  dV.dx1   =  real(ifftn(ifftshift(dVdx1Ft))); % back to real
else
  dV.dx1 = 0;
end
% coordinate 2
if dv2Flag == 1
  ik2 = diffObj.ik2rep3( ind1, ind2, ind3 );
  dVdx2Ft   = ik2 .*  vFt;
  dV.dx2   =  real(ifftn(ifftshift(dVdx2Ft)));
else
  dV.dx2 = 0;
end
% coordinate 3
if dv3Flag == 1
  ik3 = diffObj.ik3rep3( ind1, ind2, ind3 );
  dVdx3Ft = ik3 .*  vFt;
  dV.dx3 =  real(ifftn(ifftshift(dVdx3Ft)));
else
  dV.dx3 = 0;
end 