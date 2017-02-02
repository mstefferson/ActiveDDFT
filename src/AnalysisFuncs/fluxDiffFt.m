function [JxDiff, JyDiff, JphiDiff] = ...
  fluxDiffFt( rhoFT, D, diffObj)

% Derivatives in k-space
drho_dxFt = diffObj.ik1rep3 .*  rhoFT;
drho_dyFt = diffObj.ik2rep3 .*  rhoFT;
drho_dphiFt = diffObj.ik3rep3 .* rhoFT;
% Back to real
drho_dx = real( ifftn(ifftshift( drho_dxFt ) ) );
drho_dy = real( ifftn(ifftshift( drho_dyFt ) ) );
drho_dphi = real( ifftn(ifftshift( drho_dphiFt ) ) );
% Calculate fluxes
JxDiff = - ( D.xx .* drho_dx + D.xy .* drho_dy );
JyDiff = - ( D.xy .* drho_dx + D.yy .* drho_dy );
JphiDiff = - ( D.mm .* drho_dphi );
