% Course graining effects. Currently assumes isotropic diffusion 

function dRhoFt = dRhoNoise( rho, amp, n1, n2, n3, diffObj, systemObj, dt )
% take sqrt once
dx = systemObj.l1 / systemObj.n1;
dy = systemObj.l2 / systemObj.n2;
dphi = systemObj.l3 / systemObj.n3;

fluxAmp = sqrt( 24 * amp(1) * diffObj.D_pos * rho / ...
  ( dt * dx * dy * dphi ) );
% random values
rand1 =  -0.5 + rand(n1, n2, n3);
rand2 =  -0.5 + rand(n1, n2, n3);
% fluxes
j1 = fluxAmp .* rand1;
j1Ft = fftshift( fftn( j1 ) ); 
j2 = fluxAmp .* rand2;
j2Ft = fftshift( fftn( j2 ) );
dRhoFt = diffObj.ik1rep3 .* j1Ft + diffObj.ik2rep3 .* j2Ft;
% third dimension if needed
if n3 > 1
  fluxAmp = sqrt( 24 * amp(1) * diffObj.D_rot * rho / ( dt * dx * dy * dphi ) );
  rand3 = -0.5 + rand(n1, n2, n3);
  j3 = fluxAmp .* rand3;
  j3Ft = fftshift( fftn( j3 ) );
  dRhoFt = dRhoFt + diffObj.ik3rep3 .* j3Ft;
end

if 0
  ind1 = 1:n1;
  ind2 = [2:n1 1];
  dx = 10 / n1;
  dRho2 = 1 ./ dx .* ( j1(ind2, 1:n2) - j1(ind1, 1:n2) + ...
  j2(1:n1, ind2) - j2(1:n1, ind1) );
  dRho = real( ifftn( ifftshift( dRhoFt ) ) );
  subplot(2,2,1);
  imagesc( abs(j1+j2))
  colorbar
  subplot(2,2,2);
  imagesc(sqrtRho)
  colorbar
  subplot(2,2,3);
  imagesc(dRho2)
  colorbar
  subplot(2,2,4);
  imagesc(dRho)
  colorbar
  keyboard
end
end
