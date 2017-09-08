% 
function dRho = dRhoNoise( rho, amp, n1, n2, n3, diffObj )
keyboard
% take sqrt once
sqrtRho = sqrt( rho );
% random values
rand1 =  amp(1) .* ( 0.5 + rand(n1, n2, n3) );
rand2 =  amp(1) .* ( 0.5 + rand(n1, n2, n3) );
% fluxes
j1 = rand1 .* sqrtRho;
j2 = rand2 .* sqrtRho;
% third dimension if needed
if n3 > 1
  rand3 =  amp(2) .* ( 0.5 + rand(n1, n2, n3) );
  j3 = rand3 .* sqrtRho;
end
% take divergence
dRho = diffObj.ik1rep3 .* j1 + diffObj.ik2rep3 .* j2 ...
  + diffObj.ik3rep3 .* j3;
end
