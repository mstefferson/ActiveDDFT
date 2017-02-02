function [JxDr, JyDr, JphiDr] = ...
  fluxDrive( rho, vd, systemObj, cosPhi, sinPhi )

% Flux from driving
JxDr = vd .* cosPhi .* rho;
JyDr = vd .* sinPhi .* rho;
JphiDr = zeros( systemObj.n1, systemObj.n2, systemObj.n3 );
