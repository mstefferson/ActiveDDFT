function [JxDr, JyDr, JphiDr] = fluxDiff( rho, vd, systemOb, cosPhi, sinPhi )

% Flux from driving
JxDr = vd .* cosPhi .* rho;
JyDr = vd .* sinPhi .* rho;
JphiDr = zeros( systemObj.Nx, systemObj.Ny, systemObj.Nm );
