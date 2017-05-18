% [fm, fmFt] = mayerFncHrLabFrame( n1, n2, n3, l1, l2, lrod )
% 
% Description: calculates mayer function fm(x,y,u1,u2)---and its 
% fourier transform---between two rods 
% with relative distance and orientations u1, u2
%
function [fm, fmFt] = mayerFncHrLabFrame( n1, n2, n3, l1, l2, lRod )
% spatial vectors
dx1 = l1 ./ n1;
dx2 = l2 ./ n2;
maxIndX1 = ceil( n1 * lRod / l1);
maxIndX2 = ceil( n2 * lRod / l2);
indX1 = unique( [ 1:maxIndX1+1 n1-maxIndX1+1:n1] );
indX2 = unique( [ 1:maxIndX2+1 n2-maxIndX2+1:n2] );
x1 = dx1 .* [0:n1/2 -n1/2+1:1:-1];
x2 = dx2 .* [0:n2/2 -n2/2+1:1:-1];
dx1Vec = x1(indX1);
dx2Vec = x2(indX2);
% angles
phi = 0:2*pi/n3:2*pi-2*pi/n3;
cosPhi = cos(phi);  
sinPhi = sin(phi);  
% perpendicular unit vectors
[v1X, v2X] = meshgrid( -sinPhi, -sinPhi );
[v1Y, v2Y] = meshgrid( cosPhi, cosPhi );
% allocate and loop
fm = zeros( n1, n2, n3, n3 );
% loop it
for ii = 1:length( indX1 )
  for jj = 1:length( indX2 )
    indsFmX = indX1(ii);
    indsFmY = indX2(jj);
    overlap = overlapTestFrenkel( dx1Vec(ii), dx2Vec(jj), v1X, v1Y, v2X, v2Y, lRod );
    fm( indsFmX, indsFmY, :, : ) = reshape( -overlap, [1 1 n3 n3] );
  end
end
% fourier transform it
fmFt = fftshift( fftn( fm ) );
