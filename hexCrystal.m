% crystal = hexCrystal(l, n, c, aWant, sigma)
% 
% Builds a periodic hexagonal crystal
% l: box length
% n: grid points
% numPart: number of particles
% aWant: desired crystal spacing
% sigma: peak width

function crystal = hexCrystal(l, n, numPart, aWant, sigma)
% build grid
dx = l / n;
x = 0:dx:l-dx;
% mesh grid
[x2mesh, x1mesh] = meshgrid( x, x );
% get crystal lattice vectors
% get spacing
nPeaksCol = round(l/aWant);
a = l / nPeaksCol;
bWant = sqrt(3) * a / 2;
nPeaksRow = round(l/bWant);
b = l / nPeaksRow;
colPlace = 0:nPeaksCol-1;
rowPlace = 0:nPeaksRow-1;
numPeaks = nPeaksRow * nPeaksCol;
latticeInds = combvec( rowPlace, colPlace );
% create lattice basis vectors
a1 = [ 0; a ];
a2 = [ b; a/2 ];
transferMat = [ a2 a1]; % a2 first cause rowPlace in row 1
crystalPeaks = transferMat * latticeInds;
crystalPeaksWrap = mod( crystalPeaks, l );
crystalPeaksWrap= crystalPeaksWrap.';
% shift it off edges
crystalPeaksWrap(:,1) = crystalPeaksWrap(:,1) + b/4;
crystalPeaksWrap(:,2) = crystalPeaksWrap(:,2) + a/4;
% loop over centers
crystal = zeros( n, n );
for ii = 1:numPeaks
  crystal = crystal + exp( -( x1mesh - crystalPeaksWrap(ii,1) ) .^2 / (2 *sigma .^2 ) ...
  -( x2mesh - crystalPeaksWrap(ii,2) ) .^2 / (2 *sigma .^2 )  ); 
end
% get correct size
crystal = crystal(1:n,1:n);
% normalize
currNorm = trapz_periodic( x,trapz_periodic( x, crystal, 1 ), 2 );
crystal  = numPart / currNorm * crystal;
