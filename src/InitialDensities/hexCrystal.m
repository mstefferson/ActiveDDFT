% crystal = hexCrystal(l, n, c, aWant, sigma)
%
% Builds a periodic hexagonal crystal
% l: box length
% n: grid points
% numPart: number of particles
% aWant: desired crystal spacing
% sigma: peak width
%
function [crystal, gridInfo] = hexCrystal(l, n, numPart, aWant, sigma)
% get crystal lattice vectors
% get spacing
nPeaksCol = round(l/aWant);
a = l / nPeaksCol;
% na = round( a * n / l );
bWant = sqrt(3) * a / 2;
nPeaksRow = round(l/(2*bWant));
nPeaksRow = 2 * nPeaksRow;
b = l / (nPeaksRow);
% in k-space
kaExact = 4 * pi / ( sqrt(3) * a );
nka = round( kaExact * l / (2*pi) ) ;
nkaEven = 2 * round( kaExact * l / (4*pi) ) ;
kaClosest = 2 * pi / l * nka;
kaClosestEven = 2 * pi / l * nkaEven;
% make a crystal or liquid
if isinf(a)
  crystal = numPart / l^2 * ones(n,n);
else
  % build grid
  dx = l / n;
  x = 0:dx:l-dx;
  % mesh grid
  [x2mesh, x1mesh] = meshgrid( x, x );
  colPlace = 0:nPeaksCol;
  rowPlace = 0:nPeaksRow;
  numPeaks = (nPeaksCol+1) * (nPeaksRow+1);
  latticeInds = combvec( rowPlace, colPlace );
  % create lattice basis vectors
  a1 = [ 0; a ];
  a2 = [ b; a/2 ];
  transferMat = [ a2 a1]; % a2 first cause rowPlace in row 1
  crystalPeaks = transferMat * latticeInds;
  crystalPeaksWrap = mod( crystalPeaks, l+a );
  crystalPeaksWrap= crystalPeaksWrap.';
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
end
% store grid info
gridInfo.a = a;
% gridInfo.na = na;
gridInfo.kaExact = kaExact;
gridInfo.kna = nka;
gridInfo.kaClosest = kaClosest;
gridInfo.knaEven = nkaEven;
gridInfo.kaClosestEven = kaClosestEven;
