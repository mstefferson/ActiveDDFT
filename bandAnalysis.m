function [maxF, minF, aveF, diffSc, fwhm] = bandAnalysis( F ) 

Nx = systemObj.Nx;
Ny = systemObj.Ny;
Nm = systemObj.Nm;
x  = gridObj.x;
y  = gridObj.y;

% Find direction of max variance
rowVar = std( F( :, 1 ) );
colVar = std( F( 1, : ) );

% Take slices
if rowVar >= colVar
  maxVar = rowVar;
  rows = 1:Nx;
  cols = 1;
  height = 1:Nm;
  NposVar = Nx;
  posVar = x;
  posVarLab = 'x';
  shiftDim = 1;
  centerPos = Nx / 2;
  Fslice = F(:,1);
  Lvar = Lx;
else
  maxVar = colVar;
  rows = 1;
  cols = 1:Ny;
  height = 1:Nm;
  NposVar = Ny;
  posVar = y;
  posVarLab = 'y';
  shiftDim = 2;
  centerPos = Ny / 2;
  Fslice = F(1,:);
  Lvar = Ly;
end

[maxF, maxFind] = max( Fslice );
[minF, minFind] = max( Fslice );
aveF = (maxF + minF) / 2;
diffSc = ( maxF - minF ) ./ aveF;

% find left half max if max point is larger than N/2
if maxFind >= NposVar / 2
  [~, halfmaxInd]= min( Fslice(1:maxFind) - maxF );
  deltaInd = maxFind - halfmaxInd;
else %right hand max
  [~, halfmaxInd]= min( Fslice(maxFind:end) - maxF );
  deltaInd = halfmaxInd- maxFind ;
end

fwhm = 2 .* deltaInd ./ NposVar .* Lvar;
