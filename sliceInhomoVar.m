function [maxVar, rows, cols, posVarLab, Lvar] = bandAnalysis( systemObj, x, y ) 

% Find direction of max variance
rowVar = std( F( :, 1 ) );
colVar = std( F( 1, : ) );

% Take slices
if rowVar >= colVar
  maxVar = rowVar;
  rows = 1:systemObj.Nx;
  cols = 1;
  posVarLab = 'x';
  Lvar = systemObj.Lx;
else
  maxVar = colVar;
  rows = 1;
  cols = 1:systemObj.Ny;
  posVarLab = 'y';
  Lvar = systemObj.Ly;
end

