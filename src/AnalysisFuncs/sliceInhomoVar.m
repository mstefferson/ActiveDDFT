function [maxVar, rows, cols, posVarLab, NposVar, Lvar] = ...
  sliceInhomoVar( systemObj, F ) 

% Find direction of max variance
rowVar = std( F( :, 1 ) );
colVar = std( F( 1, : ) );

% Take slices
if rowVar >= colVar
  maxVar = rowVar;
  rows = 1:systemObj.n1;
  cols = 1;
  posVarLab = 'x';
  NposVar = systemObj.n1;
  Lvar = systemObj.l1;
else
  maxVar = colVar;
  rows = 1;
  cols = 1:systemObj.n2;
  posVarLab = 'y';
  NposVar = systemObj.n2;
  Lvar = systemObj.l2;
end

