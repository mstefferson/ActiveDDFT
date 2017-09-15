function lorenz = lorenzianCalc( systemObj, lorenzParams, gridObj)
% grab some parameters
x1 = gridObj.x1;
x2 = gridObj.x2;
x3 = gridObj.x3;
width1 = lorenzParams.width1;
width2 = lorenzParams.width2;
center1 = lorenzParams.center1;
center2 = lorenzParams.center2;
% get shift inds
xcgIndShift1 = mod( round(systemObj.n1 * center1 / systemObj.l1), systemObj.n1);
xcgIndShift2 = mod( round(systemObj.n2 * center2 / systemObj.l2), systemObj.n2);
% calc lorenzian in 2 or 3d
if systemObj.n3 > 1;
  width3 = lorenzParams.width3;
  center3 = lorenzParams.center3;
  [x2mesh, x1mesh, x3mesh] = meshgrid( x2, x1, x3);
  lorenz = lorenzParams.amp ./ ( pi * ( ( x1mesh / width1 ) .^ 2 ...
    + ( x2mesh / width2 ) .^ 2 ...
    + ( x3mesh / width3 ) .^ 2  + 1 ) );
  lorenz = circshift( circshift( circshift( ...
    lorenz, xcgIndShift1, 1 ), xcgIndShift2, 2 ), xcgIndShift3, 3 );
else
  [x2mesh, x1mesh] = meshgrid( x2, x1 );
  lorenz = lorenzParams.amp ./ ( pi * ( ( x1mesh / width1 ) .^ 2 ...
    + ( x2mesh  / width2 ) .^ 2 + 1 ) );
  lorenz = circshift( circshift( lorenz, xcgIndShift1, 1 ), xcgIndShift2, 2 );
end
end

