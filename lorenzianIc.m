function rho = lorenzianIc(systemObj,gridObj)

x1 = gridObj.x1;
x2 = gridObj.x2;
x3 = gridObj.x3;
width1 = 0.5;
width2 = 0.5;
width3 = 0.5;
shift1 = 0;
shift2 = 0;
shift3 = pi/4;


if systemObj.n3 > 1
[x2mesh, x1mesh, x3mesh] = meshgrid( x2, x1, x3);
rho = 1 ./ ( ( (x1mesh - shift1) / width1 ) .^ 2 ...
  + ( (x2mesh - shift2) / width2 ) .^ 2 ...
  + ( (x3mesh - shift3) / width3 ) .^ 2  + 1 );
else
[x2mesh, x1mesh] = meshgrid( x2, x1 );
rho = 1 ./ ( ( (x1mesh - shift1) / width1 ) .^ 2 ...
  + ( (x2mesh - shift2) / width2 ) .^ 2 + 1 );

end

end

