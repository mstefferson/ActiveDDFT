function rho = deltaPolarIc(systemObj, shiftAngle)
% Initialize rho
rho = systemObj.c .* ...
  ones(systemObj.n1,systemObj.n2,systemObj.n3);
% delta function in f
f = zeros(1,systemObj.n3);
indPeak = mod( round( shiftAngle * systemObj.n1 / (2*pi) ) ...
  + systemObj.n3/2, systemObj.n3 ) + 1;
f(indPeak) = systemObj.n3 / (2*pi) ;
f = reshape( f, [1 1 systemObj.n3] );
% Map distribution to a homogeneous system
rho = rho .* f;
end
