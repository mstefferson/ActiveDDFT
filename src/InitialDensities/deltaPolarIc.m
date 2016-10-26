function rho = deltaPolarIc(systemObj)

% Initialize rho
rho = systemObj.c .* ...
    ones(systemObj.Nx,systemObj.Ny,systemObj.Nm);

% delta function in f
f = zeros(1,systemObj.Nm);
f( round(systemObj.Nm/2) ) = systemObj.Nm / (2*pi) ;

% Map distribution to a homogeneous system
for i = 1:systemObj.Nm
    rho(:,:,i) = rho(:,:,i) .* f(i);
end


end
