function rho = deltaPolarIc(systemObj)

% Initialize rho
rho = systemObj.c .* ...
    ones(systemObj.n1,systemObj.n2,systemObj.n3);

% delta function in f
f = zeros(1,systemObj.n3);
f( round(systemObj.n3/2) ) = systemObj.n3 / (2*pi) ;

% Map distribution to a homogeneous system
for i = 1:systemObj.n3
    rho(:,:,i) = rho(:,:,i) .* f(i);
end


end
