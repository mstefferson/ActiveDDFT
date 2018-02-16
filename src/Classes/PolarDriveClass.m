classdef PolarDriveClass < handle
  properties
    Flag = [];
    Fd = [];
    Iota1 = 0; % flux 1 (x) before mobility
    Iota2 = 0;  % flux 2 (y) before mobility  
    IotaFac1 = 0; %factor to hit iota1
    IotaFac2 = 0; %factor to hit iota1
  end
  
  methods
    % Constructor
    function obj = DrhoPolarDriveClass( flag, fd, n3, ...
        phi, k1m2d, k2m2d )
      % set properties
      obj.Flag = flag;
      if fd == 0
        obj.Flag = 0;
        obj.Iota1 = 0; 
        obj.Iota2 = 0;
      end
      if obj.Flag
        obj.Fd = fd;
        phiReshape = reshape( phi, [1 1 n3] );
        cosPhi3 = cos( phiReshape );
        sinPhi3 = sin( phiReshape );
        obj.IotaFac1 = fd .* cosPhi3;
        obj.IotaFac2 = fd .* sinPhi3;
      end
    end % constructor
    % calc preflux iota ( j = D iota )
    function calcIota( obj, rho )
      obj.Iota1 = obj.IotaFac1 .* rho; 
      obj.Iota2 = obj.IotaFac2 .* rho;
    end    % calc dRho
  end % methods
end %class

