classdef DrhoPolarDriveClass < handle
  properties
    Flag = [];
    Fd = [];
    Dpos = [];
    KbT = [];
    Mob = [];
    Fx = [];
    Fy = [];
    Iota1 = 0; % flux 1 (x) before mobility
    Iota2 = 0;  % flux 2 (y) before mobility  
    IotaFac1 = 0; %factor to hit iota1
    IotaFac2 = 0; %factor to hit iota1
    Ikx = [];
    Iky = [];
  end
  
  methods
    % Constructor
    function obj = DrhoPolarDriveClass( flag, fd, n1, n2, n3, ...
        phi, k1m2d, k2m2d, Dpos, kbT )
      % set properties
      obj.Flag = flag;
      if fd == 0
        obj.Flag = 0;
        obj.Iota1 = 0; 
        obj.Iota2 = 0;
      end
      if obj.Flag
        obj.Fd = fd;
        obj.KbT = kbT;
        obj.Dpos = Dpos;
        obj.Mob = Dpos ./ kbT;
        obj.Ikx = sqrt(-1) * k1m2d;
        obj.Iky = sqrt(-1) * k2m2d;
        phiReshape = zeros( 1, 1, n3 );
        phiReshape(1,1,:) = phi;
        cosPhi3 = cos( repmat( phiReshape, [n1, n2, 1] ) );
        sinPhi3 = sin( repmat( phiReshape, [n1, n2, 1] ) );
        obj.IotaFac1 = fd ./ kbT .* cosPhi3;
        obj.IotaFac2 = fd ./ kbT .* sinPhi3;
      end
    end % constructor

    % calc preflux iota ( j = D iota )
    function calcIota( obj, rho )
      obj.Iota1 = obj.IotaFac1 .* rho; 
      obj.Iota2 = obj.IotaFac2 .* rho;
    end    % calc dRho

    function [dRho] = calcDrho( obj, rho )
      % calculate pre jflux iota
      obj.calcIota( rho )
      % calc dRho
      dRho = -obj.Ikx .* fftshift(fftn( obj.Dpos .* obj.Iota1 ) ) + ... ;
        -obj.Iky .* fftshift(fftn( obj.Dpos .* obj.Iota2 ) );
    end

  end % methods
end %class

