classdef DrhoPolarDriveClass
  properties
    Flag = [];
    Vd = [];
    Vx = [];
    Vy = [];
    Ikx = [];
    Iky = [];
  end
  
  methods
    % Constructor
    function obj = DrhoPolarDriveClass( flag, v, n1, n2, n3, ...
        phi, k1m2d, k2m2d )
      % set properties
      obj.Flag = flag;
      if v == 0
        obj.Flag = 0;
      end
      if obj.Flag
        obj.Ikx = sqrt(-1) * k1m2d;
        obj.Iky = sqrt(-1) * k2m2d;
        phiReshape = zeros( 1, 1, n3 );
        phiReshape(1,1,:) = phi;
        cosPhi3 = cos( repmat( phiReshape, [n1, n2, 1] ) );
        sinPhi3 = sin( repmat( phiReshape, [n1, n2, 1] ) );
        obj.Vx = v .* cosPhi3;
        obj.Vy = v .* sinPhi3;
      end
    end % constructor
    % make V
    function [dRho] = calcDrho( obj, rho )
      % calc dRho
      dRho = -obj.Ikx .* fftshift(fftn( rho .* obj.Vx ) ) + ... ;
        -obj.Iky .* fftshift(fftn( rho .* obj.Vy ) );
    end
  end % methods
end %class

