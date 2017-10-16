classdef DensityDepRotDiffClass < handle
  properties
    Flag = [];
    Dr0 = [];
    RhoMax = [];
    DrNlFact = [];
    DrNl = [];
    Ik3 = [];
    JEx = [];
    JDiff = [];
    Jft = [];
  end
  
  methods
    % Constructor
    function obj = DensityDepRotDiffClass( rhoMax, n1, n2, n3, dr0, ik3 )
      if rhoMax == 0
        obj.Flag = 0;
      else
        obj.Flag = 1;
        obj.Dr0    = dr0;
        obj.RhoMax = rhoMax;
        obj.DrNlFact = -obj.Dr0 / rhoMax;
        obj.Ik3 = reshape( ik3, [1 1 n3] );
        obj.DrNl = zeros(n1,n2,n3);
        obj.JEx = zeros(n1,n2,n3);
        obj.JDiff = zeros(n1,n2,n3);
        obj.Jft = zeros(n1,n2,n3);
      end
    end
    function [obj] = setExcessJ( obj, j )
      obj.JEx = j;
    end
    
    % Set the nl diffusion coeff
    function [obj] = calcDiffNl( obj, rho )
      obj.DrNl = obj.DrNlFact .* rho;
      obj.DrNl( obj.DrNl <  -obj.Dr0 ) = -obj.Dr0;
    end
    
    % calc d rho
    function [dRho] = calcDrho( obj, rho, rhoFt, jEx )
      obj.calcDiffNl( rho );
      % "flux" without mobility
      obj.JDiff = -real( ifftn( ifftshift( obj.Ik3 .* rhoFt ) ) );
      obj.JEx = jEx;
      obj.Jft = fftshift( fftn( obj.DrNl .* ( obj.JDiff + obj.JEx ) ) );
      dRho = - obj.Ik3 .* obj.Jft;
      keyboard
    end
  end %methods
end %class

