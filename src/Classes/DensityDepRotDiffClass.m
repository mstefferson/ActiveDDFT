classdef DensityDepRotDiffClass < handle
  properties
    Flag = []; % flag to calculate or not
    Dr0 = []; % constant diffusion coefficient
    RhoMax = []; % density, c, where rho goes to zero
    DrNlFact = []; % non linear diffusion coefficient factor
    DrNl = []; % non linear diffusion contribution
    Ik3 = []; % sqrt(-1) * k1 vec
    JEx3 = []; % flux excess (without diff coeff), in dir 3
    JDiff3 = []; % flux diffusion (without diff coeff), in dir 3
    JFt3 = []; % flux total in dir 3
  end
  
  methods
    % Constructor
    % rhoMax: density diffusion goes to zero
    % n1,2,3: number of grid points in 1, 2, 3
    % dr0: rot diffusion constant
    % ik3: sqrt(-1) * k vector
    function obj = DensityDepRotDiffClass( rhoMax, b, n1, n2, n3, dr0, ik3 )
      if rhoMax == 0
        obj.Flag = 0;
      else
        obj.Flag = 1;
        obj.Dr0    = dr0;
        % scale rho by average excluded volume and angle
        obj.RhoMax = rhoMax / b / (2*pi);
        obj.DrNlFact = -obj.Dr0 / obj.RhoMax ;
        obj.Ik3 = reshape( ik3, [1 1 n3] );
        obj.DrNl = zeros(n1,n2,n3);
        obj.JEx3 = zeros(n1,n2,n3);
        obj.JDiff3 = zeros(n1,n2,n3);
        obj.JFt3 = zeros(n1,n2,n3);
      end
    end
    function [obj] = setExcessJ( obj, j )
      obj.JEx3 = j;
    end
    
    % Set the nl diffusion coeff
    function [obj] = calcDiffNl( obj, rho )
      obj.DrNl = obj.DrNlFact .* rho;
      obj.DrNl( obj.DrNl <  -obj.Dr0 ) = -obj.Dr0;
    end
    
    % calc d rho
    function [dRho_dt] = calcDrho( obj, rho, rhoFt, jEx )
      obj.calcDiffNl( rho );
      % "flux" without mobility
      obj.JDiff3 = -real( ifftn( ifftshift( obj.Ik3 .* rhoFt ) ) );
      obj.JEx3 = jEx;
      obj.JFt3 = fftshift( fftn( obj.DrNl .* ( obj.JDiff3 + obj.JEx3 ) ) );
      dRho_dt = -obj.Ik3 .* obj.JFt3;
   end
  end %methods
end %class

