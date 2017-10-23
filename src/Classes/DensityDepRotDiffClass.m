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
    function [dRho_dt, dRhoDiff_dt, dRhoInt_dt] = calcDrho( obj, rho, rhoFt, jEx )
      obj.calcDiffNl( rho );
      % "flux" without mobility
      obj.JDiff = -real( ifftn( ifftshift( obj.Ik3 .* rhoFt ) ) );
      obj.JEx = jEx;
      obj.Jft = fftshift( fftn( obj.DrNl .* ( obj.JDiff + obj.JEx ) ) );
      dRhoDiff_dt = -obj.Ik3 .* fftshift( fftn( obj.DrNl .* obj.JDiff  ) );
      dRhoInt_dt = -obj.Ik3 .* fftshift( fftn( obj.DrNl .* obj.JEx  ) );
      dRho_dt = -obj.Ik3 .* obj.Jft;
      
      dRhoDiff_dtIf = real( ifftn( ifftshift( dRhoDiff_dt ) ) );
      dRhoInt_dtIf = real( ifftn( ifftshift( dRhoInt_dt ) ) );
      dRhoTotal_dtIf = real( ifftn( ifftshift( dRho_dt ) ) );

    end
  end %methods
end %class

