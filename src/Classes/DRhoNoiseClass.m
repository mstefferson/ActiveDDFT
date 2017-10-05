classdef DRhoNoiseClass
  properties
    Flag  = 0;
    IsotropicDiffusion = 1;
    PosFluxAmp = [];
    RotFluxAmp = [];
    N1 = [];
    N2 = [];
    N3 = [];
    Ik1 = [];
    Ik2 = [];
    Ik3 = [];
    AngleFluxFlag = 1;
  end
  
  methods
    % Constructor
    function obj = DRhoNoiseClass( amp, n1, n2, n3, l1, l2, l3,...
        Dpos, Drot, ik1, ik2, ik3,  dt )
      % set flags
      if amp == 0
        obj.Flag = 0;
      else
        obj.Flag = 1;
        % set properties
        obj.N1 = n1;
        obj.N2 = n2;
        obj.N3 = n3;
        obj.Ik1 = ik1;
        obj.Ik2 = ik2;
        % grid stuff
        dx = l1 / n1;
        dy = l2 / n2;
        dphi = l3 / n3;
        % position
        obj.PosFluxAmp = sqrt( 24 * amp(1) * Dpos  / ...
          ( dt * dx * dy * dphi ) );
        % random values
        if n3 > 1
          obj.AngleFluxFlag = 1;
          obj.RotFluxAmp = sqrt( 24 * amp(2) * Drot  / ...
            ( dt * dx * dy * dphi ) );
          obj.Ik3 = ik3;
        else
          obj.AngleFluxFlag = 0;
          obj.RotFluxAmp = 0;
        end
      end
    end % constructor
    
    %%% other methods %%
    % calculate  dRho
    function dRhoFt = calcDrho(obj,rho)
      % get all variables
      sqrtRho = sqrt(rho);
      jAmp = obj.PosFluxAmp .* sqrtRho;
      rand1 = obj.makeRandField();
      rand2 = obj.makeRandField();
      j1Ft = obj.calcFluxFt( jAmp, rand1);
      j2Ft = obj.calcFluxFt( jAmp, rand2);
      dRhoFt = obj.Ik1 .* j1Ft + obj.Ik2 .* j2Ft;
      % third dimension if needed
      if obj.RotFluxAmp
        rand3 = obj.makeRandField();
        jAmp = obj.RotFluxAmp .* sqrtRho;
        j3Ft = obj.calcFluxFt( jAmp, rand3);
        dRhoFt = dRhoFt + obj.Ik3 .* j3Ft;
      end
    end
    % make rand field
    function randField = makeRandField(obj)
      randField  = -0.5 + rand(obj.N1, obj.N2, obj.N3);
    end
  end %methods
  
  methods (Static)
    % get the ft of flux
    function jft = calcFluxFt( jAmp, randField)
      j = jAmp .* randField;
      jft = fftshift( fftn( j ) );
    end
  end %methods (Static)
end%class
