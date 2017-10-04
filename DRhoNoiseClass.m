classdef DRhoNoiseClass
  properties
    PosFluxAmp = [];
    RotFluxAmp = [];
    N1 = [];
    N2 = [];
    N3 = [];
    Ik1rep = [];
    Ik2rep = [];
    Ik3rep = [];
    AngleFluxFlag = 1;
  end
  
  methods
    % Constructor
    function obj = DRhoNoiseClass( amp, n1, n2, n3, l1, l2, Dpos, Drot, ik1, ik2, ik3,  dt )
  % set properties
  obj.N1 = n1;
  obj.N2 = n2;
  obj.N3 = n3;
  obj.Ik1rep = ik1;
  obj.Ik2rep = ik2;
  % grid stuff
  dx = l1 / n1;
  dy = l2 / n2;
  dphi = l3 / systemObj.n3;
  % position
  obj.PosFluxAmp = sqrt( 24 * amp(1) * Dpos  / ...
  ( dt * dx * dy * dphi ) );
  % random values
  if n3 > 1
    obj.AngleFluxFlag = 1;
    obj.RotFluxAmp = sqrt( 24 * amp(1) * Drot  / ...
      ( dt * dx * dy * dphi ) );
    obj.Ik3rep = ik3;
  else
    obj.AngleFluxFlag = 0;
    obj.RotFluxAmp = 0;
  end
end % constructor

    % make V
    function obj = makeV(obj)
      obj.Vv = obj.A * obj.Xv;
      obj.VvReshape = reshape( obj.Vv, obj.ReshapeInds ) ;
    end
    % make Derivative
    function obj = makeDerivative(obj)
      dv = obj.A * ones( 1, obj.Length );
      obj.Dv = dv;
      if obj.Dim == 1
        obj.DvDx1 = reshape( dv, obj.ReshapeInds );
        obj.DvDx2 = 0;
        obj.DvDx3 = 0;
      elseif obj.Dim == 2
        obj.DvDx1 = 0;
        obj.DvDx2 = reshape( dv, obj.ReshapeInds );
        obj.DvDx3 = 0;
      elseif obj.Dim == 3
        obj.DvDx1 = 0;
        obj.DvDx2 = 0;
        obj.DvDx3 = reshape( dv, obj.ReshapeInds );
      end
    end
  end %methods (Static)
end %class






