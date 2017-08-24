classdef PolarExternalVClass
  properties
    Str = '';
    Dim = 0;
    Es1 = 1;
    Phase = 0;
    Length = 1;
    Xv = [];
    Vv = [];
    ReshapeInds = [1 1 1];
    VvReshape = [];
    Dv = [];
    DvDx1 = [];
    DvDx2 = [];
    DvDx3 = [];
  end
  
  methods
    % Constructor
    function obj = PolarExternalVClass( str, es1, phase, x )
      if nargin == 4
        obj.Str = str;
        obj.Length = length(x);
        obj.Dim = 3;
        obj.Es1 = es1;
        obj.Phase = phase;
        obj.Xv = x;
        % reshape inds
        obj.ReshapeInds = [ 1 1 obj.Length ];
        % make potenials
        obj = obj.makeV();
        obj = obj.makeDerivative();
      else
        fprintf('Error: incorrect number of inputs in constructor\n');
        error('Error: incorrect number of inputs in constructor')
      end
    end
    % make V
    function obj = makeV(obj)
      obj.Vv = -obj.Es1 * cos( obj.Xv - obj.Phase ) ;
      obj.VvReshape = reshape( obj.Vv, obj.ReshapeInds ) ;
    end
    % make Derivative
    function obj = makeDerivative(obj)
      dv = obj.Es1 * sin( obj.Xv - obj.Phase );
      obj.Dv = dv;
      obj.DvDx1 = 0;
      obj.DvDx2 = 0;
      obj.DvDx3 = reshape( dv, obj.ReshapeInds );
    end
  end %methods (Static)
end %class






