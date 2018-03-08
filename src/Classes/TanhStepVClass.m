classdef TanhStepVClass
  properties
    Str = '';
    Dim = 0;
    Es1 = 1;
    Ls1 = 1;
    Ls2 = 0;
    Sig = 1;
    N = 1;
    L = [];
    Center = [];
    ShiftAmount = 0;
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
    function obj = TanhStepVClass( str, dim, es1, ls1, x0, steepness, x )
      if nargin > 6
        obj.Str = str;
        obj.N = length(x);
        obj.Dim = dim;
        obj.Es1 = es1;
        obj.Ls1 = ls1;
        obj.Xv = x;
        obj.L = x(end) - 2*x(1) + x(2);
        if isempty(steepness)
          obj.Ls2 = obj.L / 100;
        else
          obj.Ls2 = steepness;
        end
        obj.Center = mod( x0 + obj.L / 2, obj.L) - obj.L / 2;
        obj.ShiftAmount = round( obj.N * obj.Center / obj.L );
        % reshape inds
        if dim == 1
          obj.ReshapeInds = [ obj.N 1 1 ];
        elseif dim == 2
          obj.ReshapeInds = [ 1 obj.N 1 ];
        elseif dim == 3
          obj.ReshapeInds = [ 1 1 obj.N ];
        end
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
      obj.Vv = -obj.Es1 * ( tanh( (obj.Xv + obj.Ls1/2) / obj.Ls2 ) - ...
      tanh( (obj.Xv - obj.Ls1/2) / obj.Ls2 ) );
      obj.Vv = circshift(obj.Vv, obj.ShiftAmount);
      obj.VvReshape = reshape( obj.Vv, obj.ReshapeInds ) ;
    end
    % make Derivative
    function obj = makeDerivative(obj)
      dv = -obj.Es1 / obj.Ls2 * ( ...
      -tanh( (obj.Xv + obj.Ls1/2) / obj.Ls2 ) .^ 2 ...
      + tanh( (obj.Xv - obj.Ls1/2) / obj.Ls2 ) .^2 );
      dv = circshift(dv, obj.ShiftAmount);
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






