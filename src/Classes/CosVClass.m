classdef CosVClass
  properties
    Str = '';
    Dim = 0;
    A = 1;
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
    function obj = CosVClass( str, dim, a, x0, x )
      if nargin == 5
        obj.Str = str;
        obj.N = length(x);
        obj.Dim = dim;
        obj.A = a;
        obj.Xv = x;
        obj.Center = x0;
        obj.L = x(end) - 2*x(1) + x(2);
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
      x = 2 * pi / obj.N * ( -obj.N/2:obj.N/2-1 );
      obj.Vv = obj.A * cos( x/2 ) .^ 2;
      obj.Vv = circshift(obj.Vv, obj.ShiftAmount);
      obj.VvReshape = reshape( obj.Vv, obj.ReshapeInds ) ;
    end
    % make Derivative
    function obj = makeDerivative(obj)
      x = 2 * pi / obj.N * ( -obj.N/2:obj.N/2-1 );
      dv = obj.A * cos( x/2 ) .* sin( x/2 );
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






