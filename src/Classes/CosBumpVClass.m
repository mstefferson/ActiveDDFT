classdef CosBumpVClass
  properties
    Str = '';
    Dim = 0;
    Es1 = 1;
    N = 1;
    L = [];
    Ls1 = [];
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
    function obj = CosBumpVClass( str, dim, a, l, x0, x )
      if nargin == 6
        obj.Str = str;
        obj.N = length(x);
        obj.Dim = dim;
        obj.Es1 = a;
        obj.Xv = x;
        obj.L = x(end) - 2*x(1) + x(2);
        obj.Center = mod( x0 + obj.L / 2, obj.L) - obj.L / 2;
        obj.Ls1 = l;
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
      obj.Vv = -obj.Es1 * cos( pi / obj.Ls1 * obj.Xv ) .^ 2;
      obj.Vv( obj.Xv > obj.Ls1 / 2 )  = 0;
      obj.Vv( obj.Xv < -obj.Ls1 / 2 )  = 0;
      obj.Vv = circshift(obj.Vv, obj.ShiftAmount);
      obj.VvReshape = reshape( obj.Vv, obj.ReshapeInds ) ;
    end
    % make Derivative
    function obj = makeDerivative(obj)
      % take derivative
      dv = obj.Es1 * pi / obj.Ls1 * cos( pi / obj.Ls1 * obj.Xv ) ...
        .* sin( pi / obj.Ls1 * obj.Xv );
      % set outside of well to zero
      dv( obj.Xv > obj.Ls1 / 2 ) = 0;
      dv( obj.Xv < -obj.Ls1 / 2 ) = 0;
      % shift it
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






