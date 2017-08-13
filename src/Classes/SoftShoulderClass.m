classdef SoftShoulderClass
  properties
    Str = '';
    Es1 = 1;
    Es2 = 1;
    Ls1 = 1;
    Ls2 = 1;
    Length = 1;
    V = [];
    VFt = [];
    ReshapeInds = [ 1 1 1];
  end
  
  methods
    % Constructor
    function obj = SoftShoulderClass( str, es1, es2, ls1, ls2, n1, l1, n2, l2 )
      if nargin == 9
        obj.Str = str;
        obj.Es1 = es1;
        obj.Es2 = es2;
        obj.Ls1 = ls1;
        obj.Ls2 = ls2;
        obj.ReshapeInds = [ n1, n2, 1];
        % make potenials
        [obj] = obj.makeV( n1, l1, n2, l2 );
      else
        fprintf('Error: incorrect number of inputs in constructor\n');
        error('Error: incorrect number of inputs in constructor')
      end
    end
    % make V
    function [obj] = makeV( obj, n1, l1, n2, l2 )
      % Calcualte distances
      dx1    = l1/n1;
      dx2    = l2/n2;
      x1 = dx1 .* [ 0 : 1 : n1 / 2  -n1/2+1 : 1 : -1];
      x2 = dx2 .* [ 0 : 1 : n2 / 2  -n2/2+1 : 1 : -1];
      [x2m, x1m] = meshgrid( x2, x1 );
      r2 = x1m .^ 2 + x2m .^ 2;
      r8 = r2 .^ 4;
      R8 = obj.Ls1 ^ 8;
      Rs8 = obj.Ls2 ^ 8;
      % soft shoulder potential
      obj.V = obj.Es1 .* (  exp( - r8 ./ R8 ) + obj.Es2 .* exp( -r8 ./ Rs8 ) );
      % FT
      obj.VFt = fftshift( fftn( obj.V ) );
    end
  end
end %class
