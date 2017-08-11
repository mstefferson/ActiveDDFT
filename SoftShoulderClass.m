classdef SoftShoulderClass
  properties
    Str = '';
    Eps = 1;
    A = 1;
    R = 1;
    Rs = 1;
    Length = 1;
    ReshapeInds = [1 1 1];
    Vv = [];
    VvFt = [];
  end
  
  methods
    % Constructor
    function obj = SoftShoulderClass( str, epsilon, a, R, Rs, n1, l1, n2, l2 )
      if nargin == 9
        obj.Str = str;
        obj.Eps = epsilon;
        obj.A = a;
        obj.R = R;
        obj.Rs = Rs;
        % make potenials
        keyboard
        [obj.Vv, obj.VvFt] = obj.makeV( epsilon, a, R, Rs, n1, l1, n2, l2 );
      else
        fprintf('Error: incorrect number of inputs in constructor\n');
        error('Error: incorrect number of inputs in constructor')
      end
    end
    % make V
    function [v, vFt] = makeV( epsilon, a, R, Rs, n1, l1, n2, l2 )
      % Calcualte distances
      dx1    = l1/n1;
      dx2    = l2/n2;
      x1 = dx1 .* [ 0 : 1 : n1 / 2  -n1/2+1 : 1 : -1];
      x2 = dx2 .* [ 0 : 1 : n2 / 2  -n2/2+1 : 1 : -1];
      [x2m, x1m] = meshgrid( x2, x1 );
      r2 = x1m .^ 2 + x2m .^ 2;
      r8 = r2 .^ 4;
      R8 = R ^ 8;
      Rs8 = Rs ^ 8;
      % soft shoulder potential
      v = epsilon .* (  exp( - r8 ./ R8 ) + a .* exp( -r8 ./ Rs8 ) );
      % FT
      vFt = fftshift( fftn( v ) );
    end
  end
end %class






