classdef PolarAlignClass
  properties
    Str = '';
    Es1 = 1;
    V = [];
    VFt = [];
    ReshapeInds = [1, 1, 1];
 end
  
  methods
    % Constructor
    function obj = PolarAlignClass( str, es1, n3, l3 )
      if nargin == 4
        obj.Str = str;
        obj.Es1 = es1;
        obj.ReshapeInds = [1, 1, n3];
        % make potenials
        [obj] = obj.makeV( n3, l3 );
      else
        fprintf('Error: incorrect number of inputs in constructor\n');
        error('Error: incorrect number of inputs in constructor')
      end
    end
    % make V
    function [obj] = makeV( obj, n3, l3 )
    % calculate angles
    phi = l3/n3 * (0:n3-1);
    % v =  - J u1 dot u2
    obj.V = -obj.Es1 .* cos( phi );
    % Ft
    obj.VFt = fftshift( fftn( obj.V ) );
    end
  end
end %class






