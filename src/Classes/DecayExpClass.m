classdef DecayExpClass
  properties
    Str = '';
    Es1 = 1;
    Ls1 = 1;
    V = [];
    VFt = [];
    ReshapeInds = [1, 1, 1];
 end
  
  methods
    % Constructor
    function obj = DecayExpClass( str, es1, ls1, n1, n2, l1, l2 )
      if nargin == 7
        obj.Str = str;
        obj.Es1 = es1;
        obj.Ls1 = ls1;
        obj.ReshapeInds = [n1, n1, 1];
        % make potenials
        [obj] = obj.makeV( n1, n2, l1, l2 );
      else
        fprintf('Error: incorrect number of inputs in constructor\n');
        error('Error: incorrect number of inputs in constructor')
      end
    end
    % make V
    function [obj] = makeV( obj,  n1, n2, l1, l2 )
    % Calcualte distances
    dx1    = l1/n1;
    dx2    = l2/n2;
    x1 = dx1 .* [ 0 : 1 : n1 / 2  -n1/2+1 : 1 : -1]; 
    x2 = dx2 .* [ 0 : 1 : n2 / 2  -n2/2+1 : 1 : -1];
    [x2m, x1m] = meshgrid( x2, x1 );
    r = sqrt( x1m .^ 2 + x2m .^ 2 );
    % soft shoulder potential
    obj.V = obj.Es1 .* exp( - r ./ obj.Ls1 );
    % FT
    obj.VFt = fftshift( fftn( obj.V ) );
    end
  end
end %class






