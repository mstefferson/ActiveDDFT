classdef PolarAlignGaussClass
  properties
    Str = '';
    CorrStr = '';
    Beta = 1;
    Es1 = 1;
    Ls1 = 1;
    V = []; % potential
    C2 = [];% pair direct correclation function
    ReshapeInds = [1, 1, 1];
 end
  
  methods
    % Constructor
    function obj = PolarAlignGaussClass( str, corr, es1, ls1, kbt,...
    n1, n2, n3, l1, l2, l3 )
      if nargin == 11
        obj.Str = str;
        obj.CorrStr = corr;
        obj.Beta = 1 / kbt;
        obj.Es1 = es1;
        obj.Ls1 = ls1;
        obj.ReshapeInds = [n1, n2, n3];
        % make potenials
        [obj] = obj.makeV( n1, n2, n3, l1, l2, l3 );
      else
        fprintf('Error: incorrect number of inputs in constructor\n');
        error('Error: incorrect number of inputs in constructor')
      end
      % pair direct correlation
      if strcmp(corr,'mf')
        obj.C2 = -obj.Beta * obj.V;
      elseif strcmp(corr, 'vir')
        obj.C2 = exp( -obj.Beta * obj.V ) - 1;
      else
        fprintf('In soft shoulder. Do not recognize correlation str. Going MF\n')
        obj.C2 = -obj.Beta * obj.V;
      end
    end
    % make V
    function [obj] = makeV( obj,  n1, n2, n3, l1, l2, l3 )
      % calculate distances
      dx1    = l1/n1;
      dx2    = l2/n2;
      x1 = dx1 .* [ 0 : 1 : n1 / 2  -n1/2+1 : 1 : -1]; 
      x2 = dx2 .* [ 0 : 1 : n2 / 2  -n2/2+1 : 1 : -1];
      [x2m, x1m] = meshgrid( x2, x1 );
      r2 = x1m .^ 2 + x2m .^ 2;
      % calulate angle
      phi = l3/n3 * (0:n3-1);
      cosPhi = reshape( cos( phi), [ 1 1 n3 ] );
      % polar alignedment with gaussian dropoff
      obj.V = obj.Es1 .*  cosPhi .* exp( - r2 / (2 * obj.Ls1 ^2 ) );
    end
  end
end %class






