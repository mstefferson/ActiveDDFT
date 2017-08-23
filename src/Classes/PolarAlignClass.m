classdef PolarAlignClass
  properties
    Str = '';
    CorrStr = '';
    Beta = 1;
    Es1 = 1;
    V = []; % potential
    C2 = [];% pair direct correclation function
    ReshapeInds = [1, 1, 1];
  end
  
  methods
    % Constructor
    function obj = PolarAlignClass( str, corr, es1, kbt, n3, l3 )
      if nargin == 6
        obj.Str = str;
        obj.CorrStr = corr;
        obj.Beta = 1 / kbt;
        obj.Es1 = es1;
        obj.ReshapeInds = [1, 1, n3];
        % make potenials
        [obj] = obj.makeV( n3, l3 );
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
    function [obj] = makeV( obj, n3, l3 )
      % calculate angles
      phi = l3/n3 * (0:n3-1);
      % v =  - J u1 dot u2
      obj.V = -obj.Es1 .* cos( phi );
    end
  end
end %class






