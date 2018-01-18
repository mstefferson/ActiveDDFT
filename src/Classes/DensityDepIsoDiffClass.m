% DensityDepIsoDiffClass: handles non-linear density contributions to diffusion
% considering an anisotropic diffusion coefficient
% Can turn on NL diffusion in x,y,phi

classdef DensityDepIsoDiffClass < handle
  properties
    Flag = []; % flag to calculate or not
    FlagPos = []; % flag to calculate position contribution or not
    FlagRot = []; % flag to calculate rotational contribution or not
    Order = []; % 'lin' or 'quad'
    OrderId = []; % 1 = linear, 2 = quad
    N1 = []; % grid points in 1
    N2 = []; % grid points in 2
    N3 = []; % grid points in 3
    DimInclude = []; % dimension to include
    D0Pos = 0; % constant diffusion positional
    D0R = 0; % constant diffusion rotational
    K03 = []; % k = 0 in 3rd dimension
    DnlNVec = []; % vector of gridpoints denoting size of nl diffusion
    RhoMax = []; % density, c, where rho goes to zero
    DNlFactR = 0; % nonlinear factor Dr element of matrix
    DNlFactPos = 0; % nonlinear factor positional element of matrix
    DNlR = 0; % nonlinear Dr element of matrix
    DNlPos = 0; % nonlinear Dpos element of matrix
    Ik = cell(1,3); % sqrt(-1) * k1 vec
    DNlMinR = 0; % mininum value for Dr matrix element
    DNlMinPos = 0; % mininum value for Dpos matrix element
  end
  
  methods
    % Constructor
    % rhoMax: density diffusion goes to zero
    % n1,2,3: number of grid points in 1, 2, 3
    % d0: diffusion constant
    % ik1: sqrt(-1) * k1 vector
    % ik2: sqrt(-1) * k3 vector
    function obj = DensityDepIsoDiffClass( order, rhoMax, posRotFlag, ...
        d0, ik, b, n1, n2, n3)
      if rhoMax == 0
        obj.Flag = 0;
      elseif posRotFlag == 0
        obj.Flag = 0;
      else
        obj.Flag = 1;
        if length(posRotFlag) ~= 2
          posRotFlag = [ 1 1 ];
        end
        % flag
        obj.FlagPos = posRotFlag(1);
        obj.FlagRot = posRotFlag(2);
        obj.Flag = 1;
        obj.Order = order;
        obj.setOrder();
        obj.N1 = n1;
        obj.N2 = n2;
        obj.N3 = n3;
        % convert rhoMax input, in bc, to rho = c / 2pi
        obj.K03 = n3 / 2 + 1;
        obj.RhoMax = rhoMax / b;
        obj.DnlNVec = [n1, n2];
        % store constant diffusion coefficients
        obj.D0Pos = d0(1);
        obj.D0R = d0(2);
        if obj.FlagPos
          obj.DimInclude = 1:2;
          [~,obj.DNlFactPos] = obj.calcNlCoeff( obj.D0Pos );
          obj.DNlMinPos = -obj.D0Pos;
          obj.DNlPos = zeros( obj.DnlNVec );
        end
        if obj.FlagRot
          obj.DimInclude = unique( [obj.DimInclude 3] );
          [~,obj.DNlFactR] = obj.calcNlCoeff( obj.D0R );
          obj.DNlMinR = -obj.D0R;
          obj.DNlR = zeros( obj.DnlNVec );
        end
        nVec = { [n1 1 1], [1 n2 1], [1 1 n3] };
        for ii = obj.DimInclude
          % scale rho by average excluded volume and angle
          obj.Ik{ii} = reshape( ik{ii}, nVec{ii} );
        end
      end
    end % constructor
    
    % Set the nl diffusion order
    function [obj] = setOrder( obj )
      if strcmp( obj.Order, 'lin' )
        obj.OrderId = 1;
        fprintf('NL diffusion linear in density\n')
      elseif strcmp( obj.Order, 'quad' )
        obj.OrderId = 2;
        fprintf('NL diffusion quadratic in density\n')
      else
        obj.OrderId = 1;
        fprintf('Incorrect density order. Setting to linear\n')
      end
    end % setOrder
    
    % calc the nl diffusion coefficient factor based on order
    function [obj, dNlFact] = calcNlCoeff( obj, d0 )
      if obj.OrderId == 1
        dNlFact{1} = -d0 / obj.RhoMax ;
      else
        dNlFact{1} = -2 * d0 / obj.RhoMax ;
        dNlFact{2} = d0 / ( obj.RhoMax .^ 2 );
      end
    end % calcNlCoeff
    
    % Set the nl diffusion coeff
    function [obj] = calcDiffNl( obj, rho )
      obj.DNlR = zeros( obj.DnlNVec );
      obj.DNlPos = zeros( obj.DnlNVec );
      inds2calc = rho < obj.RhoMax;
      for ii = 1:obj.OrderId
        if obj.FlagRot
          obj.DNlR(inds2calc) = obj.DNlR(inds2calc) + ...
            obj.DNlFactR{ii} .* ( rho(inds2calc) .^ ii );
        end
        if obj.FlagPos
          obj.DNlPos(inds2calc) = obj.DNlPos(inds2calc) +...
            obj.DNlFactPos{ii} .* ( rho(inds2calc) .^ ii );
        end
      end
      if obj.FlagRot
        obj.DNlR(~inds2calc) = obj.DNlMinR;
      end
      if obj.FlagPos
        obj.DNlPos(~inds2calc) = obj.DNlMinPos;
      end
    end % calcDiffNL
    
    % calc d rho
    function [dRho_dt] = calcDrho( obj, rhoFt, iotaEx )
      % calculate Dnl
      c = 2*pi / obj.N3 * ifftn(ifftshift( rhoFt(:,:,obj.K03) ) );
      obj.calcDiffNl( c );
      % "flux" without mobility
      iotaTemp = cell( 1, 3 );
      iotaFtTemp = cell( 1, 3 );
      for ii = 1:3
        if any( ii == obj.DimInclude )
          iotaDiffTemp = obj.calcIotaDiff( rhoFt, obj.Ik{ii} );
          iotaTemp{ii} = iotaDiffTemp + iotaEx{ii};
        else
          iotaTemp{ii} = 0;
        end
      end
      % fourier transform of true flux
      if obj.FlagPos
        iotaFtTemp{1} = fftshift( fftn( ...
          obj.DNlPos .* iotaTemp{1} ) );
        iotaFtTemp{2} = fftshift( fftn( ...
          obj.DNlPos .* iotaTemp{2} ) );
      end
      if obj.FlagRot
        iotaFtTemp{3} = fftshift( fftn( ...
          obj.DNlR .* iotaTemp{3} ) );
      end
      % minus divergence of flux
      dRho_dt = 0;
      for ii = obj.DimInclude
        dRho_dt = dRho_dt - obj.Ik{ii} .* ( iotaFtTemp{ii} );
      end
    end % calcDrho
  end %methods
  
  methods (Static)
    % fix densities that are too low
    function dCurrent = fixNegativeDiff( dCurrent, dMin )
      logInds = dCurrent < dMin;
      if isscalar( dMin )
        dCurrent( logInds ) = dMin;
      else
        dCurrent( logInds ) = dMin( logInds );
      end
    end
    
    % calc iota diff
    function [iota] = calcIotaDiff( rhoFtTemp, ik )
      % "flux" without mobility from diffusion
      iota  = -real( ifftn( ifftshift( ik .* rhoFtTemp ) ) );
    end
    
    % calc iota total
    function [iota] = calcIotaTotal(iotaDiff, iotaOther )
      % "flux" without mobility total
      iota = iotaDiff + iotaOther;
    end
  end % static methods
end %class
