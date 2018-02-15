%rho DensityDepAnisoDiffClass: handles non-linear density contributions to diffusion
% considering an anisotropic diffusion coefficient
% Currently, only alters perpendicular and rotational diffusion
%
classdef DensityDepAnisoDiffClass < handle
  properties
    Flag = []; % flag to calculate or not
    FlagPos = []; % flag to calculate or not
    FlagRot = []; % flag to calculate or not
    Type = []; % 'lin' or 'quad'
    TypeId = []; % 1 = linear, 2 = quad
    PosTerms =[]; % number of coefficients in position terms
    RotTerms = []; % number of coefficients in position terms
    N1 = []; % grid points in 1
    N2 = []; % grid points in 2
    N3 = []; % grid points in 3
    K03 = []; % k = 0 in 3rd dimension
    DnlNVec = []; % vector of gridpoints denoting size of nl diffusion
    RhoMax = []; % density, c, where rho goes to zero
    DimInclude = []; % dimension to include
    Ik = cell(1,3); % sqrt(-1) * k1 vec
    D0Perp = 0; % constant diffusion in perpendicular
    D0R = 0; % constant diffusion rotational
    DNlFactPerp = 0; % total nonlinear contribution perp
    DNlFactR = 0; % nonlinear factor Dr element of matrix
    DNlFact11 = 0; % nonlinear factor D11 element of matrix
    DNlFact12 = 0; % nonlinear factor D12 element of matrix
    DNlFact22 = 0; % nonlinear factor D22 element of matrix
    DNlR = 0; % nonlinear Dr element of matrix
    DNl11 = 0; % nonlinear D11 element of matrix
    DNl12 = 0; % nonlinear D12 element of matrix
    DNl22 = 0; % nonlinear D22 element of matrix
    DNlMinR = 0; % mininum value for Dr matrix element
    DNlMin11 = 0; % mininum value for D11 matrix element
    DNlMin12 = 0; % mininum value for D12 matrix element
    DNlMin22 = 0; % mininum value for D22 matrix element
  end
  
  methods
    % Constructor
    % rhoMax: density diffusion goes to zero
    % posRotFlag [pos, rot]: flag for include position/rotation
    % n1,2,3: number of grid points in 1, 2, 3
    % d0: diffusion constant (perp, rotational )
    % ik: sqrt(-1) * k3 vector (cell)
    function obj = DensityDepAnisoDiffClass( type, rhoMax, posRotFlag, d0, ik, ...
        b, n1, n2, n3, phi)
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
        obj.Type = type;
        obj.setType();
        obj.N1 = n1;
        obj.N2 = n2;
        obj.N3 = n3;
        % convert rhoMax input, in bc, to rho = c / 2pi
        obj.K03 = n3 / 2 + 1;
        obj.RhoMax = rhoMax / b;
        obj.DnlNVec = [n1, n2];
        % store constant diffusion coefficients
        obj.D0Perp = d0(1);
        obj.D0R = d0(2);
        % set pre factors based on flag
        if obj.FlagPos
          obj.DimInclude = 1:2;
          [~,obj.DNlFactPerp] = obj.calcNlCoeff( obj.D0Perp );
          % build trig
          phi = reshape( phi, [1 1 n3] );
          sinPhi = sin( phi );
          cosPhi = cos( phi );
          obj.DNlFact11 = obj.buildDiffMatElement( ...
            obj.DNlFactPerp, sinPhi .* sinPhi );
          obj.DNlFact12 = obj.buildDiffMatElement( ...
            obj.DNlFactPerp, -sinPhi .* cosPhi );
          obj.DNlFact22 = obj.buildDiffMatElement( ...
            obj.DNlFactPerp, cosPhi .* cosPhi );
          obj.DNlMin11 = ...
            repmat( -obj.D0Perp * sinPhi .* sinPhi, [n1 n2 1] );
          obj.DNlMin12 = ...
            repmat( obj.D0Perp * sinPhi .* cosPhi, [n1 n2 1] );
          obj.DNlMin22 = ...
            repmat( -obj.D0Perp * cosPhi .* cosPhi, [n1 n2 1] );
          obj.PosTerms = 1:length(obj.DNlFact11);
          obj.DNl11 = zeros(obj.DnlNVec);
          obj.DNl12 = zeros(obj.DnlNVec);
          obj.DNl22 = zeros(obj.DnlNVec);
        end
        if obj.FlagRot
          obj.DimInclude = unique( [obj.DimInclude 3] );
          [~,obj.DNlFactR] = obj.calcNlCoeff( obj.D0R );
          obj.RotTerms = 1:length(obj.DNlFactR);
          obj.DNlMinR = -obj.D0R;
          obj.DNlR = zeros(obj.DnlNVec);
        end
        % handle k vec
        nVec = { [n1 1 1], [1 n2 1], [1 1 n3] };
        for ii = obj.DimInclude
          % scale rho by average excluded volume and angle
          obj.Ik{ii} = reshape( ik{ii}, nVec{ii} );
        end
      end
    end
    
    % Set the nl diffusion type
    function [obj] = setType( obj )
      if strcmp( obj.Type, 'lin' ) % linear in density dnl~c
        obj.TypeId = 1;
      elseif strcmp( obj.Type, 'quad' ) % quad in density dnl~c+c^2
        obj.TypeId = 2;
      elseif strcmp( obj.Type, 'tanh' ) % tanh in density dnl~tanh(c)
        obj.TypeId = 3;
      elseif strcmp( obj.Type, 'exp' ) % exp in density dnl~tanh(c)
        obj.TypeId = 4;
      elseif strcmp( obj.Type, 'rods' ) % thin rod form in density dnl~1/(1+c^2)
        obj.TypeId = 5;
      else
        obj.TypeId = 1;
        obj.Type = 'lin';
        fprintf('Incorrect density type. Setting to linear\n')
      end
      fprintf('Id: %d; Type: %s\n', obj.TypeId, obj.Type )
    end % setType
    
    % calc the nl diffusion coefficient factor based on type
    function [obj, dNlFact] = calcNlCoeff( obj, d0 )
      % linear
      if obj.TypeId == 1
        dNlFact{1} = -d0 / obj.RhoMax ; % linear coeff
        % quad
      elseif obj.TypeId == 2
        dNlFact{1} = -2 * d0 / obj.RhoMax ; % linear coeff
        dNlFact{2} = d0 / ( obj.RhoMax .^ 2 ); % quad coeff
        % tanh/exp/rods
      elseif obj.TypeId == 3 || obj.TypeId == 4 || obj.TypeId == 5
        dNlFact{1} = d0; % amplitude
        dNlFact{2} = obj.RhoMax; % scale
      end
    end % calcNlCoeff
    
    % build diffusion matrix elements and cells for each type
    function dNlFact = buildDiffMatElement( obj, dNl, trigFunction )
      dNlFact = cell(1,2);
      if obj.TypeId == 1 || obj.TypeId == 2
        for ii = 1:length(dNl)
          dNlFact{ii} = dNl{ii} * trigFunction;
          dNlFact{ii} = repmat( dNlFact{ii}, [obj.N1, obj.N2, 1] );
        end
      else
        dNlFact{1} = dNl{1} * trigFunction;
        dNlFact{2} = dNl{2};
      end
    end % buildDiffMatElement
    
    % Set the nl diffusion coeff
    function [obj] = calcDiffNl( obj, rho )
      % calc Dnl by type
      if obj.TypeId == 1 || obj.TypeId == 2
        if obj.FlagRot
          obj.DNlR = obj.calcDnlPower( obj.DnlNVec, obj.DNlFactR, ...
            obj.RotTerms, rho, obj.RhoMax, obj.DNlMinR );
        end
        if obj.FlagPos
          obj.DNl11 = obj.calcDnlPower( obj.DnlNVec, obj.DNlFact11, ...
            obj.PosTerms, rho, obj.RhoMax, obj.DNlMin11 );
          obj.DNl12 = obj.calcDnlPower( obj.DnlNVec, obj.DNlFact12, ...
            obj.PosTerms, rho, obj.RhoMax, obj.DNlMin11 );
          obj.DNl22 = obj.calcDnlPower( obj.DnlNVec, obj.DNlFact22, ...
            obj.PosTerms, rho, obj.RhoMax, obj.DNlMin22 );
        end
      elseif obj.TypeId == 3
        if obj.FlagRot
          obj.DNlR = obj.calcDnlTanh( obj.DNlFactR, rho );
        end
        if obj.FlagPos
          obj.DNl11 = obj.calcDnlTanh( obj.DNlFact11, rho );
          obj.DNl12 = obj.calcDnlTanh( obj.DNlFact12, rho );
          obj.DNl22 = obj.calcDnlTanh( obj.DNlFact22, rho );
        end
      elseif obj.TypeId == 4
        if obj.FlagRot
          obj.DNlR = obj.calcDnlExp( obj.DNlFactR, rho );
        end
        if obj.FlagPos
          obj.DNl11 = obj.calcDnlExp( obj.DNlFact11, rho );
          obj.DNl12 = obj.calcDnlExp( obj.DNlFact12, rho );
          obj.DNl22 = obj.calcDnlExp( obj.DNlFact22, rho );
        end
      elseif obj.TypeId == 5
        if obj.FlagRot
          obj.DNlR = obj.calcDnlRods( obj.DNlFactR, rho );
        end
        if obj.FlagPos
          obj.DNl11 = obj.calcDnlRods( obj.DNlFact11, rho );
          obj.DNl12 = obj.calcDnlRods( obj.DNlFact12, rho );
          obj.DNl22 = obj.calcDnlRods( obj.DNlFact22, rho );
        end
      end % end types
    end % calcDiffNl
    
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
          obj.DNl11 .* iotaTemp{1} + obj.DNl12 .* iotaTemp{2} ) );
        iotaFtTemp{2} = fftshift( fftn( ...
          obj.DNl12 .* iotaTemp{1} + obj.DNl22 .* iotaTemp{2} ) );
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
    end
  end %methods
  
  methods (Static)
   
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

    % calcDnl with rods function.
    function dnl = calcDnlRods( c, conc )
      dnl = c{1} .* (  1 ./ ( 1  +  conc / c{2} ) .^2 - 1 );
    end

    % calcDnl with exp function.
    function dnl = calcDnlExp( c, conc )
      dnl = c{1} .* (  exp( -conc / c{2} ) - 1 );
    end

    % calcDnl with tanh function.
    function dnl = calcDnlTanh( c, conc )
      dnl = -c{1} .* tanh( conc / c{2} );
    end

    % calcDnl with power series function.
    function dnl = calcDnlPower( nlNVec, c, numTerm, conc, conc_max, dnlMin )
        dnl = zeros( nlNVec );
        % only calculate it for density less than max
        inds2calc = conc < conc_max;
        % calc positional Dnl
        for ii = 1:numTerm
          dnl(inds2calc) = dnl(inds2calc) +...
            c{ii} .* ( conc(inds2calc) .^ ii );
        end
        % for densities about density max, set to min value
        dnl( ~inds2calc ) = dnlMin;
    end
   end % static methods
end %class
