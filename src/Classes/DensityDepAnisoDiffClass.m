classdef DensityDepAnisoDiffClass < handle
  properties
    Flag = []; % flag to calculate or not
    FlagPos = []; % flag to calculate or not
    FlagRot = []; % flag to calculate or not
    N1 = []; % grid points in 1
    N2 = []; % grid points in 2
    N3 = []; % grid points in 3
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
    DNlRMin = 0; % mininum value for D11 matrix
    DNl11Min = 0; % mininum value for D11 matrix
    DNl12Min = 0; % mininum value for D12 matrix
    DNl22Min = 0; % mininum value for D22 matrix
  end
  
  methods
    % Constructor
    % rhoMax: density diffusion goes to zero
    % posRotFlag [pos, rot]: flag for include position/rotation
    % n1,2,3: number of grid points in 1, 2, 3
    % d0: diffusion constant (perp, rotational )
    % ik: sqrt(-1) * k3 vector (cell)
    function obj = DensityDepAnisoDiffClass( rhoMax, posRotFlag, d0, ik, ...
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
        obj.N1 = n1;
        obj.N2 = n2;
        obj.N3 = n3;
        % convert rhoMax input, in bc, to rho = c / 2pi
        obj.RhoMax = rhoMax / b / (2*pi);
        obj.D0Perp = d0(1);
        obj.D0R = d0(2);
        % set pre factors based on flag
        if obj.FlagPos
          obj.DimInclude = [ 1:2];
          % build trig
          phi = reshape( phi, [1 1 n3] );
          sinPhi = sin( phi );
          cosPhi = sin( phi );
          obj.DNlFactPerp = -obj.D0Perp / obj.RhoMax ;
          obj.DNlFact11 = obj.DNlFactPerp * sinPhi .* sinPhi;
          obj.DNlFact12 = -obj.DNlFactPerp * sinPhi .* cosPhi;
          obj.DNlFact22 = obj.DNlFactPerp * cosPhi .* cosPhi;
          obj.DNl11Min = repmat( -obj.D0Perp * sinPhi .* sinPhi, [n1 n2 1] );
          obj.DNl12Min = repmat( obj.D0Perp * sinPhi .* cosPhi, [n1 n2 1] );
          obj.DNl22Min = repmat( -obj.D0Perp * cosPhi .* cosPhi, [n1 n2 1] );
          obj.DNl11 = zeros(n1,n2,n3);
          obj.DNl12 = zeros(n1,n2,n3);
          obj.DNl22 = zeros(n1,n2,n3);
        else
          obj.DimInclude = [];
          obj.DNlFactPerp = 0;
          obj.DNlFact11 = 0;
          obj.DNlFact12 = 0;
          obj.DNlFact22 = 0;
        end
        if obj.FlagRot
          obj.DimInclude = 3;
          obj.DimInclude = [obj.DimInclude 3];
          obj.DNlFactR = -obj.D0R / obj.RhoMax ;
          obj.DNlRMin = -obj.D0R;
          obj.DNlR = zeros(n1,n2,n3);
        else
          obj.DNlFactR = 0;
        end
        % handle k vec
        nVec = { [n1 1 1], [1 n2 1], [1 1 n3] };
        for ii = obj.DimInclude
          % scale rho by average excluded volume and angle
          obj.Ik{ii} = reshape( ik{ii}, nVec{ii} );
        end
      end
    end
    
    % Set the nl diffusion coeff
    function [obj] = calcDiffNl( obj, rho )
      if obj.FlagRot
        obj.DNlR = obj.DNlFactR .* rho;
        %obj.DNlR( obj.DNlR < -obj.D0R ) = -obj.D0R;
        obj.DNlR = obj.fixNegativeDiff( obj.DNlR, obj.DNlRMin );
      end
      if obj.FlagPos
        obj.DNl11 = obj.DNlFact11 .* rho;
        obj.DNl12 = obj.DNlFact12 .* rho;
        obj.DNl22 = obj.DNlFact22 .* rho;
        % I am unsure about this
        %obj.DNl11( obj.DNl11 < -obj.D0Perp ) = -obj.D0Perp;
        %obj.DNl12( obj.DNl12 < -obj.D0Perp ) = -obj.D0Perp;
        %obj.DNl22( obj.DNl22 < -obj.D0Perp ) = -obj.D0Perp;
        obj.DNl11 = obj.fixNegativeDiff( obj.DNl11, obj.DNl11Min );
        obj.DNl12 = obj.fixNegativeDiff( obj.DNl12, obj.DNl12Min );
        obj.DNl22 = obj.fixNegativeDiff( obj.DNl22, obj.DNl22Min );
      end
    end
    
    % fix densities that are too low
    function dCurrent = fixNegativeDiff( obj, dCurrent, dMin )
      logInds = dCurrent < dMin;
      if isscalar( dMin )
        dCurrent( logInds ) = dMin;
      else
        dCurrent( logInds ) = dMin( logInds );
      end
    end
    
    % calc d rho
    function [dRho_dt] = calcDrho( obj, rho, rhoFt, jEx )
      obj.calcDiffNl( rho );
      % "flux" without mobility
      jTemp = cell( 1, 3 );
      jftTemp = cell( 1, 3 );
      for ii = 1:3
        if any( ii == obj.DimInclude )
          jDiffTemp  = -real( ifftn( ifftshift( obj.Ik{ii} .* rhoFt ) ) );
          jTemp{ii} = jDiffTemp + jEx{ii};
        else
          jTemp{ii} = 0;
        end
      end
      % fourier transform of true flux
      if obj.FlagPos
        jftTemp{1} = fftshift( fftn( ...
          obj.DNl11 .* jTemp{1} + obj.DNl12 .* jTemp{2} ) );
        jftTemp{2} = fftshift( fftn( ...
          obj.DNl12 .* jTemp{1} + obj.DNl22 .* jTemp{2} ) );
      end
      if obj.FlagRot
        jftTemp{3} = fftshift( fftn( ...
          obj.DNlR .* jTemp{3} ) );
      end
      % minus divergence of flux
      dRho_dt = 0;
      for ii = obj.DimInclude
        dRho_dt = dRho_dt - obj.Ik{ii} .* ( jftTemp{ii} );
      end
    end
  end %methods
end %class

