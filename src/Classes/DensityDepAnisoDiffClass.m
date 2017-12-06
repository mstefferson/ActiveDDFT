classdef DensityDepAnisoDiffClass < handle
  properties
    Flag = []; % flag to calculate or not
    N1 = []; % grid points in 1
    N2 = []; % grid points in 2
    N3 = []; % grid points in 3
    RhoMax = []; % density, c, where rho goes to zero
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
  end
  
  methods
    % Constructor
    % rhoMax: density diffusion goes to zero
    % n1,2,3: number of grid points in 1, 2, 3
    % d0: diffusion constant (perp, rotational )
    % ik: sqrt(-1) * k3 vector (cell)
    function obj = DensityDepAnisoDiffClass( rhoMax, d0, ik, ...
        b, n1, n2, n3, phi)
      if rhoMax == 0
        obj.Flag = 0;
      else
        obj.Flag = 1;
        obj.N1 = n1;
        obj.N2 = n2;
        obj.N3 = n3;
        % build trig
        phi = reshape( phi, [1 1 n3] );
        sinPhi = sin( phi );
        cosPhi = sin( phi );
        % convert rhoMax input, in bc, to rho = c / 2pi
        obj.RhoMax = rhoMax / b / (2*pi);
        obj.D0Perp = d0(1);
        obj.D0R = d0(2);
        obj.DNlFactPerp = -obj.D0Perp / obj.RhoMax ;
        obj.DNlFactR = -obj.D0R / obj.RhoMax ;
        obj.DNlFact11 = obj.DNlFactPerp * sinPhi .* sinPhi;
        obj.DNlFact12 = -obj.DNlFactPerp * sinPhi .* cosPhi;
        obj.DNlFact22 = obj.DNlFactPerp * cosPhi .* cosPhi;
        % handle k vec
        nVec = { [n1 1 1], [1 n2 1], [1 1 n3] };
        for ii = 1:3
          % scale rho by average excluded volume and angle
          obj.Ik{ii} = reshape( ik{ii}, nVec{ii} );
        end
      end
    end
    
    % Set the nl diffusion coeff
    function [obj] = calcDiffNl( obj, rho )
      obj.DNlR = obj.DNlFactR .* rho;
      obj.DNl11 = obj.DNlFact11 .* rho;
      obj.DNl12 = obj.DNlFact12 .* rho; 
      obj.DNl22 = obj.DNlFact22 .* rho; 
      obj.DNlR( obj.DNlR < -obj.D0R ) = -obj.D0R;
      % I am unsure about this
      obj.DNl11( obj.DNl11 < -obj.D0Perp ) = -obj.D0Perp;
      obj.DNl12( obj.DNl12 < -obj.D0Perp ) = -obj.D0Perp;
      obj.DNl22( obj.DNl22 < -obj.D0Perp ) = -obj.D0Perp;
    end
    
    % calc d rho
    function [dRho_dt] = calcDrho( obj, rho, rhoFt, jEx )
      obj.calcDiffNl( rho );
      % "flux" without mobility
      jTemp = cell(1, 3 );
      for ii = 1:3
        jDiffTemp  = -real( ifftn( ifftshift( obj.Ik{ii} .* rhoFt ) ) );
        jTemp{ii} = jDiffTemp + jEx{ii};
      end
      % fourier transform of true flux
      jftTemp1 = fftshift( fftn( ...
       obj.DNl11 .* jTemp{1} + obj.DNl12 .* jTemp{2} ) );
      jftTemp2 = fftshift( fftn( ...
       obj.DNl12 .* jTemp{1} + obj.DNl22 .* jTemp{2} ) );
      jftTemp3 = fftshift( fftn( ...
       obj.DNlR .* jTemp{3} ) );      
      % minus divergence of flux
      dRho_dt = ...
        -obj.Ik{1} .* ( jftTemp1 ) + ...
        -obj.Ik{2} .* ( jftTemp2 ) + ...
        -obj.Ik{3} .* ( jftTemp3 );
   end
  end %methods
end %class

