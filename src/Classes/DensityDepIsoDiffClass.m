classdef DensityDepIsoDiffClass < handle
  properties
    Flag = []; % flag to calculate or not
    N1 = []; % grid points in 1
    N2 = []; % grid points in 2
    N3 = []; % grid points in 3
    D0 = zeros(1,3); % constant diffusion coefficient
    RhoMax = []; % density, c, where rho goes to zero
    DNlFact = cell(1,3); % non linear diffusion coefficient factor
    DNl = cell(1,3); % non linear diffusion contribution
    Ik = cell(1,3); % sqrt(-1) * k1 vec
    NlDiffComponents = []; % compentes (1,2,3) we want nl diffusion in
    DNlMin = [];
  end
  
  methods
    % Constructor
    % rhoMax: density diffusion goes to zero
    % n1,2,3: number of grid points in 1, 2, 3
    % d0: diffusion constant
    % ik1: sqrt(-1) * k1 vector
    % ik2: sqrt(-1) * k3 vector
    function obj = DensityDepIsoDiffClass( rhoMax, indsWant, d0, ik, ...
        b, n1, n2, n3)
      if rhoMax == 0
        obj.Flag = 0;
      else
        obj.Flag = 1;
        obj.N1 = n1;
        obj.N2 = n2;
        obj.N3 = n3;
        % convert rhoMax input, in bc, to rho = c / 2pi
        obj.RhoMax = rhoMax / b / (2*pi);
        obj.NlDiffComponents = indsWant;
        obj.D0 = d0;
        obj.DNlFact = -obj.D0 / obj.RhoMax ;
        obj.DNlMin = -obj.D0;
        nVec = { [n1 1 1], [1 n2 1], [1 1 n3] };
        for ii = obj.NlDiffComponents
          % scale rho by average excluded volume and angle
          obj.Ik{ii} = reshape( ik{ii}, nVec{ii} );
          obj.DNl{ii} = zeros(n1,n2,n3);
        end
      end
    end
    
    % Set the nl diffusion coeff
    function [obj] = calcDiffNl( obj, rho )
      for ii = obj.NlDiffComponents
        obj.DNl{ii} = obj.DNlFact(ii) .* rho;
        obj.DNl{ii} = obj.fixNegativeDiff( obj.DNl{ii}, obj.DNlMin(ii) );
        %obj.DNl{ii}( obj.DNl{ii} <  -obj.D0(ii) ) = -obj.D0(ii);
      end
    end
    
    % fix densities that are too low
    function dCurrent = fixNegativeDiff( obj, dCurrent, dMin )
      logInds = dCurrent < dMin;
      if isscalar( dMin  )
        dCurrent( logInds ) = dMin;
      else
        dCurrent( logInds ) = dMin( logInds );
      end
    end

    % calc d rho
    function [dRho_dt] = calcDrho( obj, rho, rhoFt, jEx )
      obj.calcDiffNl( rho );
      % "flux" without mobility
      dRho_dt = zeros( obj.N1, obj.N2, obj.N3 );
      for ii = obj.NlDiffComponents
        jDiffTemp  = -real( ifftn( ifftshift( obj.Ik{ii} .* rhoFt ) ) );
        jftTemp = fftshift( fftn( ...
          obj.DNl{ii} .* ( jDiffTemp + jEx{ii} ) ) );
        dRho_dt = dRho_dt - obj.Ik{ii} .* jftTemp;
      end
   end
  end %methods
end %class

