% DensityDepIsoDiffClass: handles non-linear density contributions to diffusion
% considering an anisotropic diffusion coefficient
% Can turn on NL diffusion in x,y,phi

classdef DensityDepIsoDiffClass < handle
  properties
    Flag = []; % flag to calculate or not
    FlagPos = []; % flag to calculate position contribution or not
    FlagRot = []; % flag to calculate rotational contribution or not
    Type = []; % 'lin' or 'quad'
    TypeId = []; % 1 = linear, 2 = quad 3 = tanh
    PosTerms =[]; % number of coefficients in position terms
    RotTerms = []; % number of coefficients in position terms
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
        obj.Type = order;
        obj.setType();
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
          obj.PosTerms = 1:length(obj.DNlFactPos);
          obj.DNlMinPos = -obj.D0Pos;
          obj.DNlPos = zeros( obj.DnlNVec );
        end
        if obj.FlagRot
          obj.DimInclude = unique( [obj.DimInclude 3] );
          [~,obj.DNlFactR] = obj.calcNlCoeff( obj.D0R );
          obj.RotTerms = 1:length(obj.DNlFactR);
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
    function [obj] = setType( obj )
      if strcmp( obj.Type, 'lin' )
        obj.TypeId = 1;
      elseif strcmp( obj.Type, 'quad' )
        obj.TypeId = 2;
      elseif strcmp( obj.Type, 'tanh' )
        obj.TypeId = 3;
      else
        obj.TypeId = 1;
        obj.Type = 'lin';
        fprintf('Incorrect density order. Setting to linear\n')
      end
      fprintf('Id: %d; Type: %s\n', obj.TypeId, obj.Type )
    end % setType
    
    % calc the nl diffusion coefficient factor based on order
    function [obj, dNlFact] = calcNlCoeff( obj, d0 )
      % linear
      if obj.TypeId == 1
        dNlFact{1} = -d0 / obj.RhoMax ;
      % quad
      elseif obj.TypeId == 2
        dNlFact{1} = -2 * d0 / obj.RhoMax ;
        dNlFact{2} = d0 / ( obj.RhoMax .^ 2 );
      % tanh
      elseif obj.TypeId == 3
        dNlFact{1} = -d0;
        dNlFact{2} = obj.RhoMax;
      end
    end % calcNlCoeff
    
    % Set the nl diffusion coeff
    function [obj] = calcDiffNl( obj, rho )
      obj.DNlR = zeros( obj.DnlNVec );
      obj.DNlPos = zeros( obj.DnlNVec );
      if obj.TypeId == 1 || obj.TypeId == 2
        inds2calc = rho < obj.RhoMax;
      else
        try
          inds2calc = logical( rho );
        catch err
          keyboard
        end
      end
      % calc rotational Dnl
      if obj.TypeId == 1 || obj.TypeId == 2
        for ii = 1:obj.RotTerms
          obj.DNlR(inds2calc) = obj.DNlR(inds2calc) + ...
            obj.DNlFactR{ii} .* ( rho(inds2calc) .^ ii );
        end
        % calc positional Dnl
        for ii = 1:obj.PosTerms
          obj.DNlPos(inds2calc) = obj.DNlPos(inds2calc) +...
            obj.DNlFactPos{ii} .* ( rho(inds2calc) .^ ii );
        end
      elseif obj.TypeId == 3
        if obj.FlagRot
          obj.DNlR(inds2calc) = obj.DNlR(inds2calc) + ...
            obj.DNlFactR{1} .* tanh( rho(inds2calc) / obj.DNlFactR{2} );
        end
        if obj.FlagPos
          obj.DNlPos(inds2calc) = obj.DNlPos(inds2calc) +...
            obj.DNlFactPos{1} .* tanh( rho(inds2calc) / obj.DNlFactPos{2} );
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
      c = 2*pi / obj.N3 * real( ifftn(ifftshift( rhoFt(:,:,obj.K03) ) ) );
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
      %dRho_dt_Ft = cell(1,2);
      %dRho_dt_real = cell(1,2);
      for ii = obj.DimInclude
        dRho_dt = dRho_dt - obj.Ik{ii} .* ( iotaFtTemp{ii} );
       % dRho_dt_Ft{ii} = - obj.Ik{ii} .* ( iotaFtTemp{ii} );
       % dRho_dt_real{ii} = real( ifftn( ifftshift( dRho_dt_Ft{ii} ) ) );
       % disp(ii)
      end
      %dRho_dt_sum = real( ifftn( ifftshift( dRho_dt ) ) );
      %% plot flux
      %figure()
      %subplot(2,2,1)
      %sum1 = sum( iotaTemp{1}, 3 );
      %imagesc( sum1 ); colorbar
      %subplot(2,2,2)
      %plot( sum1(:,obj.K03) )
      %hold
      %plot( sum1(obj.K03,:) )
      %hold
      %title( 'flux 1' )
      %subplot(2,2,3)
      %sum2 = sum( iotaTemp{2}, 3 );
      %imagesc( sum2 ); colorbar
      %subplot(2,2,4)
      %plot( sum2(:,obj.K03) )
      %hold
      %plot( sum2(obj.K03,:) )
      %hold
      %title( 'flux 2' )      
      %%% plot dRho
      %figure()
      %subplot(2,2,1)
      %sum1 = sum( dRho_dt_real{1}, 3 );
      %imagesc( sum1 ); colorbar
      %title( 'dRho1' )
      %subplot(2,2,2)
      %plot( sum1(:,obj.K03) )
      %hold
      %plot( sum1(obj.K03,:) )
      %hold
      %subplot(2,2,3)
      %sum2 = sum( dRho_dt_real{2}, 3 );
      %imagesc( sum2 ); colorbar
      %title( 'dRho1' )
      %subplot(2,2,4)
      %plot( sum2(:,obj.K03) )
      %hold
      %plot( sum2(obj.K03,:) )
      %hold
      %title( 'dRho2' )
      %%% Quick and dirty derivative
      %figure()
      %dIota1 = obj.N1 / 10 * reshape( ...
        %iotaTemp{1}(2:end,:,:) - iotaTemp{1}(1:end-1,:,:), [obj.N1-1 obj.N2 obj.N3] ) ;
      %dIota2 = obj.N2 / 10 * reshape( ...
        %iotaTemp{2}(:,2:end,:) - iotaTemp{2}(:,1:end-1,:), [obj.N1 obj.N2-1 obj.N3] );
      %subplot(2,2,1)
      %sum1 = sum( dIota1,3 );
      %imagesc( sum1 ); colorbar
      %title( 'dRho1 manual' )
      %subplot(2,2,2)
      %plot( sum1(:,obj.K03) )
      %hold
      %plot( sum1(obj.K03,:) )
      %hold
      %title( 'dRho1 manual' )
      %subplot(2,2,3)
      %sum2 = sum( dIota2,3 );
      %imagesc( sum2 ); colorbar
      %title( 'dRho2 manual' )
      %subplot(2,2,4)
      %plot( sum2(:,obj.K03) )
      %hold
      %plot( sum2(obj.K03,:) )
      %hold
      %title( 'dRho2 manual' )
      %%% Check in k-space
      %ind = 1;
      %figure()
      %% 1
      %temp1 = iotaTemp{1}(:,:,ind);
      %subplot(2,3,1)
      %imagesc( temp1 ); colorbar
      %temp1Ft = fftshift( fftn( temp1 ) );
      %dtemp1Ft = obj.Ik{1} .* temp1Ft;
      %dtemp1 = real( ifftn( ifftshift( dtemp1Ft ) ) );
      %subplot(2,3,2)
      %imagesc( dtemp1 ); colorbar
      %dtemp1_v2 = real( ifftn( ifftshift( obj.Ik{1} .* iotaFtTemp{1}(:,:,obj.K03) ) ) );
      %subplot(2,3,3)
      %imagesc( dtemp1_v2 ); colorbar 
      %% 2
      %temp2 = iotaTemp{2}(:,:,ind);
      %subplot(2,3,4)
      %imagesc( temp2 ); colorbar
      %temp2Ft = fftshift( fftn( temp2 ) );
      %dtemp2Ft = obj.Ik{2} .* temp2Ft;
      %dtemp2 = real( ifftn( ifftshift( dtemp2Ft ) ) );
      %subplot(2,3,5)
      %imagesc( dtemp2 ); colorbar  
      %dtemp2_v2 = real( ifftn( ifftshift( obj.Ik{2} .* iotaFtTemp{2}(:,:,obj.K03) ) ) );
      %subplot(2,3,6)
      %imagesc( dtemp2_v2 ); colorbar 
      %%% compare ft
      %figure()
      %subplot( 2,2,1 )
      %imagesc( imag( temp1Ft ) ); colorbar
      %title('temp1')
      %subplot( 2,2,2 )
      %imagesc( imag( iotaFtTemp{1}(:,:,obj.K03) ) ); colorbar
      %title('iotaFt(1)')
      %subplot( 2,2,3 )
      %imagesc( imag( temp2Ft ) ); colorbar
      %title('temp2')
      %subplot( 2,2,4 )
      %imagesc( imag( iotaFtTemp{2}(:,:,obj.K03) ) ); colorbar
      %title('iotaFt(2)')
      %%% ifft temp
      %%%
      %figure()
      %plot( reshape( imag( iotaFtTemp{1}(obj.K03+1,obj.K03+1,:) ), [1 obj.N3] ) )
      %%%
      %keyboard
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
