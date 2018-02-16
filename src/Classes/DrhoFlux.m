classdef PolarDriveClass < handle
  properties
    FlagJ1= [];
    FlagJ2= [];
    FlagJ3= [];
    AnisoFlag= [];
    N1 = 0;
    N2 = 0;
    N3 = 0;
    Mob11 = 0; % flux 1 (x) before mobility
    Mob12 = 0;  % flux 2 (y) before mobility  
    Mob22 = 0; %factor to hit iota1
  end
  
  methods
    % Constructor
    function obj = DrhoFlux( diffObj, polarDriveFlag, anisoFlag ...
     dv1Flag, dv2Flag, dv3Flag, n1, n2, n3 )
      % set properties
      % handle flags
      if polarDrive.Flag || interObj.dv1Flag
        FlagJ1 = 1;
      end
      if polarDrive.Flag || interObj.dv2Flag
        FlagJ2 = 1;
      end
      if interObj.dv3Flag
        FlagJ3 = 1;
      end
      obj.N1 = n1;
      obj.N2 = n2;
      obj.N3 = n3;
      % Anisoflag
      AnisoFlag = anisoFlag;
      % store mobility
      Mob11 = diffObj.mob11;
      Mob12 = diffObj.mob12;
      Mob22 = diffObj.mob22;
      Mob33 = diffObj.mob33;
    end % constructor
    % calc preflux iota ( j = D iota )
    function [negDivFluxEx_FT] = calcDrho( obj, iota1, iota2, iota3 )
      % initialize
      negDivFluxEx_FT = zeros( obj.N1, obj.N2, obj.N3 );
      % calc j1
      if obj.FlagJ1
        % handle mixing if there is any
        if obj.AnisoFlag == 0
          j1 = diffObj.mob11 .* iota1 + diffObj.mob12 .* iota2;
        else
          j1 = diffObj.mob11 .* iota1;
        end
        j1Ft = fftshift( fftn( j1 ) );
        % calculate 
        negDivFluxEx_FT = - ( diffObj.ik1 .* j1Ft );
      end
      % calc j2
      if obj.FlagJ2
        % handle mixing if there is any
        if obj.AnisoFlag == 0
          j2 = diffObj.mob12 .* iota1 + diffObj.mob22 .* iota2;
        else
          j2 = diffObj.mob22 .* iota2;
        end
        j2Ft = fftshift( fftn( j2 ) );
        % calculate 
        negDivFluxEx_FT = negDivFluxEx_FT - ( diffObj.ik2 .* j2Ft );
      end
      % calc j3
      if obj.FlagJ3
        j3 = diffObj.mob33 .* iota3;
        j3Ft = fftshift( fftn( j3 ) );
        % calculate 
        negDivFluxEx_FT = negDivFluxEx_FT - ( diffObj.ik3 .* j3Ft );
      end
    end    % calc dRho
  end % methods
end %class

