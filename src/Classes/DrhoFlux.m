classdef DrhoFlux < handle
  properties
    % flags
    FlagJ1= [];
    FlagJ2= [];
    FlagJ3= [];
    AnisoFlag= [];
    N1 = 0;
    N2 = 0;
    N3 = 0;
    % mobility matrix elements
    Mob11 = 0;
    Mob12 = 0;
    Mob22 = 0;
    Mob33 = 0;
    % k vectors. must be correct shape (from mobility obj)
    Ik1 = 0; 
    Ik2 = 0;
    Ik3 = 0;
  end
  
  methods
    % Constructor
    function obj = DrhoFlux( diffObj, polarDriveFlag, anisoFlag, ...
        dv1Flag, dv2Flag, dv3Flag, n1, n2, n3, ik1, ik2, ik3 )
      % set properties
      % handle flags
      if polarDriveFlag || dv1Flag
        obj.FlagJ1 = 1;
        obj.Ik1 = ik1;
      else
        obj.FlagJ1 = 0;
        obj.Ik1 = 0;
      end
      if polarDriveFlag || dv2Flag
        obj.FlagJ2 = 1;
        obj.Ik2 = ik2;
      else
        obj.FlagJ2 = 0;
        obj.Ik2 = 0;
      end
      if dv3Flag
        obj.FlagJ3 = 1;
        obj.Ik3 = ik3;
      else
        obj.FlagJ3 = 0;
        obj.Ik3 = 0;
      end
      obj.N1 = n1;
      obj.N2 = n2;
      obj.N3 = n3;
      % Anisoflag
      obj.AnisoFlag = anisoFlag;
      % store mobility
      obj.Mob11 = diffObj.Mob11;
      obj.Mob12 = diffObj.Mob12;
      obj.Mob22 = diffObj.Mob22;
      obj.Mob33 = diffObj.Mob33;
    end % constructor
    % calc preflux iota ( j = D iota )
    function [negDivFluxEx_FT] = calcDrho( obj, iota1, iota2, iota3 )
      % initialize
      negDivFluxEx_FT = zeros( obj.N1, obj.N2, obj.N3 );
      % calc j1
      if obj.FlagJ1
        % handle mixing if there is any
        if obj.AnisoFlag == 0
          j1 = obj.Mob11 .* iota1 + obj.Mob12 .* iota2;
        else
          j1 = obj.Mob11 .* iota1;
        end
        j1Ft = fftshift( fftn( j1 ) );
        % calculate
        negDivFluxEx_FT = - ( obj.Ik1 .* j1Ft );
      end
      % calc j2
      if obj.FlagJ2
        % handle mixing if there is any
        if obj.AnisoFlag == 0
          j2 = obj.Mob12 .* iota1 + obj.Mob22 .* iota2;
        else
          j2 = obj.Mob22 .* iota2;
        end
        j2Ft = fftshift( fftn( j2 ) );
        % calculate
        negDivFluxEx_FT = negDivFluxEx_FT - ( obj.Ik2 .* j2Ft );
      end
      % calc j3
      if obj.FlagJ3
        j3 = obj.Mob33 .* iota3;
        j3Ft = fftshift( fftn( j3 ) );
        % calculate
        negDivFluxEx_FT = negDivFluxEx_FT - ( obj.Ik3 .* j3Ft );
      end
    end    % calc dRho
  end % methods
end %class

