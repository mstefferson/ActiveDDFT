classdef PairDistIntegratorClass < handle
  properties
    CalcG1G2Flag = []; % Flag to calculate g1/g2
    N1 = []; % number of grid points in n1
    N2 = []; % number of grid points in n2
    N3 = []; % number of grid points in n3
    MayerInds1 = []; % Mayer inds
    MayerInds2 = []; % Mayer inds
    MayerInds3 = []; % Mayer inds
    AllInds1 = []; % All inds
    AllInds2 = []; % All inds
    AllInds3 = []; % All inds
    CosPhiMayInds = []; % cos for g1 calc
    CosSqrPhiMayInds = []; % cos^2 for g1 calc
    X1 = [];
    X2 = [];
    Phi = [];
    ExpV = [];
    FixedInd1 = [];
    FixedInd2 = [];
    ShiftInds1 = [];
    ShiftInds2 = [];
    ShiftInds3 = [];
    IntOverRprime0 = [];
    IntOverRprime1 = [];
    IntOverRprime2 = [];
    Delta0 = [];
    Delta1 = [];
    Delta2 = [];
  end
  
  methods
    % Constructor
    function obj = PairDistIntegratorClass( n1, n2, n3, calcG1G2Flag, ...
        x1, x2, phi, mayer )
      obj.N1 = n1;
      obj.N2 = n2;
      obj.N3 = n3;
      obj.CalcG1G2Flag =  calcG1G2Flag;
      obj.AllInds1 = 1:n1;
      obj.AllInds2 = 1:n2;
      obj.AllInds3 = 1:n3;
      obj.MayerInds1 = [0:n1/2 -n1/2+1:1:-1];
      obj.MayerInds2 = [0:n2/2 -n2/2+1:1:-1];
      obj.MayerInds3 = [0:n3/2 -n3/2+1:1:-1];
      obj.X1 = x1;
      obj.X2 = x2;
      obj.Phi = phi;
      obj.ExpV = mayer + 1;
      % initialize
      obj.IntOverRprime0= zeros( n3, n3 );
      if calcG1G2Flag
        % take the cos once
        obj.CosPhiMayInds = cos(mayerInds3);
        obj.CosSqrPhiMayInds = cos(mayerInds3).^2;
        % initialize
        obj.IntOverRprime1 = zeros( n3, n3 );
        obj.IntOverRprime2 = zeros( n3, n3 );
      else
        obj.IntOverRprime1 = 0;
        obj.IntOverRprime2 = 0;
      end
    end
    
    % update shifted inds for intergral
    function updateShiftInds( obj, ind, dim )
      % update based on dimension
      if dim == 1
        % set x
        obj.FixedInd1 = ind;
        r1Temp = obj.MayerInds1( ind );
        % set x+x'
        obj.ShiftInds1 = mod( (r1Temp-1) + obj.AllInds1 - 1, obj.N1 ) + 1;
      end
      if dim == 2
        % set y
        obj.FixedInd2 = ind;
        r2Temp = obj.MayerInds2( ind );
        % set y+y'
        obj.ShiftInds2 = mod( (r2Temp-1) + obj.AllInds2 - 1, obj.N2 ) + 1;
      end
    end % updateShiftInds
    % calculate integrals
    function [obj] = calcDeltaIntegrals( obj, rho )
      % calculate integrands
      obj.calcIntegrand(rho);
      obj.Delta0 = trapz_periodic( obj.Phi, ...
        trapz_periodic( obj.Phi, obj.IntOverRprime0, 1 ), 2 );
      if obj.CalcG1G2Flag
        obj.Delta1 = trapz_periodic( obj.Phi, ...
          trapz_periodic( obj.Phi, obj.IntOverRprime1, 1 ), 2 );
        obj.Delta2 = trapz_periodic( obj.Phi, ...
          trapz_periodic( obj.Phi, obj.IntOverRprime2, 1 ), 2 );
      end
    end % calc DeltaIntergral
    
    % calculate integrands
    function [obj] = calcIntegrand( obj, rho )
      % integrate of r' each u and u+u'.
      for mm = 1:obj.N3
        phi1Temp = obj.MayerInds3(mm);
        allInds3mm = obj.AllInds3(mm);
        shiftInds3 = mod( (phi1Temp-1) +  obj.AllInds3 - 1, obj.N3 ) + 1;
        shiftInds3mm = shiftInds3(mm);
        if obj.CalcG1G2Flag
          cosPhiTemp = obj.CosPhiMayInds(mm);
          cosSqrPhiTemp = obj.CosSqrPhiMayInds(mm);
        end
        for nn = 1:obj.N3
          allInds3nn = obj.AllInds3(nn);
          shiftInds3nn = shiftInds3(nn);
          %g0
          mat2Intdelta0 =  obj.ExpV( obj.FixedInd1, obj.FixedInd2, allInds3nn,  shiftInds3mm ) .* ...
            rho( obj.ShiftInds1, obj.ShiftInds2, shiftInds3nn ) .* rho(obj.AllInds1, obj.AllInds2, allInds3mm );
          obj.IntOverRprime0(nn,mm) = trapz_periodic( obj.X1, ...
            trapz_periodic( obj.X2, mat2Intdelta0, 2 ), 1);
          if obj.CalcG1G2Flag
            %g1
            mat2Intdelta1 =  cosPhiTemp * mat2Intdelta0;
            obj.IntOverRprime1(nn,mm)= trapz_periodic( obj.X1, ...
              trapz_periodic( obj.X2, mat2Intdelta1, 2 ), 1);
            %g2
            mat2Intdelta2 =  cosSqrPhiTemp * mat2Intdelta0;
            obj.IntOverRprime2(nn,mm) = trapz_periodic( obj.X1, ...
              trapz_periodic( obj.X2, mat2Intdelta2, 2 ), 1);
          end
        end % nn
      end % mm
    end % calcInts
  end %methods
end %class
