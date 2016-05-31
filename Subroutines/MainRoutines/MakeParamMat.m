% [paramMat] = MakeParamMat( ParamObj, RhoInit, Flags )
% 
% Description: Creates a parameter matrix that is used by RunHardRod

function [paramMat, numRuns] = MakeParamMat( ParamObj, RhoInit, Flags )
%Find number parameters
numbc = length(ParamObj.bc);
numvD = length(ParamObj.vD);
numIC = length(RhoInit.IntCond);
numSM = length(Flags.StepMeth);

% Handls all grid points equal and square box different
if Flags.AllNsSame == 1
  Nvec = unique( [ParamObj.Nx ParamObj.Ny ParamObj.Nm] );
  numN = length(Nvec);
  numNx = numN;
  numNy = 1;
  numNm = 1;
else
  numNx = length(ParamObj.Nx);
  numNy = length(ParamObj.Ny);
  numNm = length(ParamObj.Nm);
  numN = numNx * numNy * numNm;
end

if Flags.SquareBox == 1
  Lvec = unique( [ParamMaster.Lx ParamMaster.Ly] );
  numL = length(Lvec);
  numLx = numL;
  numLy = 1;
else
  numLx = length(ParamObj.Lx);
  numLy = length(ParamObj.Ly);
  numL = numLx * numLy;
end

% number of parameters
numNxNy = numNx*numNy;
numNxNyNm = numNxNy * numNm;
numNxNyNmLx = numNxNyNm * numLx;
numNxNyNmLxLy = numNxNyNmLx * numLy;
numNxNyNmLxLyvD = numNxNyNmLxLy * numvD;
numParams = numNxNyNmLxLyvD * numbc;
numParamsIC = numParams * numIC;
numRuns = numParamsIC * numSM;

% Create Paramater matrix
% paramMat columns: (Nx, Ny, Nm, Lx, Ly, vD, bc, IC, SM, runID)
paramMat = zeros( numRuns, 10);
for i = 1:numNx
  for j = 1:numNy
    for k = 1:numNm
      for l = 1:numLx
        for m = 1:numLy
          for n = 1:numvD
            for o = 1:numbc
              for p = 1:numIC
                for q = 1:numSM
                  rowInd =  1 + (i-1) + (j-1) * numNx + ( k-1) * numNxNy + ...
                    (l-1) * numNxNyNm + (m-1) * numNxNyNmLx + ...
                    (n-1) * numNxNyNmLxLy + (o-1) * numNxNyNmLxLyvD + ...
                    (p-1) * numParams + (q-1) * numParamsIC;
                  
                  % Special cases
                  if Flags.AllNsSame
                    paramMat(rowInd,1:3) = Nvec(i);
                  else
                    paramMat(rowInd,1) = ParamObj.Nx(i);
                    paramMat(rowInd,2) = ParamObj.Ny(j);
                    paramMat(rowInd,3) = ParamObj.Nm(k);
                  end
                  
                  if Flags.SquareBox
                    paramMat(rowInd,4:5) = Lvec(l);
                  else
                    paramMat(rowInd,4) = ParamObj.Lx(l);
                    paramMat(rowInd,5) = ParamObj.Ly(m);
                  end
                  % Everything else
                  paramMat(rowInd,6) = ParamObj.vD(n);
                  paramMat(rowInd,7) = ParamObj.bc(o);
                  paramMat(rowInd,8) = RhoInit.IntCond(p);
                  paramMat(rowInd,9) = Flags.StepMeth(q);
                  paramMat(rowInd,10) = ParamObj.runID - 1 + rowInd;
                end
              end
            end
          end
        end
      end
    end
  end
end


