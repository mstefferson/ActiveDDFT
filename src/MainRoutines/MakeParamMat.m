% [paramMat] = MakeParamMat( ParamObj, rhoInit, flags )
%
% Description: Creates a parameter matrix that is used by RunHardRod

function [paramMat, numRuns] = MakeParamMat( systemObj, particleObj, runObj, rhoInit, flags )
%Find number parameters
numbc = length(systemObj.bc);
numvD = length(particleObj.vD);
numIC = length(rhoInit.IntCond);
numSM = length(flags.StepMeth);
numTr = runObj.num_trial;

% Handls all grid points equal and square box different
if flags.AllNsSame == 1
  Nvec = unique( [systemObj.Nx systemObj.Ny systemObj.Nm] );
  numN = length(Nvec);
  numNx = numN;
  numNy = 1;
  numNm = 1;
else
  numNx = length(systemObj.Nx);
  numNy = length(systemObj.Ny);
  numNm = length(systemObj.Nm);
end

if flags.SquareBox == 1
  Lvec = unique( [ParamMaster.Lx ParamMaster.Ly] );
  numL = length(Lvec);
  numLx = numL;
  numLy = 1;
else
  numLx = length(systemObj.Lx);
  numLy = length(systemObj.Ly);
end

% number of parameters
numNxNy = numNx*numNy;
numNxNyNm = numNxNy * numNm;
numNxNyNmLx = numNxNyNm * numLx;
numNxNyNmLxLy = numNxNyNmLx * numLy;
numNxNyNmLxLyvD = numNxNyNmLxLy * numvD;
numParams = numNxNyNmLxLyvD * numbc;
numParamsIC = numParams * numIC;
numParamsICSM = numParamsIC * numSM;
numRuns = numParamsICSM * numTr;

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
                  for r = 1:numTr
                    rowInd =  1 + (i-1) + (j-1) * numNx + ( k-1) * numNxNy + ...
                      (l-1) * numNxNyNm + (m-1) * numNxNyNmLx + ...
                      (n-1) * numNxNyNmLxLy + (o-1) * numNxNyNmLxLyvD + ...
                      (p-1) * numParams + (q-1) * numParamsIC +...
                      (r-1) * numParamsICSM;
                    
                    % Special cases
                    if flags.AllNsSame
                      paramMat(rowInd,1:3) = Nvec(i);
                    else
                      paramMat(rowInd,1) = systemObj.Nx(i);
                      paramMat(rowInd,2) = systemObj.Ny(j);
                      paramMat(rowInd,3) = systemObj.Nm(k);
                    end
                    
                    if flags.SquareBox
                      paramMat(rowInd,4:5) = Lvec(l);
                    else
                      paramMat(rowInd,4) = systemObj.Lx(l);
                      paramMat(rowInd,5) = systemObj.Ly(m);
                    end
                    % Everything else
                    paramMat(rowInd,6) = particleObj.vD(n);
                    paramMat(rowInd,7) = systemObj.bc(o);
                    paramMat(rowInd,8) = rhoInit.IntCond(p);
                    paramMat(rowInd,9) = flags.StepMeth(q);
                    paramMat(rowInd,10) = runObj.runID + (r - 1);
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


