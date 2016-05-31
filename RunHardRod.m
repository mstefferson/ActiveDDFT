
% Date
DateTimeStart =  datestr(now);
fprintf('Starting RunHardRod: %s\n', DateTimeStart);

% Print what you are doing
if Flags.AnisoDiff  == 1;
  fprintf('Anisotropic Hard Rod \n')
else
  fprintf('Isotropic Hard Rod \n')
end

% Grab initial parameters
if exist('Params.mat','file') == 0;
  if exist('InitParams.m','file') == 0;
    cpmatparams
  end;
  InitParams
end
load Params.mat;
% Copy the master parameter list to ParamObj
ParamObj = ParamMaster;
RhoInit  = RhoInitMaster;
Flags    = FlagMaster;


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

% For some reason, param_mat gets "sliced". Create vectors to get arround
paramNx  = paramMat(:,1);
paramNy  = paramMat(:,2);
paramNm  = paramMat(:,3);
paramLx  = paramMat(:,4);
paramLy  = paramMat(:,5);
paramvD  = paramMat(:,6);
parambc  = paramMat(:,7);
paramIC  = paramMat(:,8);
paramSM  = paramMat(:,9);
paramrun = paramMat(:,10);

% Loops over all run

for ii = 1:numRuns
  % Assign parameters
  ParamObj.Nx = paramNx(ii);
  ParamObj.Ny = paramNy(ii);
  ParamObj.Nm = paramNm(ii);
  ParamObj.Lx = paramLx(ii);
  ParamObj.Ly = paramLy(ii);
  ParamObj.vD = paramvD(ii);
  ParamObj.bc = parambc(ii);
  RhoInit.IntCond = paramIC(ii);
  Flags.StepMeth = paramSM(ii);
  ParamObj.runID = paramrun(ii);
  ParamObj.Norm  = ParamObj.bc / ( ParamObj.L_rod ^ 2 / pi ) *... % number of particles
    ParamObj.Lx * ParamObj.Ly;
  if Flags.SquareBox
    ParamObj.L_box = ParamObj.Lx;
  else
    ParamObj.L_box = [ParamObj.Lx ParamObj.Ly];
  end
  
  % Name the file
  filename = [ 'Hr_Ani' num2str( Flags.AnisoDiff ) ...
    '_N' num2str( ParamObj.Nx ) num2str( ParamObj.Ny ) num2str( ParamObj.Nm )  ...
    '_Lx' num2str( ParamObj.Lx ) 'Ly' num2str( ParamObj.Ly )...
    '_vD' num2str( ParamObj.vD ) '_bc' num2str( ParamObj.bc ) ...
    '_IC' num2str( RhoInit.IntCond ) '_SM' num2str( Flags.StepMeth ) ...
    '_t' num2str( ParamObj.trial ) '.' num2str( ParamObj.runID ) '.mat' ];
  
  disp(filename);
  disp(ParamObj);
  
end
%%
%if  strcmp( IntConcStr,'PlaneWaveEq' ) || strcmp( IntConcStr,'SepPWeq' )
%if 1.499 < bc && bc < 1.501
%bc = 1.502;
%end
%end


%% Make the output directory string and input file
%if SaveMe
%if ~exist('Outputs', 'dir'); mkdir('Outputs'); end;
%DiaryStr = sprintf('DiarySingRunt%d.log',trial);
%diary(DiaryStr);
%runfile  = 'AnisRunLog.log';
%rlId = fopen(runfile,'a+');
%%      fprintf(rlId, 'Anisotropic Hard Rod Hard Rod \n');
%end

%FileDir = ...
%sprintf('HrAnD%i%i%i_bc%.2f_Int%i_v%.1f_IC%dsm%dt%d',...
%Nx,Ny,Nm,bc,Interactions,vD,IntCond,StepMeth,trial);
%FileInpt = ...
%sprintf('InptAnD_N%i%i%i_bc%.2f_Int%i_v%.1f_IC%dsm%dt%d.txt', ...
%Nx,Ny,Nm,bc,Interactions,vD,IntCond,StepMeth,trial);

%Where2SavePath    = sprintf('%s/%s/%s',pwd,'Outputs',FileDir);

%% Figure out Total Runs
%tic
%%     keyboard
%[DenFinal, DenFTFinal, GridObj, ParamMaster,TimeObj,...
%DidIBreak,SteadyState,MaxReldRho] = ...
%HR2DrotMain(FileInpt);
%toc
%disp('Params');disp(ParamMaster);disp('Time');disp(TimeObj);
%fprintf('Break = %d Steady = %d Max dRho/Rho = %.2e\n',...
%DidIBreak,SteadyState,MaxReldRho)
%DateTimeEnd =  datestr(now);
%fprintf('End: %s\n\n', DateTimeEnd);

%if SaveMe
%diary off;
%if ~exist(Where2SavePath, 'dir'); mkdir(Where2SavePath); end
%movefile('*.mat', Where2SavePath)
%movefile('*.txt', Where2SavePath)
%movefile('Diary*', Where2SavePath)
%if MakeMovies
%movefile('*.avi', Where2SavePath)
%end
%if MakeOP
%movefile('*.fig', Where2SavePath)
%movefile('*.jpg', Where2SavePath)
%end

%rlId = fopen(runfile,'a+');
%fprintf(rlId,'%s\n', FileInpt(1:end-4));
%fprintf(rlId,'( dt t_rec t_tot ss_eps ): ( %.2e %.2e %.2e %.2e )\n',...
%Timetmp);
%fprintf(rlId,'( Nx, Ny, Nm ): ( %d, %d, %d )\n', Nx, Ny, Nm);
%fprintf(rlId,'( Lx, Ly, Lrod ) = ( %.1f, %.1f, %.1f )\n',...
%Lx, Ly, L_rod);
%fprintf(rlId,...
%'( bc, vD, Mob_par, Mob_perp, Mob_rot ) = ( %.2f, %.2f, %.1f, %.1f, %.1f )\n',...
%bc, vD, Mob_par, Mob_perp, Mob_rot);
%fprintf(rlId,'( MdX, MdY, MdM, rand) = ( %d, %d, %d, %d )\n',...
%NumModesX, NumModesY, NumModesM, RandomAmp);
%fprintf(rlId,'Break = %d Steady = %d Max dRho/Rho = %.2e\n',...
%DidIBreak,SteadyState,MaxReldRho);
%fprintf(rlId, 'Start: %s\n', DateTimeStart);
%fprintf(rlId,'End: %s\n\n', DateTimeEnd);
%fclose('all');
%else
%remove('*.txt')
%end

%end

