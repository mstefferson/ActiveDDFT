
% Date
DateTimeStart =  datestr(now);
fprintf('Starting RunHardRod: %s\n', DateTimeStart);

% Print what you are doing
if FlagObj.AnisoDiff  == 1;
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

NumNx = length(ParamObj.Nx);
NumNy = length(ParamObj.Ny);
NumNm = length(ParamObj.Nm);
NumLx = length(ParamObj.Lx);
NumLy = length(ParamObj.Ly);
Numbc = length(ParamObj.bc);
NumvD = length(ParamObj.vD);
NumParams = NumNx * NumNy * NumNm * NumLx * NumLy * Numbc * NumvD;


if  strcmp( IntConcStr,'PlaneWaveEq' ) || strcmp( IntConcStr,'SepPWeq' )
    if 1.499 < bc && bc < 1.501
        bc = 1.502;
    end
end

if strcmp( IntConcStr,'Loaded' )
    IntDenName = sprintf('DenBlow2');
else
    IntDenName = sprintf('Irrelevant');
end

% Make the output directory string and input file
 if SaveMe
     if ~exist('Outputs', 'dir'); mkdir('Outputs'); end;
     DiaryStr = sprintf('DiarySingRunt%d.log',trial);
     diary(DiaryStr);
     runfile  = 'AnisRunLog.log';
     rlId = fopen(runfile,'a+'); 
%      fprintf(rlId, 'Anisotropic Hard Rod Hard Rod \n');
    
 end

FileDir = ...
    sprintf('HrAnD%i%i%i_bc%.2f_Int%i_v%.1f_IC%dsm%dt%d',...
    Nx,Ny,Nm,bc,Interactions,vD,IntCond,StepMeth,trial);
FileInpt = ...
    sprintf('InptAnD_N%i%i%i_bc%.2f_Int%i_v%.1f_IC%dsm%dt%d.txt', ...
    Nx,Ny,Nm,bc,Interactions,vD,IntCond,StepMeth,trial);

Where2SavePath    = sprintf('%s/%s/%s',pwd,'Outputs',FileDir);

% Figure out Total Runs
    tic
    %     keyboard
    [DenFinal, DenFTFinal, GridObj, ParamMaster,TimeObj,...
        DidIBreak,SteadyState,MaxReldRho] = ...
        HR2DrotMain(FileInpt);
    toc
    disp('Params');disp(ParamMaster);disp('Time');disp(TimeObj);      
    fprintf('Break = %d Steady = %d Max dRho/Rho = %.2e\n',...
        DidIBreak,SteadyState,MaxReldRho)
    DateTimeEnd =  datestr(now);
    fprintf('End: %s\n\n', DateTimeEnd);
    
    if SaveMe
        diary off; 
        if ~exist(Where2SavePath, 'dir'); mkdir(Where2SavePath); end
        movefile('*.mat', Where2SavePath)
        movefile('*.txt', Where2SavePath)
        movefile('Diary*', Where2SavePath)
        if MakeMovies
            movefile('*.avi', Where2SavePath)
        end
        if MakeOP
            movefile('*.fig', Where2SavePath)
            movefile('*.jpg', Where2SavePath)
        end       
        
        rlId = fopen(runfile,'a+');
        fprintf(rlId,'%s\n', FileInpt(1:end-4));
        fprintf(rlId,'( dt t_rec t_tot ss_eps ): ( %.2e %.2e %.2e %.2e )\n',...
            Timetmp);
        fprintf(rlId,'( Nx, Ny, Nm ): ( %d, %d, %d )\n', Nx, Ny, Nm);
        fprintf(rlId,'( Lx, Ly, Lrod ) = ( %.1f, %.1f, %.1f )\n',...
         Lx, Ly, L_rod);
        fprintf(rlId,...
            '( bc, vD, Mob_par, Mob_perp, Mob_rot ) = ( %.2f, %.2f, %.1f, %.1f, %.1f )\n',...
         bc, vD, Mob_par, Mob_perp, Mob_rot);
        fprintf(rlId,'( MdX, MdY, MdM, rand) = ( %d, %d, %d, %d )\n',...
         NumModesX, NumModesY, NumModesM, RandomAmp);
        fprintf(rlId,'Break = %d Steady = %d Max dRho/Rho = %.2e\n',...
        DidIBreak,SteadyState,MaxReldRho);
        fprintf(rlId, 'Start: %s\n', DateTimeStart);
        fprintf(rlId,'End: %s\n\n', DateTimeEnd);
        fclose('all');
    else
        remove('*.txt')
    end

end

