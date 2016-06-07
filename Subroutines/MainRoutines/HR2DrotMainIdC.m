% HR2DrotMainIdC
%
% Program is the main() for running the diffusion of 2D hard rods with
% orientation. Program handles interactions using DDFT. 
% 
% Isotropic diffusion cube form.
%
% Angle-indepent diffusion matrix. Approximate interactions.

function [DenFinal, DenFTFinal, GridObj, ParamObj,TimeObj,...
    DidIBreak,SteadyState,MaxReldRho] = ...
    HR2DrotMainIdC(InputFile)
% Add paths (this should already be added, but just to be careful)
% Save error messages in file
try
    EvolvedDen = 0;DenFinal = 0;DenFTFinal = 0;GridObj = 0;ParamObj = 0;
    TimeObj = 0;DidIBreak = 0;SteadyState = 0;MaxReldRho = 0;
    %         keyboard
    tMainID  = tic;
    
    %Grab the parameters
    tParamID = tic;
    % keyboard
    DataTemp    = importdata(InputFile);
    ParamVec    = DataTemp.data(1,:);
    TimeVec     = DataTemp.data(2,~isnan(DataTemp.data(2,:)));    %Pull out the time stuff
    FileNameMat = DataTemp.textdata(1);
    Path2Save   = DataTemp.textdata(2);
    IntDenType  = DataTemp.textdata(3);
    LoadName    = DataTemp.textdata(4);

    % Make some  objects
    ParamNmVec = {'trial' 'Interactions' 'Drive'  'StepMeth' 'IntCond' ...
        'MakeOP' 'MakeMovies' 'SaveMe'...
        'Nx' 'Ny' 'Nm' 'Lx' 'Ly' 'L_rod' 'Diam' 'Eta_visc'...
        'kB' 'Tmp' 'Norm' 'WeightPos' 'WeightAng' 'NumModesX' 'NumModesY' 'NumModesM' 'bc'...
        'Mob_pos','Mob_rot'  };
    TimeNmVec  = {'dt' 't_record' 't_tot' 'ss_epsilon'};
    
    ParamObj   = struct('NmVec',{ParamNmVec},'ValVec',...
        ParamVec,'trial',ParamVec(1),...
        'Interactions',ParamVec(2), 'Drive',ParamVec(3),...
        'StepMeth', ParamVec(4), 'IntCond', ParamVec(5), ...
        'MakeOP',ParamVec(6),'MakeMovies',ParamVec(7),'SaveMe',ParamVec(8),...
        'Nx', ParamVec(9),'Ny', ParamVec(10),'Nm', ParamVec(11),...
        'Lx', ParamVec(12),'Ly', ParamVec(13),'L_rod', ParamVec(14), ...
        'Tmp',ParamVec(15), 'Norm',ParamVec(16), 'WeightPos',ParamVec(17), ...
        'WeightAng',ParamVec(18), 'Random',ParamVec(19),...
        'NumModesX',ParamVec(20), 'NumModesY',ParamVec(21), ...
        'NumModesM', ParamVec(22),'bc',ParamVec(23), ...
        'c', ParamVec(24), ...
        'Mob_pos',ParamVec(25),'Mob_rot',ParamVec(26), 'vD', ParamVec(27) );
   
    ParamRunTime = toc(tParamID);
    fprintf('Read input and made ParamObj: %.3g \n',ParamRunTime)
    %disp(ParamRunTime);

    % Create a file that holds warning print statements
    WarningStmtString = sprintf('WarningStmts_%i.txt',ParamObj.trialID);
    wfid              = fopen(WarningStmtString,'a+');    % a+ allows to append data
    
    LocString = sprintf('Location_%i.txt',ParamObj.trialID);
    lfid      = fopen(LocString,'a+');    % a+ allows to append data
    fprintf(lfid,'Starting main, current code\n');
    
    % Time Recording
    tTimeID = tic;
    N_time   = ceil(TimeVec(3)/TimeVec(1)); %number of time steps
    N_record = ceil(TimeVec(3)/TimeVec(2)); %number of time points to record. Does not include initial density
    N_count  = ceil(TimeVec(2)/TimeVec(1)); %spacing between times to record
    
    TimeObj = struct('TimeNmVec',{TimeNmVec},'TimeVecOrg',{TimeVec},...
        'dt',TimeVec(1), 't_record',TimeVec(2), 't_tot',TimeVec(3), 'ss_epsilon',TimeVec(4),...
        'N_time', N_time, 'N_record',N_record,'N_count',N_count);
    
    % Fix the time
    [TimeObj.t_tot,TimeObj.N_time,TimeObj.t_rec,TimeObj.N_rec,TimeObj.N_count]= ...
        TimeStepRecMaker(TimeObj.dt,TimeObj.t_tot,TimeObj.t_record);
    fprintf(lfid,'Made Time Obj\n');
    TimeRunTime = toc(tTimeID);
    fprintf('Made Time Obj: %.3g \n', TimeRunTime)
    %disp(TimeRunTime);

    %%%Make all the grid stuff%%%%%%%%%%%%%%
    tGridID = tic;
    [GridObj] = GridMakerPBCxk(...
        ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);
    fprintf(lfid,'Made grid\n');
    GridRunTime = toc(tGridID);
    fprintf('Made grid: %.3g \n', GridRunTime);
    %disp(GridRunTime);
    
    %Make diffusion coeff (send smallest dx dy for stability
    tDiffID = tic;
    [DiffMobObj] = DiffMobCoupCoeffCalcIsoDiff(...
        ParamObj.Tmp,ParamObj.Mob_pos,ParamObj.Mob_rot);
    fprintf(lfid,'Made diffusion object\n');
    DiffRunTime = toc(tDiffID);
    fprintf('Made diffusion object: %.3g\n', DiffRunTime);
    %disp(DiffRunTime);
    
    %Initialze density
    tIntDenID = tic;
    [rho] = MakeConc(GridObj,ParamObj);
    Nc    = 20;
    % Equilib distribution
    [Coeff_best,~] = CoeffCalcExpCos2D(Nc,GridObj.phi,ParamObj.bc); % Calculate coeff
    feq = DistBuilderExpCos2Dsing(Nc,GridObj.phi,Coeff_best);        % Build equil distribution
    fprintf(lfid,'Made initial density\n');
    IntDenRunTime = toc(tIntDenID);
    fprintf('Made initial density: %.3g \n', IntDenRunTime);
   %disp(IntDenRunTime);

    % Run the main code
    tBodyID      = tic;
    [DenRecObj]  = HR2DrotDenEvolverFTBodyIdC(...
        wfid,lfid,rho,ParamObj, TimeObj,GridObj,DiffMobObj,feq);
    EvolvedDen = 1;
    fprintf(lfid,'Ran Main Body\n');
    BodyRunTime  = toc(tBodyID);
    fprintf('Ran Main Body: %.3g \n', BodyRunTime);
    %disp(BodyRunTime);
    fprintf(lfid,'Body Run Time = %f\n\n', BodyRunTime);
    
    % Store final density and transform
    DenFinal   = DenRecObj.Density_rec(:,:,:,end);
    DenFTFinal = DenRecObj.DensityFT_rec(:,:,:,end);
    DidIBreak  = DenRecObj.DidIBreak;
    SteadyState = DenRecObj.SteadyState;
    MaxReldRho  = DenRecObj.MaxReldRho;
    
    % Run movies if you want
    if ParamObj.MakeOP  == 1
        tOpID           = tic ;
        %                 keyboard
        if  DenRecObj.DidIBreak == 0
            [OrderParamObj] = CPNrecMaker(...
                ParamObj.Nx,ParamObj.Ny,DenRecObj.TimeRecVec,...
                GridObj,DenRecObj.Density_rec,feq);
        else %Don't incldue the blowed up denesity for movies. They don't like it.
            TimeRecVecTemp = DenRecObj.TimeRecVec(1:end-1);
            [OrderParamObj] = CPNrecMaker(ParamObj.Nx,ParamObj.Ny,...
                TimeRecVecTemp,GridObj,...
                DenRecObj.Density_rec(:,:,:,1:length(TimeRecVecTemp)),...
                feq);
        end
        fprintf(lfid,'Made interaction order paramater object\n');
        OpRunTime = toc(tOpID);
        fprintf('Made interaction order paramater object: %.3g \n', OpRunTime);
        %disp(OpRunTime);
        fprintf(lfid,'OrderParam Run time = %f\n', OpRunTime);
        
        if ParamObj.MakeMovies == 1
            % Build OP records
            
            % Make matlab movies
            tMovID       = tic;
            %         keyboard
            HoldX = ParamObj.Nx /2 + 1;
            HoldY = ParamObj.Ny /2 + 1;
                    if DenRecObj.DidIBreak == 0
          
            OPMovieMakerTgtherAvi(ParamObj.trialID,GridObj.x,GridObj.y, GridObj.phi, ...
                OrderParamObj.C_rec, OrderParamObj.NOP_rec,OrderParamObj.POP_rec,...
                reshape( DenRecObj.Density_rec(HoldX, HoldY, : , :), ...
                [ParamObj.Nm length(DenRecObj.TimeRecVec)] ),...
                DenRecObj.TimeRecVec)
            
            else
                
                OPMovieMakerTgtherAvi(ParamObj.trialID,GridObj.x,GridObj.y, GridObj.phi, ...
                OrderParamObj.C_rec, OrderParamObj.NOP_rec,OrderParamObj.POP_rec,...
                reshape( DenRecObj.Density_rec(HoldX, HoldY, : ,1 :end - 1), ...
                [ParamObj.Nm length(DenRecObj.TimeRecVec) - 1] ),...
                DenRecObj.TimeRecVec(1:end - 1 ) )
                
            end
            
            fprintf(lfid,'Made movies\n');
            MovRunTime   = toc(tMovID);
            fprintf('Made movies: %.3g \n', MovRunTime);
            %disp(MovRunTime);
            % Record how long it took
            fprintf(lfid,'Make Mov Run Time = %f\n',  MovRunTime);      
        end % End if movies
      
             % Plot amps
        kx0 = ParamObj.Nx / 2 + 1;
        ky0 = ParamObj.Ny / 2 + 1;
        km0 = ParamObj.Nm / 2 + 1;
        Nrec = length( DenRecObj.TimeRecVec);

        FTind2plot = zeros( 8, 3 );
        FTmat2plot = zeros( 8, Nrec );
        
        FTind2plot(1,:) = [kx0     ky0     km0 + 1];
        FTind2plot(2,:) = [kx0 + 1 ky0     km0 + 1];
        FTind2plot(3,:) = [kx0     ky0 + 1 km0 + 1];
        FTind2plot(4,:) = [kx0 + 1 ky0 + 1 km0 + 1];
        FTind2plot(5,:) = [kx0     ky0     km0 + 2];
        FTind2plot(6,:) = [kx0 + 1 ky0     km0 + 2];
        FTind2plot(7,:) = [kx0     ky0 + 1 km0 + 2];
        FTind2plot(8,:) = [kx0 + 1 ky0 + 1 km0 + 2];
        
        for i = 1:8 
            FTmat2plot(i,:) =  reshape(... 
            DenRecObj.DensityFT_rec( FTind2plot(i,1), FTind2plot(i,2), FTind2plot(i,3),: ),...   
            [ 1, Nrec ]  );   
        end
        
%         keyboard
        ampPlotterFT(FTmat2plot, FTind2plot, DenRecObj.TimeRecVec, ParamObj.Nx, ParamObj.Ny,...
            ParamObj.Nm, DenRecObj.bc,ParamObj.vD,  ParamObj.SaveMe, ParamObj.trialID)
        
    end % if OP
    
    if ParamObj.SaveMe
        MemObj = 0;
        % Save all parameters
        
        % Save everything. Save seperately for big files
        DenStr = sprintf('DenRec_%i',ParamObj.trialID);
        TimeStr = sprintf('TimeObj_%i',ParamObj.trialID);
        ParamStr = sprintf('ParamObj_%i',ParamObj.trialID);
        GridStr = sprintf('GridObj_%i',ParamObj.trialID);
        
        save(DenStr,'DenRecObj','-v7.3')
        save(TimeStr,'GridObj','-v7.3')
        save(ParamStr,'ParamObj','-v7.3')
        save(GridStr,'GridObj','-v7.3')
        
        if ParamObj.MakeOP
            OpStr = sprintf('OP_%i',ParamObj.trialID);
            save(OpStr,'OrderParamObj','-v7.3')
        end
    end
    % Save how long everything took
    fprintf(lfid,'Everything saved. Run finished\n');
    TotRunTime = toc(tMainID);
    fprintf('Everything saved. Run Finished: %.3g \n', TotRunTime);
    %disp(TotRunTime);
    fprintf(lfid,'Total Run time = %f\n', TotRunTime);
    
    fclose('all');
    
catch err %Catch errors
    
    
    ErrFileNmStr = sprintf('errFile%i.txt',ParamObj.trialID);
    efid         = fopen(ErrFileNmStr,'a+');
    % write the error to file and to screen
    % first line: message
    %     fprintf(efid,'%s', err.getReport('extended', 'hyperlinks','off')) ;
    fprintf('%s', err.getReport('extended')) ;
    disp(err.message);
    fclose(efid);
    fclose('all');
    
    keyboard
    %    keyboard
    if ParamObj.SaveMe
        
        TimeStr = sprintf('TimeObj_%i',ParamObj.trialID);
        ParamStr = sprintf('ParamObj_%i',ParamObj.trialID);
        GridStr = sprintf('GridObj_%i',ParamObj.trialID);
        
        save(TimeStr,'GridObj','-v7.3')
        save(ParamStr,'ParamObj','-v7.3')
        save(GridStr,'GridObj','-v7.3')
        if EvolvedDen
            DenStr = sprintf('DenRec_%i',ParamObj.trialID);
            save(DenStr,'DenRecObj','-v7.3');
        end
    end
    
end %End try and catch

% clc
%close all
end % End HR2DrotVgrExeMain.m
