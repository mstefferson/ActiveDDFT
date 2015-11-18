% HR2DrotMainDrIDCube
%
% Program is the main() for running the diffusion of 2D hard rods with
% orientation. Program handles interactions using DDFT
%
% Angle-indepent diffusion matrix. Approximate interactions.

function [DenFinal, DenFTFinal, GridObj, ParamObj,TimeObj,...
    DidIBreak,SteadyState,MaxReldRho] = ...
    HR2DrotMainDrIDCube(InputFile)
% Add paths (this should already be added, but just to be careful)
% Save error messages in file
try
    EvolvedDen = 0;DenFinal = 0;DenFTFinal = 0;GridObj = 0;ParamObj = 0;
    TimeObj = 0;DidIBreak = 0;SteadyState = 0;MaxReldRho = 0;
    %         keyboard
    tMainID  = tic;
    
    %Grab the parameters
    % keyboard
    DataTemp    = importdata(InputFile);
    ParamVec    = DataTemp.data(1,:);
    TimeVec     = DataTemp.data(2,~isnan(DataTemp.data(2,:)));    %Pull out the time stuff
    FileNameMat = DataTemp.textdata(1);
    Path2Save   = DataTemp.textdata(2);
    IntDenType  = DataTemp.textdata(3);
    %     keyboard
    % Make some  objects
    ParamNmVec = {'trial' 'Interactions' 'Drive' 'MakeOP' 'MakeMovies' 'SaveMe'...
        'Nx' 'Ny' 'Nm' 'Lx' 'Ly' 'L_rod' 'Diam' 'Eta_visc'...
        'kB' 'Tmp' 'Norm' 'WeightPos' 'WeightAng' 'NumModesX' 'NumModesY' 'NumModesM' 'bc'...
        'Mob_pos','Mob_rot'  };
    TimeNmVec  = {'delta_t' 't_record' 't_tot' 'ss_epsilon'};
    
    ParamObj   = struct('NmVec',{ParamNmVec},'ValVec',...
        ParamVec,'trial',ParamVec(1),...
        'Interactions',ParamVec(2), 'Drive',ParamVec(3),...
        'MakeOP',ParamVec(4),'MakeMovies',ParamVec(5),'SaveMe',ParamVec(6),...
        'Nx', ParamVec(7),'Ny', ParamVec(8),'Nm', ParamVec(9),...
        'Lx', ParamVec(10),'Ly', ParamVec(11),'L_rod', ParamVec(12), ...
        'Tmp',ParamVec(13), 'Norm',ParamVec(14), 'WeightPos',ParamVec(15), ...
        'WeightAng',ParamVec(16), 'Random',ParamVec(17),...
        'NumModesX',ParamVec(18), 'NumModesY',ParamVec(19), ...
        'NumModesM', ParamVec(20),'bc',ParamVec(21), ...
        'c', ParamVec(22), ...
        'Mob_pos',ParamVec(23),'Mob_rot',ParamVec(24), 'v0', ParamVec(25) );
    %     keyboard
    % Create a file that holds warning print statements
    WarningStmtString = sprintf('WarningStmts_%i.txt',ParamObj.trial);
    wfid              = fopen(WarningStmtString,'a+');    % a+ allows to append data
    
    LocString = sprintf('Location_%i.txt',ParamObj.trial);
    lfid      = fopen(LocString,'a+');    % a+ allows to append data
    fprintf(lfid,'Starting main, current code\n');
    
    %ime Recording
    N_time   = ceil(TimeVec(3)/TimeVec(1)); %number of time steps
    N_record = ceil(TimeVec(3)/TimeVec(2)); %number of time points to record. Does not include initial density
    N_count  = ceil(TimeVec(2)/TimeVec(1)); %spacing between times to record
    
    TimeObj = struct('TimeNmVec',{TimeNmVec},'TimeVecOrg',{TimeVec},...
        'delta_t',TimeVec(1), 't_record',TimeVec(2), 't_tot',TimeVec(3), 'ss_epsilon',TimeVec(4),...
        'N_time', N_time, 'N_record',N_record,'N_count',N_count);
   
    % Fix the time
   [TimeObj.t_tot,TimeObj.N_time,TimeObj.t_rec,TimeObj.N_rec,TimeObj.N_count]= ...
      TimeStepRecMaker(TimeObj.delta_t,TimeObj.t_tot,TimeObj.t_record);

    %%%Make all the grid stuff%%%%%%%%%%%%%%
    [GridObj] = GridMakerPBCxk(...
        ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);
    fprintf(lfid,'Made grid\n');
    
    %Make diffusion coeff (send smallest dx dy for stability
    [DiffMobObj] = DiffMobCoupCoeffCalcIsoDiff(...
        ParamObj.Tmp,ParamObj.Mob_pos,ParamObj.Mob_rot);
    fprintf(lfid,'Made diffusion object\n');
    
    %Initialze density
    [rho] = MakeConcFromInd(GridObj,ParamObj,IntDenType);
    Nc    = 20;
 % Equilib distribution
    [Coeff_best,~] = CoeffCalcExpCos2D(Nc,GridObj.phi,ParamObj.bc); % Calculate coeff
    feq = DistBuilderExpCos2Dsing(Nc,GridObj.phi,Coeff_best);        % Build equil distribution

%     keyboard
    % Run the main code
    tBodyID      = tic;
    
    [DenRecObj]  = HR2DrotDenEvolverFTBodyIDCube(...
        wfid,lfid,rho,ParamObj, TimeObj,GridObj,DiffMobObj,feq);
    EvolvedDen = 1;
    BodyRunTime  = toc(tBodyID);
    fprintf(lfid,'Made density object\n');
    fprintf(lfid,'Body Run Time = %f\n\n', BodyRunTime);
    
    % Store final density and transform
    DenFinal   = DenRecObj.Density_rec(:,:,:,end);
    DenFTFinal = DenRecObj.DensityFT_rec(:,:,:,end);
    DidIBreak  = DenRecObj.DidIBreak;
    SteadyState = DenRecObj.SteadyState;
    MaxReldRho  = DenRecObj.MaxReldRho;
    
%         keyboard

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
        OpRunTime       = toc(tOpID);
        fprintf(lfid,'Made interaction order paramater object\n');
        
        if ParamObj.MakeMovies == 1
            % Build OP records
            
            % Make matlab movies
            tMovID       = tic;
%             keyboard
        OPMatMovieMakerTgthr(GridObj,ParamObj,OrderParamObj,...
            DenRecObj.Density_rec,feq);
            %                 keyboard
            MovRunTime   = toc(tMovID);
            %         keyboard
            fprintf(lfid,'Made movies\n');
            % Record how long it took
            fprintf(lfid,'OrderParam Run time = %f\n', OpRunTime);
            fprintf(lfid,'Make Mov Run Time = %f\n',  MovRunTime);
        end % End if movies        
    end % if OP
    
    if ParamObj.SaveMe
        MemObj = 0;
        % Save all parameters
        
        % Save everything. Save seperately for big files
        DenStr = sprintf('DenRec_%i',ParamObj.trial);
        TimeStr = sprintf('TimeObj_%i',ParamObj.trial);
        ParamStr = sprintf('ParamObj_%i',ParamObj.trial);
        GridStr = sprintf('GridObj_%i',ParamObj.trial);
        
        save(DenStr,'DenRecObj','-v7.3')
        save(TimeStr,'GridObj','-v7.3')
        save(ParamStr,'ParamObj','-v7.3')
        save(GridStr,'GridObj','-v7.3')
       
        if ParamObj.MakeOP
            OpStr = sprintf('OP_%i',ParamObj.trial);
          save(OpStr,'OrderParamObj','-v7.3')
        end
    end
    % Save how long everything took
    fprintf(lfid,'Everything saved\n');
    TotRunTime = toc(tMainID);
    fprintf(lfid,'Total Run time = %f\n', TotRunTime);
    
    fclose('all');
    % keyboard
    %     if ParamObj.SaveMe
    % Move everything
    %         MoveStrTxt = sprintf('*%i.txt', ParamObj.trial);
    %         MoveStrMat = sprintf('*%i.mat', ParamObj.trial);
    %         movefile(MoveStrTxt,Path2Save{1});
    %         movefile(MoveStrMat,Path2Save{1});
    %         cd /home/mws/Documents/Research/BG/DDFT/Outputs
    %     end
    
catch err %Catch errors 
            keyboard
    
    ErrFileNmStr = sprintf('errFile%i.txt',ParamObj.trial);
    efid         = fopen(ErrFileNmStr,'a+');
    % write the error to file and to screen
    % first line: message
%     fprintf(efid,'%s', err.getReport('extended', 'hyperlinks','off')) ;
    fprintf('%s', err.getReport('extended')) ;
    disp(err.message);
    fclose(efid);
    fclose('all');
    
%    keyboard
    if ParamObj.SaveMe
      
        TimeStr = sprintf('TimeObj_%i',ParamObj.trial);
        ParamStr = sprintf('ParamObj_%i',ParamObj.trial);
        GridStr = sprintf('GridObj_%i',ParamObj.trial);
       
        save(TimeStr,'GridObj','-v7.3')
        save(ParamStr,'ParamObj','-v7.3')
        save(GridStr,'GridObj','-v7.3')
        if EvolvedDen
           DenStr = sprintf('DenRec_%i',ParamObj.trial);     
           save(DenStr,'DenRecObj','-v7.3');
        end
    end
    
end %End try and catch

% clc
close all
end % End HR2DrotVgrExeMain.m
