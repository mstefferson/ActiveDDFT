% HR2DrotMain.m
%
% Program is the main() for running the diffusion of 2D hard rods with
% orientation. Program handles interactions using DDFT

function  [ DidIBreak,SteadyState,MaxReldRho] = ...
    HR2DrotMain( filename, ParamObj, TimeObj, RhoInit, Flags ) 
% Add paths (this should already be added, but just to be careful)
% Save error messages in file
try
    EvolvedDen = 0;DenFinal = 0;DenFTFinal = 0;GridObj = 0;
    DidIBreak = 0;SteadyState = 0;MaxReldRho = 0;
    
    %     keyboard
    tMainID  = tic;

    % Create a file that holds warning print statements
    WarningStmtString = sprintf('WarningStmts_%i.txt',ParamObj.trial);
    wfid  = fopen(WarningStmtString,'a+');    % a+ allows to append data
    
    LocString = sprintf('Location_%i.txt',ParamObj.trial);
    lfid      = fopen(LocString,'a+');    % a+ allows to append data
    fprintf(lfid,'Starting main, current code\n');
    
    % Make remaining objects
    % Make all the grid stuff % 
    tGridID = tic;
    [GridObj] = GridMakerPBCxk(...
        ParamObj.Nx,ParamObj.Ny,ParamObj.Nm,ParamObj.Lx,ParamObj.Ly);
    fprintf(lfid,'Made Grid\n');
    GridRunTime = toc(tGridID);
    fprintf('Made grid: %.3g \n', GridRunTime);
  
    %Make diffusion coeff (send smallest dx dy for stability
    tDiffID = tic;
    [DiffMobObj] =  DiffMobCoupCoeffCalc( wfid,ParamObj.Tmp,...
        ParamObj.Mob_par,ParamObj.Mob_perp,ParamObj.Mob_rot,...
        TimeObj.delta_t, min(GridObj.dx,GridObj.dy),...
        GridObj.dphi,GridObj.kx2D, GridObj.ky2D,ParamObj.vD);
    
    fprintf(lfid,'Made diffusion object\n');
    DiffRunTime = toc(tDiffID);
    fprintf('Made diffusion object: %.3g\n', DiffRunTime);
 
    %Initialze density
    tIntDenID = tic;
    [rho] = MakeConc(GridObj,ParamObj,RhoInit);
    Nc    = 20;
    % Equilib distribution. Don't let bc = 1.5
    if 1.499 < ParamObj.bc && ParamObj.bc < 1.501
      RhoInit.bc = 1.502;
    else
      RhoInit.bc = ParamObj.bc;
    end
    [Coeff_best,~] = CoeffCalcExpCos2D(Nc,GridObj.phi,RhoInit.bc); % Calculate coeff
    RhoInit.feq = DistBuilderExpCos2Dsing(Nc,GridObj.phi,Coeff_best);        % Build equil distribution
    fprintf(lfid,'Made initial density\n');
    IntDenRunTime = toc(tIntDenID);
    fprintf('Made initial density: %.3g \n', IntDenRunTime);

    % Run the main code
    tBodyID      = tic;
    
    [DenRecObj]  = HR2DrotDenEvolverFTBody(...
    wfid, lfid, rho, ParamObj, TimeObj, GridObj,DiffMobObj, Flags, RhoInit.feq);
    %     keyboard
    EvolvedDen = 1;
    fprintf(lfid,'Ran Main Body\n');
    BodyRunTime  = toc(tBodyID);
    fprintf('Ran Main Body: %.3g \n', BodyRunTime);
    
    fprintf(lfid,'Body Run Time = %f\n\n', BodyRunTime);
    
    % Store final density and transform
    DenFinal   = DenRecObj.Density_rec(:,:,:,end);
    DenFTFinal = DenRecObj.DensityFT_rec(:,:,:,end);
    DidIBreak  = DenRecObj.DidIBreak;
    SteadyState = DenRecObj.SteadyState;
    MaxReldRho  = DenRecObj.MaxReldRho;
    
    % Run movies if you want
    if Flags.MakeOP  == 1
        tOpID           = tic ;
        %                 keyboard
        if  DenRecObj.DidIBreak == 0
            [OrderParamObj] = CPNrecMaker(...
                ParamObj.Nx,ParamObj.Ny,DenRecObj.TimeRecVec,...
                GridObj,DenRecObj.Density_rec,RhoInit.feq);
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
        
        if Flags.MakeMovies == 1
            % Build OP records
            
            % Make matlab movies
            tMovID       = tic;
            %         keyboard
            HoldX = ParamObj.Nx /2 + 1;
            HoldY = ParamObj.Ny /2 + 1;
            
            DistRec =  reshape( DenRecObj.Density_rec(HoldY, HoldY, : , :),...
                    [ParamObj.Nm length(DenRecObj.TimeRecVec)] );
            
            if DenRecObj.DidIBreak == 0
                
                OPMovieMakerTgtherDirAvi(ParamObj.trial,...
                    GridObj.x,GridObj.y,GridObj.phi,OrderParamObj,...
                    DistRec,OrderParamObj.TimeRec);
            
            else
                
                OPMovieMakerTgtherDirAvi(ParamObj.trial,...
                    GridObj.x,GridObj.y,GridObj.phi,OrderParamObj,...
                    DistRec,OrderParamObj.TimeRec(1:end-1) );
            end

            fprintf(lfid,'Made movies\n');
            MovRunTime   = toc(tMovID);
            fprintf('Made movies: %.3g \n', MovRunTime);
            %disp(MovRunTime);
            % Record how long it took
            fprintf(lfid,'Make Mov Run Time = %f\n',  MovRunTime);      

        end % End if movies
        
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
            ParamObj.Nm, DenRecObj.bc,ParamObj.vD,  Flags.SaveMe, ParamObj.trial)
   
    end % if OP
    
    if Flags.SaveMe
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
        
        if Flags.MakeOP
            OpStr = sprintf('OP_%i',ParamObj.trial);
            save(OpStr,'OrderParamObj','-v7.3')
        end
    end
    % Save how long everything took
    % Save how long everything took
    fprintf(lfid,'Everything saved. Run finished\n');
    TotRunTime = toc(tMainID);
    fprintf('Everything saved. Run Finished: %.3g \n', TotRunTime);
    %disp(TotRunTime);
    fprintf(lfid,'Total Run time = %f\n', TotRunTime);
    
    fclose('all');
    
catch err %Catch errors
    
    
    ErrFileNmStr = sprintf('errFile%i.txt',ParamObj.trial);
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
    if Flags.SaveMe
        
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
%close all
end % End HR2DrotVgrExeMain.m
