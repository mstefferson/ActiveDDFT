% PhaseTransFinder2D

try
    
    clear  
    
    % Add paths just in case
    CurrentDir = pwd;
    addpath( genpath( CurrentDir) );
    
    PlotMe = 0;
    SaveMe = 1;
    
    Nbc   = 5;   %Number of concentractions
    Ngr   = 2;   %Number of grid points to loop
    trial = 11;
    
    % Concentration vector % including bc  = 1.5 is ok because we perturb
    % about a nematic
    bc_start = 1.0;
    bc_end   = 1.6;
    
    % L_rod
    Ngr_end   = 48;
    N_approxstart = 32; % can change to make all Ns even and fit
    
    
    NxFix = 32;
    NyFix = 32;
    NmFix = 32;
    Nc    = 10; % Number of coefficients in fit
    % Driving
    v0 = 0;
    
    % L_rod sets length scale
    L_box = 1; % Should be a multiple of L_rod = 1;
    
    % Time
    t_tot = 1;
    
    % Build vecs
    bcVec = linspace(bc_start,bc_end,Nbc);
    bcVec = [1.4 1.4375 1.475 1.5125 1.5875]; 
   
    dN = 2 * floor( (Ngr_end-N_approxstart) / ( 2*(Ngr-1) ) );
    NgrVec = Ngr_end - dN*(Ngr-1) : dN : Ngr_end;
    
    disp('NgrVec');disp(NgrVec);
    disp('bcVec');disp(bcVec);
    % Make sure they are all even
    if sum(mod(NgrVec,2)) ~= 0
        error('Only even grid numbers')
    end
    
    % Just in case...
    if length(NgrVec) ~= Ngr
        Ngr = length(NgrVec);
    end
    if length(bcVec) ~= Nbc
        bc = length(bcVec);
    end
     
    % Initialize some things
    AveNemMatPde         = zeros( Ngr, Nbc );
    AveNemMatFit         = zeros( Ngr, Nbc );
    StdNemMatPde         = zeros( Ngr, Nbc );
    GotPhaseTrans = 0;
    
    tic
    for i = 1:Ngr
        Nvec = [NgrVec(i) NgrVec(i) NgrVec(i) ]; Nx = Nvec(1); Ny = Nvec(2);Nm = Nvec(3);
        %     Nvec = [NgrVec(i) NgrVec(i) NmFix ]; Nx = Nvec(1); Ny = Nvec(2);Nm = Nvec(3);
        % Nvec = [NxFix NyFix  NgrVec(i)]; Nx = Nvec(1); Ny = Nvec(2);Nm = Nvec(3);
        FinalDenRec = zeros(Nx,Ny,Nm,Nbc);
        %     keyboard
        ticTemp = tic;
        for j = 1:Nbc
            % Calculate the distribution, rho, and sigma at a specific concentration
            [FileInpt] = InptMkrHRICDfunc(bcVec(j),v0,L_box, Nvec,t_tot,trial);
            [DenFinal_new, DenFTFinal_new,GridObj, ParamObj, TimeObj, ...
                DidIBreak,SteadyState,MaxReldRho] = ....
                HR2DrotMainDrIDCube(FileInpt);
            fprintf('Ngr = %.2f bc = %.2f\n',NgrVec(i),bcVec(j))
            fprintf('DidIBreak = %d SteadyState = %d max dRho/Rho = %.2e \n',...
                DidIBreak,SteadyState,MaxReldRho);
            
            
            %Store
            FinalDenRec(:,:,:,j) = DenFinal_new;
            
            %Update
            DenFinal = DenFinal_new;
            
            %Build OP Matrix From PDE solution
            [~,~,~,~,NOP,~,~] = OpCPNCalc(Nx,Ny,DenFinal_new,GridObj.phi,GridObj.x,GridObj.y,GridObj.phi3D);
            AveNem = mean(mean(NOP));
            StdNem = std2(NOP);
            
            AveNemMatPde(i,j) = AveNem;
            StdNemMatPde(i,j) = StdNem;
            
            [Coeff_best, CoeffMat] = CoeffCalcExpCos2D(Nc,GridObj.phi,bcVec(j));
            f = DistBuilderExpCos2Dsing(Nc,GridObj.phi,Coeff_best);
            [~,~,~,~,NOP,~,~] = OpCPNCalc(1,1,f,GridObj.phi,0,0,GridObj.phi);
            AveNemMatFit (i,j) = NOP;
            
            if SaveMe
                if i == 1 && j == 1
                    mkdir Outputs
                    DiaryStr = sprintf('DiarySingRunt%d.txt',trial);
                    diary(DiaryStr);
                    disp('NgrVec');disp(NgrVec);disp('bcVec');disp(bcVec);
                    disp('Params');disp(ParamObj);disp('Time');disp(TimeObj);
                    diary off
                    SaveStr = sprintf('OPvsGridBCmatsBcNum%dNnum%dtr%d',...
                    Nbc,Ngr,trial);
                end
                
                save(SaveStr,'AveNemMatFit', 'AveNemMatPde','StdNemMatPde',...
                    'bcVec','NgrVec','L_box','v0', 'SaveMe','trial')
            end %end SaveME
            
            NloopRunTime = toc(ticTemp);
            
        end % end bc loop
        
        fprintf('N = %d took %f sec\n',Nx,NloopRunTime)
        % subtrial   = subtrial + 1;  % Just for input file
    end %end Ngr loop
    
    if SaveMe
        SaveStr = sprintf('OPvsGridBCmatsBcNum%dNnum%d',Nbc,Ngr);
        save(SaveStr,'AveNemMatFit', 'AveNemMatPde','StdNemMatPde',...
            'bcVec','NgrVec','L_box','v0', 'SaveMe','trial')
    end
    
    if PlotMe
        AveNemPlotterBCandNgr(AveNemMatFit,AveNemMatPde,StdNemMatPde,...
            bcVec,NgrVec,L_box,v0, SaveMe,trial)
    end
    
    if SaveMe
        
        FileDir  = sprintf('AveNemBcNum%dNnum%dt%d',Nbc,Ngr,trial);
        Where2SavePath    = sprintf('%s/%s/%s',pwd,'Outputs',FileDir);
        
        mkdir(Where2SavePath)
        movefile('*.mat', Where2SavePath)
        movefile('*.txt', Where2SavePath)
        if PlotMe
            movefile('*.fig', Where2SavePath)
        end
        
    end % End Save Me
    
    toc
    
catch ME
    disp(ME.message)
%     exit;
end


