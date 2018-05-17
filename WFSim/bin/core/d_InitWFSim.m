function [Wp,sol,sys] = d_InitWFSim(Wp,options)
%INITWFSIM  Initializes the WFSim model

    % Import simulation scenario (meshing, atmospheric properties, turbine settings)
    [Wp]  = d_meshing(Wp.name,options.plotMesh,options.plotMesh, options.exportPressures);
    tur         = Wp{1}.tur;
    for i = 1:tur
        % Create empty structs
        sys{i} = struct; % This struct will contain all the system matrices at time k
        sol{i} = struct; % This struct will contain the solution (flowfields, power, ...) at time k    

        % Initialize time vector for sol at time k = 0
        sol{i} = struct('k',0,'time',Wp{i}.sim.time(1));

        Nx{i} = Wp{i}.mesh.Nx;
        Ny{i} = Wp{i}.mesh.Ny;
        % Initialize flow fields as uniform ('no turbines present yet')
        [sol{i}.u,sol{i}.uu] = deal( Wp{i}.site.u_Inf *  ones(Nx{i},Ny{i}) );  
        [sol{i}.v,sol{i}.vv] = deal( Wp{i}.site.v_Inf *  ones(Nx{i},Ny{i}) );  
        [sol{i}.p,sol{i}.pp] = deal( Wp{i}.site.p_init * ones(Nx{i},Ny{i}) ); 

        % Initialize the linearized solution variables, if necessary
        if options.Linearversion
            sol{i}.ul = sol{i}.u;
            sol{i}.vl = sol{i}.v;
            sol{i}.pl = sol{i}.p;
            [sol{i}.du,sol{i}.dv,sol{i}.dp]  = deal(zeros(Nx{i},Ny{i}));
        end;    
%         Compute boundary conditions and system matrices B1, B2.
        [sys{i}.B1,sys{i}.B2,sys{i}.bc]  = d_Compute_B1_B2_bc(Wp{i});
        sys{i}.pRCM                = []; % Load empty vector
    end

%     if options.Projection
%         [sys.Qsp, sys.Bsp]  = Solution_space(sys.B1,sys.B2,sys.bc);
%         Wp.Nalpha           = size(sys.Qsp,2);
%     end
end

