function [Wp,sol,sys, d_Wp,d_sol,d_sys] = InitWFSim(Wp,options)
%INITWFSIM  Initializes the WFSim model

    % Import simulation scenario (meshing, atmospheric properties, turbine settings)
%     [Wp, d_Wp]  = meshing(Wp.name,options.plotMesh,options.plotMesh);
%     tur         = d_Wp{1}.tur;
%     for i = 1:tur
%         % Create empty structs
%         d_sys{i} = struct; % This struct will contain all the system matrices at time k
%         d_sol{i} = struct; % This struct will contain the solution (flowfields, power, ...) at time k    
% 
%         % Initialize time vector for sol at time k = 0
%         d_sol{i} = struct('k',0,'time',d_Wp{i}.sim.time(1));
% 
%         d_Nx{i} = d_Wp{i}.mesh.d_Nxe - d_Wp{i}.mesh.d_Nxb + 1;
%         d_Ny{i} = d_Wp{i}.mesh.d_Nye - d_Wp{i}.mesh.d_Nyb + 1;
%         % Initialize flow fields as uniform ('no turbines present yet')
%         [d_sol{i}.u,d_sol{i}.uu] = deal( d_Wp{i}.site.u_Inf *  ones(d_Nx{i},d_Ny{i}) );  
%         [d_sol{i}.v,d_sol{i}.vv] = deal( d_Wp{i}.site.v_Inf *  ones(d_Nx{i},d_Ny{i}) );  
%         [d_sol{i}.p,d_sol{i}.pp] = deal( d_Wp{i}.site.p_init * ones(d_Nx{i},d_Ny{i}) ); 
% 
%         % Initialize the linearized solution variables, if necessary
%         if options.Linearversion
%             d_sol{i}.ul = d_sol{i}.u;
%             d_sol{i}.vl = d_sol{i}.v;
%             d_sol{i}.pl = d_sol{i}.p;
%             [d_sol{i}.du,d_sol{i}.dv,d_sol{i}.dp]  = deal(zeros(d_Nx{i},d_Ny{i}));
%         end;    
% %         Compute boundary conditions and system matrices B1, B2.
%         [d_sys{i}.B1,d_sys{i}.B2,d_sys{i}.bc]  = d_Compute_B1_B2_bc(d_Wp{i});
%         d_sys{i}.pRCM                = []; % Load empty vector
%     end
    
    % Create empty structs
    sys = struct; % This struct will contain all the system matrices at time k
    sol = struct; % This struct will contain the solution (flowfields, power, ...) at time k
    
    % Import simulation scenario (meshing, atmospheric properties, turbine settings)
    [Wp, d_Wp] = meshing(Wp.name,options.plotMesh,options.plotMesh); 

    % Initialize time vector for sol at time k = 0
    sol = struct('k',0,'time',Wp.sim.time(1));
    
    % Initialize flow fields as uniform ('no turbines present yet')
    [sol.u,sol.uu] = deal( Wp.site.u_Inf *  ones(Wp.mesh.Nx,Wp.mesh.Ny) );  
    [sol.v,sol.vv] = deal( Wp.site.v_Inf *  ones(Wp.mesh.Nx,Wp.mesh.Ny) );  
    [sol.p,sol.pp] = deal( Wp.site.p_init * ones(Wp.mesh.Nx,Wp.mesh.Ny) ); 

    % Initialize the linearized solution variables, if necessary
    if options.Linearversion
        sol.ul = sol.u;
        sol.vl = sol.v;
        sol.pl = sol.p;
        [sol.du,sol.dv,sol.dp]  = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny));
    end;

    % Compute boundary conditions and system matrices B1, B2.
    [sys.B1,sys.B2,sys.bc]  = Compute_B1_B2_bc(Wp);
    sys.pRCM                = []; % Load empty vector
    
    % Compute projection matrices Qsp and Bsp. These are only necessary if
    % the continuity equation is projected away (2015 ACC paper, Boersma).
    if options.Projection
        [sys.Qsp, sys.Bsp]  = Solution_space(sys.B1,sys.B2,sys.bc);
        Wp.Nalpha           = size(sys.Qsp,2);
    end
end

