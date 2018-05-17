function [ Wp ] = d_meshing( scenarioName, plotMesh, PrintGridMismatch, exportPressures )
%MESHING Meshing and settings function for the WFSim code
% This code includes all the topology information, atmospheric
% information, turbine properties and turbine control settings for any
% simulation scenario. Basically, all wind farm-related information is
% defined here.

% Default settings
if nargin <= 0; error('Please specify a meshing case.'); end;
if nargin <= 1; plotMesh = true;                         end;
if nargin <= 2; PrintGridMismatch = true;                end;

% Some pre-processing to determine correct input file location
binloc  = which('meshing.m');
if ispc; slashSymbol = '\'; else; slashSymbol = '/'; end; % Compatible with both UNIX and PC
strbslash   = strfind(binloc ,[slashSymbol 'WFSim']);
WFSimfolder = binloc(1:strbslash(end)+6);

switch lower(scenarioName)    
    
    % Wind farms for which PALM data is available        
    case{lower('2turb_adm_turb'),lower('2turb_adm_noturb')}
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = 1.0;        % Turbine power scaling
        forcescale = 1.7;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = true;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice
        lmu        = 0.45;       % Mixing length in x-direction (m)
        mu         = 0*18e-5;    % Dynamic flow viscosity        
        m          = 1;          % Turbulence model gridding property        
        n          = 4;          % Turbulence model gridding property
                         
        % Tuning notes '2turb_adm_turb', '2turb_adm_noturb' (Sep 6th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
        % Note:   Gridsearched multiple ways, should be very reliable now.
        % Results for both cases are so similar that they have been merged.
        
    case lower('2turb_yaw_adm_noturb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = .95;        % Turbine power scaling
        forcescale = 1.6;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = true;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 0.6;        % Mixing length in x-direction (m)
        mu         = 0.0;        % Dynamic flow viscosity
        m          = 1;          % Turbulence model gridding property        
        n          = 4;          % Turbulence model gridding property
        
        % Tuning notes '2turb_yaw_adm_noturb' (Sep 7th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
    
    case lower('6turb_adm_turb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = .95;        % Turbine power scaling
        forcescale = 1.5;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = false;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 0.6;        % Mixing length in x-direction (m)
        mu         = 0.0;        % Dynamic flow viscosity
        m          = 4;          % Turbulence model gridding property        
        n          = 2;          % Turbulence model gridding property        
        
        % Tuning notes '6turb_adm_turb' (Oct 06th, 2017): 
        % Ranges: lmu= xxx, f = xxx, m = xxx, n = xxx
     
    case lower('6turb_adm_partial')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = .95;        % Turbine power scaling
        forcescale = 1.5;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = false;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 0.6;        % Mixing length in x-direction (m)
        mu         = 0.0;        % Dynamic flow viscosity
        m          = 4;          % Turbulence model gridding property        
        n          = 2;          % Turbulence model gridding property        
        
        % Tuning notes '6turb_adm_turb' (Oct 06th, 2017): 
        % Ranges: lmu= xxx, f = xxx, m = xxx, n = xxx
        
    case lower('apc_9turb_adm_noturb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = .90;        % Turbine power scaling
        forcescale = 1.6;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = true;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 0.3;        % Mixing length in x-direction (m)
        mu         = 0.0;        % Dynamic flow viscosity
        m          = 4;          % Turbulence model gridding property        
        n          = 2;          % Turbulence model gridding property        
        
        % Tuning notes 'apc_9turb_adm_noturb' (Sep 10th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
        
    % Wind farms for which SOWFA data is available
    case lower('2turb_alm_noturb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = 0.95;       % Turbine power scaling
        forcescale = 1.40;       % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = true;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 1.7;        % Mixing length in x-direction (m)
        mu         = 0*18e-5;    % Dynamic flow viscosity        
        m          = 1;          % Turbulence model gridding property        
        n          = 4;          % Turbulence model gridding property

        % Tuning notes '2turb_alm_noturb' (Sep 5th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
                  
    case lower('2turb_alm_turb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = 0.95;       % Turbine power scaling
        forcescale = 1.4;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = true;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 1.4;        % Mixing length in x-direction (m)
        mu         = 0*18e-5;    % Dynamic flow viscosity
        m          = 1;          % Turbulence model gridding property        
        n          = 3;          % Turbulence model gridding property

        % Tuning notes '2turb_alm_turb' (Sep 6th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
         
    case lower('apc_9turb_alm_turb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);             % Load the LES meshing file
        Drotor      = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale  = 0.97;        % Turbine power scaling
        forcescale  = 2.0;        % Turbine force scaling
        p_init   = 0.0;           % Initial values for pressure terms (Pa)
        turbul   = true;          % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';    % Turbulence model of choice   
        lmu      = 1.20;          % Mixing length in x-direction (m)
        mu       = 0*18e-5;       % Dynamic flow viscosity        
        m        = 7;             % Turbulence model gridding property          
        n        = 1;             % Turbulence model gridding property
        
        % Tuning notes 'apc_9turb_alm_turb' (Sep 11th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.5, m = 1:8, n = 1:4
        
    case lower('yaw_2turb_alm_noturb')
        error('This case has not yet been tuned/validated. Remove this message manually in meshing.m to continue.');
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);             % Load the LES meshing file
        Drotor      = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale  = 0.95;        % Turbine power scaling
        forcescale  = 1.4;        % Turbine force scaling
        p_init   = 0.0;           % Initial values for pressure terms (Pa)
        turbul   = true;          % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';    % Turbulence model of choice   
        lmu      = 1.0;           % Mixing length in x-direction (m)
        mu       = 0*18e-5;       % Dynamic flow viscosity        
        m        = 2;             % Turbulence model gridding property        
        n        = 1;             % Turbulence model gridding property
                
    case lower('yaw_2turb_alm_turb')
        error('This case has not yet been tuned/validated. Remove this message manually in meshing.m to continue.');
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);             % Load the LES meshing file
        Drotor      = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale  = 1.0;        % Turbine power scaling
        forcescale  = 1.2;        % Turbine force scaling
        p_init   = 0.0;           % Initial values for pressure terms (Pa)
        turbul   = true;          % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';    % Turbulence model of choice   
        lmu      = 1.0;           % Mixing length in x-direction (m)
        mu       = 0*18e-5;       % Dynamic flow viscosity        
        m        = 8;             % Turbulence model gridding property        
        n        = 2;             % Turbulence model gridding property
       
    otherwise
        error('No valid meshing specified. Please take a look at Wp.name.');
end;

% L = 300;

% calculate time
time = 0:h:L;          % time vector for simulation
NN   = length(time)-1; % Total number of simulation steps

% compute input to linear model
for kk=2:length(time)
    turbInput(kk-1).dCT_prime = turbInput(kk).CT_prime - turbInput(kk-1).CT_prime; 
end
turbInput(end).dCT_prime = turbInput(end-1).dCT_prime;

% construct grid
ldx  = linspace(0,Lx,Nx);
ldy  = linspace(0,Ly,Ny);

if strcmp(lower(gridType),'lin')
    % linear gridding
    ldx  = linspace(0,Lx,Nx);
    ldy  = linspace(0,Ly,Ny);
elseif strcmp(lower(gridType),'exp')
    error('Exponential gridding currently not supported.');
%     % exponential gridding
%     ldx  = expvec( Lx, cellSize(1), uniquetol(Crx,1e-2), R_con(1), N_con(1), dx_min(1) );
%     ldy  = expvec( Ly, cellSize(2), uniquetol(Cry,1e-2), R_con(2), N_con(2), dx_min(2) );
%     Nx   = length(ldx);
%     Ny   = length(ldy);
%    
else
    error('Wrong meshing type specified in "type".');
end;
ldx2  = 0.5*(ldx(1:end-1)+ldx(2:end));
ldx2  = [ldx2 2*ldx2(end)-ldx2(end-1)]; % add extra cells
ldy2  = 0.5*(ldy(1:end-1)+ldy(2:end));
ldy2  = [ldy2 2*ldy2(end)-ldy2(end-1)]; % add extra cells

ldxx = repmat(ldx',1,Ny);
ldyy = repmat(ldy,Nx,1);
ldxx2 = repmat(ldx2',1,Ny);
ldyy2 = repmat(ldy2,Nx,1);

tur     = length(Crx);     % Number of turbines
sysLen  = 4*Drotor;
for i = 1:tur
    ID{i}       = find( ( abs(Crx-Crx(i))<=sysLen )&( abs(Cry-Cry(i))<=sysLen ) );
    d_ldy{i}    = ldy( abs(ldy - Cry(i) )<= sysLen );
    Lye{i}      = max( d_ldy{i} );     % Ending Position
    Lyb{i}      = min( d_ldy{i} );     % Beginning Position
    Nyb{i}      = find( ldy == Lyb{i} );
    Nye{i}      = find( ldy == Lye{i} );
    
    d_ldx{i}    = ldx( abs(ldx - Crx(i))<= sysLen );
    Lxe{i}      = max( d_ldx{i} );     % Ending Position
    Lxb{i}      = min( d_ldx{i} );     % Beginning Position
    Nxb{i}      = find( ldx == Lxb{i} );
    Nxe{i}      = find( ldx == Lxe{i} );
        
    d_ldy2{i}   = ldy2( abs(ldy2 - Cry(i))<= sysLen );
    Lsec{i}     = find( ldy2 == min(d_ldy2{i}) );
    Rsec{i}     = find( ldy2 == max(d_ldy2{i}) );
    d_ldx2{i}   = ldx2( abs(ldx2 - Crx(i))<= sysLen );
    Bsec{i}     = find( ldx2 == min(d_ldx2{i}) );
    Tsec{i}     = find( ldx2 == max(d_ldx2{i}) );
    
    d_ldxx{i}   = repmat(d_ldx{i}',1,(Nye{i} - Nyb{i} + 1));
    d_ldyy{i}   = repmat(d_ldy{i},(Nxe{i} - Nxb{i} + 1),1);

    % Create secondary grid from primary grid
    d_ldx2{i}  = 0.5*(d_ldx{i}(1:end-1)+d_ldx{i}(2:end));
    d_ldx2{i}  = [d_ldx2{i} 2*d_ldx2{i}(end)-d_ldx2{i}(end-1)]; % add extra cells
    d_ldy2{i}  = 0.5*(d_ldy{i}(1:end-1)+d_ldy{i}(2:end));
    d_ldy2{i}  = [d_ldy2{i} 2*d_ldy2{i}(end)-d_ldy2{i}(end-1)]; % add extra cells
    d_ldxx2{i} = repmat(d_ldx2{i}',1,(Nye{i} - Nyb{i} + 1));
    d_ldyy2{i} = repmat(d_ldy2{i},(Nxe{i} - Nxb{i} + 1),1);

    % Calculate cell dimensions
    dx{i}   = diff(d_ldx{i});
    dxx{i}  = repmat([dx{i}'; dx{i}(end)],1,(Nye{i} - Nyb{i} + 1));
    dx2{i}  = diff(d_ldx2{i});
    dxx2{i} = repmat([dx2{i}'; dx2{i}(end)],1,(Nye{i} - Nyb{i} + 1));
    dy{i}   = diff(d_ldy{i});
    dyy{i}  = repmat([dy{i}, dy{i}(end)],(Nxe{i} - Nxb{i} + 1),1);
    dy2{i}  = diff(d_ldy2{i});
    dyy2{i} = repmat([dy2{i}, dy2{i}(end)],(Nxe{i} - Nxb{i} + 1),1);

    xline{i} = [];
    kk = 0;
    for j = (ID{i})
        kk = kk + 1;
    % Calculate cells relevant for turbine (x-dir) on primary grid
    [~,tmp] = min( abs(d_ldx{i}-Crx(j)) ); % automatically picks earliest entry in vector
    xline{i} = [xline{i};tmp];
    
    % Calculate cells closest to turbines (y-dir) on both grids
    [ML_prim{i}, L_prim{i} ] = min(abs(d_ldy{i}- (Cry(j)-Drotor/2)));
    [ML_sec{i} , L_sec{i}  ] = min(abs(d_ldy2{i}-(Cry(j)-Drotor/2)));
    [MR_prim{i}, R_prim{i} ] = min(abs(d_ldy{i}- (Cry(j)+Drotor/2)));
    [MR_sec{i} , R_sec{i}  ] = min(abs(d_ldy2{i}-(Cry(j)+Drotor/2)));
    
    yline{i}{kk}  = L_prim{i}:1: R_prim{i}; % turbine cells for primary grid
    % ylinev{i} = L_sec{i} :1: R_sec{i} ; % turbine cells for secondary grid
    ylinev{i}{kk} = L_prim{i}:1: R_prim{i}+1; % JWs code -- Bart: I feel like this needs fixing
    
    if PrintGridMismatch
        % Calculate turbine-grid mismatch
        disp([' TURBINE ' num2str(j) ' GRID MISMATCH:']);
        disp(['                    Primary           Secondary ']);
        disp(['       center:   (' num2str(min(abs(Crx(j)-d_ldx{i})),'%10.2f/n') ',' num2str(min(abs(Cry(j)-d_ldy{i})),'%10.2f/n') ') m.      (' num2str(min(abs(Crx(j)-d_ldx2{i})),'%10.2f/n') ',' num2str(min(abs(Cry(j)-d_ldy2{i})),'%10.2f/n') ') m.']);
        disp(['   left blade:   (' num2str(min(abs(Crx(j)-d_ldx{i})),'%10.2f/n') ',' num2str(ML_prim{i},'%10.2f/n')             ') m.      (' num2str(min(abs(Crx(j)-d_ldx2{i})),'%10.2f/n') ',' num2str(ML_sec{i},'%10.2f/n')                ') m.']);
        disp(['  right blade:   (' num2str(min(abs(Crx(j)-d_ldx{i})),'%10.2f/n') ',' num2str(MR_prim{i},'%10.2f/n')             ') m.      (' num2str(min(abs(Crx(j)-d_ldx2{i})),'%10.2f/n') ',' num2str(MR_sec{i},'%10.2f/n')                ') m.']);
        disp(' ');
    end
    end
    
    % Calculate state sizes
    d_Nx{i} = Nxe{i} - Nxb{i} + 1;
    d_Ny{i} = Nye{i} - Nyb{i} + 1;
    Nu{i} = ( d_Nx{i} - 3 )*( d_Ny{i} - 2 );   % Number of u velocities in state vector
    Nv{i} = ( d_Nx{i} - 2 )*( d_Ny{i} - 3);   % Number of v velocities in state vector
    Np{i} = ( d_Nx{i} - 2 )*( d_Ny{i} - 2 ) - 2; % Number of pressure terms in state vector
    
    if plotMesh
        clf;
        Z1{i} = -2*ones(size(d_ldxx{i}));
        mesh(d_ldyy{i},d_ldxx{i},Z1{i},'FaceAlpha',0,'EdgeColor','black','LineWidth',0.1,'DisplayName','Primary mesh')
        hold on;
        Z2{i} = -1*ones(size(d_ldxx2{i}));
        mesh(d_ldyy2{i},d_ldxx2{i},Z2{i},'FaceAlpha',0,'EdgeColor','blue','LineWidth',0.1,'DisplayName','Secondary mesh')
        hold on;
        for j = (ID{i})
            plot([Cry(j)-Drotor/2;Cry(j)+Drotor/2],[Crx(j);Crx(j)],'LineWidth',3.0,'DisplayName',['Turbine ' num2str(j)])
            hold on
            text(Cry(j),Crx(j),['T ' num2str(j)])
        end;
        axis equal
        xlim([-0.1*Lyb{i} 1.2*Lye{i}]);
        ylim([-0.1*Lxb{i} 1.1*Lxe{i}]);
        xlabel('y (m)');
        ylabel('x (m)');
        view(0,90); % view from the top
        legend('-dynamicLegend');
        drawnow;
        if i ~= tur
            figure;
        end
    end;
end

% States of the un-decomposed system
x_u   = [1:(Nx-3)*(Ny-2)];
u = reshape(x_u,(Ny-2),(Nx-3))';
x_v   = [(Nx-3)*(Ny-2)+1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)];
v = reshape(x_v,(Ny-3),(Nx-2))';
if exportPressures
    x_p = [(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2)];
    p   = reshape(x_p,(Ny-2),(Nx-2))';
    p(isinf(p)) = 0;
end
for k = 1:tur
    k1 = 0;
    for i = Nxb{k}:(Nxe{k}-3)
        k1 = k1+1;
        k2 = 0;
        for j = Nyb{k}:(Nye{k}-2)
            k2 = k2+1;
            d_u{k}(k1,k2) = u(i,j);
        end
    end
    k1 = 0;
    for i = Nxb{k}:(Nxe{k}-2)
        k1 = k1+1;
        k2 = 0;
        for j = Nyb{k}:(Nye{k}-3)
            k2 = k2+1;
            d_v{k}(k1,k2) = v(i,j);
        end
    end
    d_x{k} = [vec( d_u{k}' );...
              vec( d_v{k}' )];
    if exportPressures
        k1 = 0;
        for i = Nxb{k}:(Nxe{k}-2)
            k1 = k1+1;
            k2 = 0;
            for j = Nyb{k}:(Nye{k}-2)
                k2 = k2+1;
                d_p{k}(k1,k2) = p(i,j);
            end
        end
        d_x{k} = [d_x{k}; vec( d_p{k}' )];
        d_x{k} = d_x{k}(1:end-2);
    end
end

%% Export to Wp and input
for i = 1:tur
    N = length(ID{i});
    Wp{i}         = struct('tur',tur,'Nu',Nu{i},'Nv',Nv{i},'Np',Np{i},'name',scenarioName);
    Wp{i}.sim     = struct('h',h,'time',time,'L',L,'NN',NN);
    Wp{i}.turbine = struct('Drotor',Drotor,'powerscale',powerscale,'forcescale',forcescale, ...
        'N',N,'Crx',Crx,'Cry',Cry,'ID',ID{i});
    Wp{i}.turbine.input = turbInput; % System inputs too
    Wp{i}.site    = struct('mu',mu,'Rho',Rho,'u_Inf',u_Inf,'v_Inf',v_Inf,'p_init',p_init, ...
        'lmu',lmu,'turbul',turbul,'m',m,'n',n,'Turbulencemodel',turbModel);
    Wp{i}.mesh    = struct('Lxb',Lxb{i},'Lxe',Lxe{i},'Lyb',Lyb{i},'Lye',Lye{i},...
                    'Nxb',Nxb{i},'Nxe',Nxe{i},'Nyb',Nyb{i},'Nye',Nye{i},...
                    'ldxx',d_ldxx{i},'ldyy',d_ldyy{i},'ldxx2',d_ldxx2{i},'ldyy2',d_ldyy2{i},...
                    'dxx',dxx{i},'dyy',dyy{i},'dxx2',dxx2{i},'dyy2',dyy2{i},'type',gridType);
    Wp{i}.actualmesh = struct('Lx',Lx,'Ly',Ly,'Nx',Nx,'Ny',Ny);
    Wp{i}.actualmesh.ldxx = ldxx;  
    Wp{i}.actualmesh.ldyy = ldyy;
    Wp{i}.actualmesh.ldxx2 = ldxx2;
    Wp{i}.actualmesh.ldyy2 = ldyy2;
    Wp{i}.mesh.xline = xline{i};
    Wp{i}.mesh.yline = yline{i};
    Wp{i}.mesh.ylinev= ylinev{i};
    Wp{i}.mesh.Nx = Wp{i}.mesh.Nxe - Wp{i}.mesh.Nxb + 1;
    Wp{i}.mesh.Ny = Wp{i}.mesh.Nye - Wp{i}.mesh.Nyb + 1;
    Wp{i}.state = struct('u',u,'v',v,'x_u',x_u,'x_v',x_v,...
                         'd_u',d_u{i},'d_v',d_v{i},'d_x',d_x{i});
	if exportPressures
        Wp{i}.state.p = p;
        Wp{i}.state.x_p = x_p;
        Wp{i}.state.d_p = d_p{i};
    end
        
    if exist('measurementFn') == 1
        Wp{i}.sim.measurementFile = measurementFn;
    end
    % Wp{i}.turbine.forcescale           = 5;
    Wp{i}.turbine.actual_forcescale    = forcescale;
    % % Wp{i}.turbine.powerscale           = 5;
    Wp{i}.turbine.actual_powerscale    = powerscale;
    % Wp{i}.site.lmu                   = 5;
    Wp{i}.site.actual_lmu              = lmu;
    % Wp{i}.site.u_Inf                 = 5;
    Wp{i}.site.actual_u_Inf            = u_Inf;
    % Wp{i}.site.v_Inf                   = 5;
    Wp{i}.site.actual_v_Inf            = v_Inf;

end
%% Construct mu if no turbulence
% if turbul==0; mu = ConstructMu(Wp); Wp.site(:).mu = mu; end

end