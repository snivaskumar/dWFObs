function [ Wp,strucObs,scriptOptions,LESData,hFigs, d_Wp,d_sol,d_sys,d_strucObs,d_scriptOptions ]...
    = d_WFObs_s_initialize( scriptOptions,configName )
% WFOBS_S_INITIALIZE  Initialize the WFSim model and the estimator settings
%
%   SUMMARY
%    This code does the necessary initializations for the WFSim model, for
%    the estimator, and for the relevant script settings.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - configName: name of the simulation case that is to be simulated.
%     All simulation scenarios can be found in the '/configurations/'
%     folder as seperate files. The default case is 'YawCase3.m'.
%
%     - scriptOptions: this struct contains all simulation settings, not
%     related to the wind farm itself (solution methodology, outputs, etc.)
%
%     - LESData: this struct contains all the flow fields and turbine data
%                from the LES data.
%
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%
%     - sol: this struct contains the system states at a certain timestep.
%         sol.u:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.v:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.uu:    Same as sol.u, used for convergence
%         sol.vv:    Same as sol.v, used for convergence
%
%     - sys: this struct contains the system matrices at a certain timestep.
%         sys.pRCM:  Reverse Cuthill-McKee algorithm for solving A*x=b faster.
%         sys.B1:    Important matrix in the boundary conditions.
%         sys.B2:    Important matrix in the boundary conditions.
%         sys.bc:    Important vector in the boundary conditions.
%
%     - hFigs: cell array of Figures to (re)plot figures into.
%

% Load configuration file from the 'configurations' folder
run(configName);    

% Create destination folder for output files
if (scriptOptions.savePlots + scriptOptions.saveWorkspace > 0)
    mkdir(scriptOptions.savePath);
end;

% Save simulation & filter settings
if (scriptOptions.savePlots + scriptOptions.saveWorkspace > 0)
    save([scriptOptions.savePath '/' strucObs.filtertype '_settings.mat']);
end;

if scriptOptions.printProgress
    disp(' WindFarmObserver (WFObs)');
    disp([' Case:  ' configName ]);
    disp(' ');
end;

% load a default random seed for consistency
if strucObs.loadRandomSeed; load('randomseed'); rng(randomseed); clear randomseed; end;

% Default settings: following WFSim options are never used in WFObs
scriptOptions.Projection      = 0;    % Use projection
scriptOptions.exportLinearSol = 0;    % Export linear solution
scriptOptions.Derivatives     = 0;    % Calculate derivatives/gradients for system

if scriptOptions.printProgress
    disp([datestr(rem(now,1)) ' __  Initializing simulation model.']);
end;

[d_Wp,d_sol,d_sys]      = d_InitWFSim(Wp,scriptOptions); % Initialize model
tur = d_Wp{1}.tur; 
for i = 1:tur
    d_strucObs{i} = strucObs;
    d_scriptOptions{i} = scriptOptions;
    % Add noise to initial conditions
    [d_sol{i}.u,d_sol{i}.uu]  = deal(d_sol{i}.u + randn(d_Wp{i}.mesh.Nx,d_Wp{i}.mesh.Ny)*strucObs.noise_init);
    [d_sol{i}.v,d_sol{i}.vv]  = deal(d_sol{i}.v + randn(d_Wp{i}.mesh.Nx,d_Wp{i}.mesh.Ny)*strucObs.noise_init);
    d_strucObs{i}.size_state = d_Wp{i}.Nu + d_Wp{i}.Nv + d_Wp{i}.Np;
    if d_scriptOptions{i}.exportPressures == 0
        d_strucObs{i}.size_output = d_Wp{i}.Nu + d_Wp{i}.Nv;
    else
        d_strucObs{i}.size_output = d_Wp{i}.Nu + d_Wp{i}.Nv + d_Wp{i}.Np;
    end;
end

% Define measurement locations
if strucObs.measFlow
    sensorsfile        = load(strucObs.sensorsPath);
    strucObs.obs_array = unique([sensorsfile.sensors{1}.obsid; sensorsfile.sensors{2}.obsid]);
    
    % Calculate obs_array locations
    for i = 1:tur
        d_strucObs{i}.obs_array_locu = struct('x',{},'y',{});
        d_strucObs{i}.obs_array_locv = struct('x',{},'y',{});
    end
    for j = 1:length(strucObs.obs_array)
        [ ~,d_locSensor,d_typeFlow,d_strucObs ] = d_WFObs_s_sensors_nr2grid( strucObs.obs_array(j), d_Wp{1}.actualmesh, d_Wp, d_strucObs);
        for i = 1:tur
            if length( d_locSensor{i} )
                if strcmp(d_typeFlow,'u')
                    d_strucObs{i}.obs_array_locu(end+1) = d_locSensor{i};
                else
                    d_strucObs{i}.obs_array_locv(end+1) = d_locSensor{i};
                end
            end
        end
    end
else
    for i = 1:tur
        d_strucObs{i}.obs_array = [];
    end
end;

% Load measurements from LES simulation (*.mat file)
LESData    = load(d_Wp{1}.sim.measurementFile); % Load measurements
LESData.ud = LESData.u + strucObs.noise_obs*randn(size(LESData.u)); % Add noise
LESData.vd = LESData.v + strucObs.noise_obs*randn(size(LESData.v)); % Add noise

for i = 1:tur
% Setup blank figure windows
hFigs{i} = {};
end

% Create global RCM vector
for i = 1:tur
    [~, d_sysRCM{i}] = d_WFSim_timestepping( d_sol{i}, d_sys{i}, d_Wp{i}, d_scriptOptions{i} );
    d_sys{i}.pRCM    = d_sysRCM{i}.pRCM;
    d_scriptOptions{i}.klen = length(num2str(d_Wp{i}.sim.NN));        % used for proper spacing in cmd output window
    d_scriptOptions{i}.tlen = length(num2str(d_Wp{i}.sim.time(end))); % length

    if d_scriptOptions{i}.printProgress
        disp([datestr(rem(now,1)) ' __  Finished initialization sequence.']);
    end;
end
end