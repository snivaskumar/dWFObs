function [ outputData ] = d_WFObs_core( scriptOptions, configName, WpOverwrite )
% WFOBS_CORE  Perform a complete time simulation with state estimation
%
%   SUMMARY
%    This code will complete a full wind farm simulation including state
%    estimation, as set up in configurations/*configName*.m.  It will use
%    measurements from data_SOWFA/* or data_PALM/* to improve the flow
%    estimations compared to the WFSim model.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - configName: name of the simulation case that is to be simulated.
%     All simulation scenarios can be found in the '/configurations/'
%     folder as seperate files. The default case is 'YawCase3.m'.
%
%     - scriptOptions: this struct contains all simulation settings, not
%     related to the wind farm itself (solution methodology, outputs, etc.)
%
%     - OutputData (*): this struct contains all interesting outputs. It packs
%       several substructs, namely:
%       - *.Wp: this struct contains all the simulation settings related to
%         the wind farm, the turbine inputs, the atmospheric properties, etc.
%         See WFSim.m for a more elaborate description of 'Wp'.
%
%       - *.sol_array: this cell array contains the system states at every
%         simulated time instant. Each cell entry contains a sol struct.
%         See WFSim.m for a more elaborate description of 'sol'. In
%         addition to the entries described in WFSim.m, each 'sol' struct
%         contains in addition:
%           *.sol.score: a struct containing estimation performance scores
%           such as the maximum estimation error, the RMS error, and the
%           computational cost (CPU time).
%           *.sol.measuredData: a struct containing the true (to be
%           estimated) values, and the measurement data fed into the
%           estimation algorithm.
%
%       - *.strucObs: this struct contains all the observer settings and
%       (temporary) files used for updates, such as covariance matrices,
%       ensemble/sigma point sets, measurement noise, etc.
%

%% Pre-processing
timerScript = tic;       % Start script timer

% Initialize model and observer variables
% [Wp,sol,sys,strucObs,scriptOptions,LESData,hFigs] = ...
%     WFObs_s_initialize(scriptOptions,configName);
[Wp,strucObs,scriptOptions,LESData,d_hFigs, d_Wp,d_sol,d_sys,d_strucObs,d_scriptOptions] = ...
    d_WFObs_s_initialize(scriptOptions,configName);

max_it   = scriptOptions.max_it;    % Convergence constraints
conv_eps = scriptOptions.conv_eps;  % Convergence constraints

tur = d_Wp{1}.tur;
% Overwrite variables if WpOverwrite is specified
if nargin > 2
    if scriptOptions.printProgress
        disp('Overwriting variables in Wp...');
    end
%     Wp = mergeStruct(Wp,WpOverwrite);
%     [sys.B1,sys.B2,sys.bc] = Compute_B1_B2_bc(Wp); % Update boundary conditions
%     tur = d_Wp{1}.tur;
    for i = 1:tur
        d_Wp{i} = mergeStruct(d_Wp{i},WpOverwrite);
        [d_sys{i}.B1,d_sys{i}.B2,d_sys{i}.bc]  = d_Compute_B1_B2_bc(d_Wp{i});
    end
end

%% Core: time domain simulations
d_sol_array = cell(tur,1);
while d_sol{1}.k < d_Wp{1}.sim.NN
    timerCPU = tic;                 % Start iteration timer
%     d_sol_array = cell(tur,1);
    for i = 1:tur
        d_sol{i}.k    = d_sol{i}.k + 1;           % Timestep forward
        d_sol{i}.time = d_Wp{i}.sim.time(d_sol{i}.k+1);% Timestep forward

        % Load measurement data
        d_sol{i}.measuredData = d_WFObs_s_loadmeasurements(LESData,d_sol{i}.k,d_Wp{i});

        % Determine freestream inflow properties from SCADA data
    %     [ Wp,sol,sys,strucObs ] = WFObs_s_freestream(Wp,sol,sys,strucObs);

        % Calculate optimal solution according to filter of choice
        [d_Wp{i},d_sol{i},d_strucObs{i}] = d_WFObs_o(d_strucObs{i},d_Wp{i},d_sys{i},d_sol{i},d_scriptOptions{i});
        
%         if i == 1
%             Nxb = d_Wp{2}.mesh.NNxb;
%             Nxe = d_Wp{2}.mesh.NNxe;
%             Nyb = d_Wp{2}.mesh.NNyb;
%             Nye = d_Wp{2}.mesh.NNye;
%             [ru,cu] = size(d_sol{2}.u);
%             [rv,cv] = size(d_sol{2}.v);
%             tmpu = d_sol{2}.u;        tmpuu = d_sol{2}.uu;
%             tmpv = d_sol{2}.v;        tmpvv = d_sol{2}.vv;
%             u = d_sol{1}.u;     uu = d_sol{1}.uu;
%             v = d_sol{1}.v;     vv = d_sol{1}.vv;
%             for j = 1:ru
%                 for k = 1:cu
%                     if (( Nxb{1}+j-1 )<=Nxe{1})&&(( Nyb{1}+k-1 )<=Nye{1})
%                         tmpu(j,k) = u( Nxb{1}+j-1,Nyb{1}+k-1 );
%                         tmpuu(j,k) = uu( Nxb{1}+j-1,Nyb{1}+k-1 );  
%                     end
%                 end
%             end
%             for j = 1:rv
%                 for k = 1:cv
%                     if (( Nxb{1}+j-1 )<=Nxe{1})&&(( Nyb{1}+k-1 )<=Nye{1})
%                         tmpv(j,k) = v( Nxb{1}+j-1,Nyb{1}+k-1 );
%                         tmpvv(j,k) = vv( Nxb{1}+j-1,Nyb{1}+k-1 );  
%                     end
%                 end
%             end
%             d_sol{2}.u = tmpu;        d_sol{2}.uu = tmpuu;
%             d_sol{2}.v = tmpv;        d_sol{2}.vv = tmpvv;
%         end
    end
    filter          = strucObs.filtertype; 
    fusion          = upper( scriptOptions.fusion );
    fusion_type     = upper( strucObs.fusionDomain );
    fusion_weight   = upper( strucObs.fusion_weight );
    constant        = strucObs.fusion_CIconstant; 
    iteration       = strucObs.fusion_CIiteration; 
    
    if ( strcmp(strucObs.filtertype,'dexkf')||...
            strcmp(strucObs.filtertype,'exkf') )&&...
            strcmp(fusion,'YES')
        d_x{1} = d_Wp{1}.state.d_x;        d_x{2} = d_Wp{2}.state.d_x;
        x = unique( union(d_x{1},d_x{2}) );
        n = length(x);
        z{1} = d_sol{1}.x;      Z{1} = d_strucObs{1}.Pk;
        z{2} = d_sol{2}.x;      Z{2} = d_strucObs{2}.Pk;
        if strcmp(fusion_type,'IFAC') 
            [xe,Ce] = dfuze(z,Z,d_x,tur,n,0,0,[1:n]','C','CONSTANT');
            tmp1 = ismember(x,d_x{1});
            tmp2 = ismember(x,d_x{2});
            d_sol{1}.x = xe(tmp1);    d_strucObs{1}.Pk = Ce(tmp1,tmp1);
            d_sol{2}.x = xe(tmp2);    d_strucObs{2}.Pk = Ce(tmp2,tmp2);
        elseif strcmp(fusion_type,'D_IFAC') 
            [d_sol{1}.x,d_strucObs{1}.Pk] = d_fuze(z,Z,d_x,tur,1,0,0,[1:n]','C',1,fusion_weight);
            [d_sol{2}.x,d_strucObs{2}.Pk] = d_fuze(z,Z,d_x,tur,2,0,0,[1:n]','C',1,fusion_weight);
        elseif strcmp(fusion_type,'AVG') 
            zz=z;
            l1 = length(d_x{1});
            l2 = length(d_x{2});
            tmp1 = zeros(l2,1);    tmp2 = zeros(l1,1);
            for i = 1:l1
                loc = d_x{1}(i)==d_x{2};
                if sum( double( loc ) )
                    tmp1 = tmp1 + double( loc );
                    tmpz = z{2}( loc );
                    tmpZ = Z{2}( loc,loc );
                    zz{1}(i) = ( z{1}(i)+z{2}(loc) )/2;
                else
                    zz{1}(i) = z{1}(i);
                end
            end
            for i = 1:l2
                if sum( d_x{2}(i)==d_x{1} )
                    tmp2 = tmp2 + double( d_x{2}(i)==d_x{1} );
                    zz{2}(i) = ( z{2}(i)+z{1}(d_x{2}(i)==d_x{1}) )/2;
                else
                    zz{2}(i) = z{2}(i);
                end
            end
            z = zz;
            d_sol{1}.x = z{1};  d_sol{2}.x = z{2};
        elseif strcmp(fusion_type,'BAR')
            a1 = d_x{1};
            a2 = d_x{2};
            l1 = length(a1);
            l2 = length(a2);
            
            k1 = 0;
            H1 = [];
            for iii = 1:l1
                if sum(a1(iii)==a2)
                    k1 = k1+1;
                    H1(k1,:) = double(a1(iii)==a2);
                end
            end
            H1 = sparse(H1);
            k2 = 0;
            H2 = [];
            for iii = 1:l2
                if sum(a2(iii)==a1)   
                    k2 = k2+1;
                    H2(k2,:) = double(a2(iii)==a1);
                end
            end
            H2 = sparse(H2);

%             Z{1} = 2*Z{1};
%             Z{2} = 2*Z{2};
            
            zftmp = z;
            Pftmp = Z;
            Pff = pinv( H2*Z{1}*H2' + H1*Z{2}*H1' );
            K1 = Z{1}*H2'*Pff;
            K2 = Z{2}*H1'*Pff;
            zftmp{1} = zftmp{1} + K1*(H1*z{2}-H2*z{1});
            zftmp{2} = zftmp{2} + K2*(H2*z{1}-H1*z{2});
            ll1 = length(K1*H2);
            Pftmp{1} = ( eye(ll1,ll1)-K1*H2 )*Z{1} ;
            ll2 = length(K2*H1);
            Pftmp{2} = ( eye(ll2,ll2)-K2*H1 )*Z{2} ;
            z = zftmp;
            Z = Pftmp;
            
            d_sol{1}.x = z{1};    d_strucObs{1}.Pk = Z{1};
            d_sol{2}.x = z{2};    d_strucObs{2}.Pk = Z{2};
        else
%             [zf, Zf] = d_fuze2(z{1},z{2},Z{1},Z{2},d_x{1},d_x{2},fusion_type,fusion_weight,constant,iteration);
%             d_sol{1}.x = zf{1};         d_strucObs{1}.Pk = Zf{1};
%             d_sol{2}.x = zf{2};         d_strucObs{2}.Pk = Zf{2};
            
            [zf, Zf] = fuse2(z{1},z{2},Z{1},Z{2},d_x{1},d_x{2},fusion_type,fusion_weight,constant,iteration,filter,d_sol{1}.k);
            tmp1 = ismember(x,d_x{1});
            tmp2 = ismember(x,d_x{2});
            d_sol{1}.x = zf(tmp1);    d_strucObs{1}.Pk = Zf(tmp1,tmp1);
            d_sol{2}.x = zf(tmp2);    d_strucObs{2}.Pk = Zf(tmp2,tmp2);
        end
    end
    for i = 1:tur
        [d_sol{i},~]  = MapSolution(d_Wp{i},d_sol{i},Inf,d_scriptOptions{i}); % Map solution to flowfields
        [~,d_sol{i}]  = d_Actuator(d_Wp{i},d_sol{i},d_scriptOptions{i});        % Recalculate power after analysis update

        % Display progress in the command window
        d_sol{i} = WFObs_s_reporting(timerCPU,d_Wp{i},d_sol{i},d_strucObs{i},d_scriptOptions{i});

        % Save reduced-size solution to an array
        d_sol{i}.measuredData = rmfield(d_sol{i}.measuredData,{'u','v','sol'});
        d_sol{i}.site         = d_Wp{i}.site; % Save site info too
        if nnz(strcmp(fieldnames(d_sol{i}),'uk')) >= 1
            d_sol_array{i}(d_sol{i}.k) = rmfield(d_sol{i},{'uu','vv','pp','uk','vk'});
        else
            d_sol_array{i}(d_sol{i}.k) = rmfield(d_sol{i},{'uu','vv','pp'});
        end
        % Display animations on screen
        [d_hFigs{i},d_scriptOptions{i}] = d_WFObs_s_animations(d_Wp{i},d_sol_array{i},d_sys{i},LESData,d_scriptOptions{i},d_strucObs{i},d_hFigs{i});
    end
end

%% Post-processing
% save workspace variables, if necessary
if scriptOptions.saveWorkspace
    save([scriptOptions.savePath '/workspace.mat'],'configName',...
        'd_Wp','d_sys','d_sol_array','d_scriptOptions','d_strucObs');
end

% Put all relevant outputs in one structure
if nargout >= 1
    outputData.d_sol_array     = d_sol_array;
    outputData.d_Wp            = d_Wp;
    outputData.d_strucObs      = d_strucObs;
    outputData.d_scriptOptions = d_scriptOptions;
    outputData.configName    = configName;
end;

% Print end of simulation statement
if scriptOptions.printProgress
    disp([datestr(rem(now,1)) ' __  Completed simulations. ' ...
        'Total CPU time: ' num2str(toc(timerScript)) ' s.']);
end
end