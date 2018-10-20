clear all
close all
clc

% Time settings
% k_now_vec  = [10 150  300 750]; % Starting point of forecast
k_now_vec  = [1997]; % Starting point of forecast

% Data sources
data = {};

load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_Queue/turb2axialm_4D_SIM_nofus_turb2_uinf11_noest_Qu0p1Qv0p1R0p1/workspace.mat')
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_Queue/turb2axialm_4D_ExKF_nofus_turb2_uinf11_noest_Qu0p1Qv0p1R0p1/workspace.mat')
path            = '/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_Queue'; 

scriptOptions   = d_scriptOptions{1};
sol_array       = d_sol_array{1};
strucObs        = d_strucObs{1};  
sys             = d_sys{1};
Wp              = d_Wp{1};
save([path '/tmp1.mat'],'strucObs', 'scriptOptions',...
                                'sol_array', 'sys', 'Wp');

scriptOptions   = d_scriptOptions{2};
sol_array       = d_sol_array{2};
strucObs        = d_strucObs{2};  
sys             = d_sys{2};
Wp              = d_Wp{2};
save([path '/tmp2.mat'],'strucObs', 'scriptOptions',...
                                'sol_array', 'sys', 'Wp');
                            
% % load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/turb2axialm_ExKF_uinf11_noest_ICI0p5_turbC_Qu0p1Qv0p1R0p1/workspace.mat')
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/turb2axialm_DExKF_2D_uinf11_noest_ICI0p5_turbC_Qu0p1Qv0p1R0p1/workspace.mat')
% scriptOptions   = d_scriptOptions{2};
% sol_array       = d_sol_array{2};
% strucObs        = d_strucObs{2};  
% sys             = d_sys{2};
% Wp              = d_Wp{2};
% save([path '/tmp.mat'],'strucObs', 'scriptOptions',...
%                                 'sol_array', 'sys', 'Wp');
                            
data{end+1} = struct(...
    'name','Subsystem 1, Without Fusion',...
    'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_Queue/tmp1.mat');
data{end+1} = struct(...
    'name','Subsystem 2, Without Fusion',...
    'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_Queue/tmp2.mat');
% data{end+1} = struct(...
%     'name','',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/tmp.mat');
outputFigName = ['k' strrep(num2str(k_now_vec),' ','_') '_flowContours_SOWFAComparison'];


%% Core operations
addpath('../../WFSim/libraries/export_fig'); % Add export_fig library
addpath('../../dev_tools/ResultAnalysis/libraries'); % Add libraries for VAF calculations

for di = 1:length(data)
    % Load workspace file
    disp(['' num2str(di) '. Loading workspace.mat for ''' data{di}.name '''.']);
    WS{di} = load(data{di}.path);
    
    for kn = 1:length(k_now_vec)
        k_now = k_now_vec(kn);
        sol   = WS{di}.sol_array(k_now);
        
        out(kn,di).u  = sol.u;
        out(kn,di).uq = sol.measuredData.uq;
        out(kn,di).e  = abs(sol.u-sol.measuredData.uq);
    end
end

%% Plot figures
% Meshing settings
% x = WS{1}.Wp.mesh.ldyy;
% y = WS{1}.Wp.mesh.ldxx2;

% Plot velocities in a contourf figure
close all; h = figure; 
if length(k_now_vec) == 3
    h.Position = [448.2000 52.2000 778.4000 726.4000];
end
if length(k_now_vec) == 4
    h.Position = [1.5666e+03 -19.8000 802.4000 967.2000];
end
% set(h,'defaultTextInterpreter','latex')
% set(0,'defaultTextInterpreter','latex')
close all;
climits_u = [0 12];
climits_e = [0 3];
nF_vert = length(data); % number of figures vertically
nF_horz = 3;  % number of figures horizontally

% applied correction for yaw angle: wake was forming at wrong side
for j = 1:length(data)
    x = WS{j}.Wp.mesh.ldyy;
    y = WS{j}.Wp.mesh.ldxx2;
	rotorRotation = -.5*WS{j}.Wp.turbine.Drotor*exp(1i*-WS{j}.Wp.turbine.input(k_now_vec).phi'*pi/180);
        subaxis(nF_vert,nF_horz,3*(j-1)+1,'SpacingVert',0.00,'SpacingHoriz',0.01);
        contourf(x,y,out(1,j).uq,[0:0.1:12],'Linecolor','none');
        caxis(climits_u);
        if j ~= length(data)
            set(gca,'XTickLabel',[]);
        end
        set(gca,'YTick',0:300:max(y(:)));
        ylabel({data{j}.name;'x-direction (m)'},'FontSize',16,'FontWeight','bold')
        if j == length(data)
            xlabel('y-dir. (m)','FontSize',16,'FontWeight','bold')
        end
        axis equal tight;
        if j == 1
            title('SOWFA','FontSize',16,'FontWeight','bold');
        end
%         turb(WS,rotorRotation);

        hold all;
        Wp = WS{j}.Wp;
        for kk=Wp.turbine.ID
            Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
            Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
            rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
                0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
                'FaceColor','w')
            plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
            plot(Qy,Qx,'w','linewidth',2)
        end
        set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
        
        
        subaxis(nF_vert,nF_horz,3*(j-1)+2,'SpacingVert',0.00,'SpacingHoriz',0.01);
        contourf(x,y,out(1,j).u,[0:0.1:12],'Linecolor','none');
        caxis(climits_u);
        if j ~= length(data)
            set(gca,'XTickLabel',[]);
        end
        set(gca,'YTickLabel',[]);
        axis equal tight;
        if j == 1
%             title('EnKF');
%             title('Architecture 4: ExKF','FontSize',16,'FontWeight','bold');
            title('Architecture 4: Open-Loop','FontSize',16,'FontWeight','bold');
        end
        if j == length(data)
            xlabel('y-dir. (m)','FontSize',16,'FontWeight','bold')
        end
%         turb(WS,rotorRotation);

        hold all;
        Wp = WS{j}.Wp;
        for kk = Wp.turbine.ID
                Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
                Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
                rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
                    0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
                    'FaceColor','w')
                plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
                plot(Qy,Qx,'w','linewidth',2)
        end
        set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
        
        if j == length(data)
            xlabel('y-dir. (m)')
            clb_u = colorbar('south');
            clb_u.Position = [0.1625 0.0475 0.4000 0.0100];
            clb_u.Label.String = 'Flow speed (m/s)';
            caxis(climits_u);
        end        
        subaxis(nF_vert,nF_horz,3*(j-1)+3,'SpacingVert',0.00,'SpacingHoriz',0.01);
        contourf(x,y,out(1,j).e,[0:0.1:3],'Linecolor','none');
        caxis(climits_e);
        if j ~= length(data)
            set(gca,'XTickLabel',[]);
        end
        set(gca,'YTickLabel',[]);
        axis equal tight;
        if j == 1
            title('Error','FontSize',16,'FontWeight','bold');
        end
        if j == length(data)
%             xlabel('y-dir. (m)')
%             clb_u = colorbar('south');
%             clb_u.Position = [0.1625 0.0475 0.4000 0.0100];
%             clb_u.Label.String = 'Flow speed (m/s)';
%             caxis(climits_u);
            clb_e = colorbar('south');
            clb_e.Position = [0.6700 0.0475 0.2300 0.0100];
            clb_e.Label.String = 'Error (m/s)';
            clb_e.Limits = climits_e;
        end        

        % Turbines
%         turb(WS,rotorRotation);
        hold all;
        Wp = WS{j}.Wp;
       for kk = Wp.turbine.ID
                Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
                Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
                rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
                    0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
                    'FaceColor','w')
                plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
                plot(Qy,Qx,'w','linewidth',2)
       end
        set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
end
colormap(jet)

% function turb(WS,rotorRotation)
%         hold all;
%         Wp = WS{1}.Wp;
%         for kk=1:Wp.turbine.N
%             Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
%             Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
%             rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
%                 0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
%                 'FaceColor','w')
%             plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
%             plot(Qy,Qx,'w','linewidth',2)
%         end
%         set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
% end