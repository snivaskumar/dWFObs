clear all
close all
clc

%ExKF_nofus_turb2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_Queue/axi_2turb_alm_turb_ExKF_nofus_turb2/workspace.mat')

for j = 1:2
    for i = 1:d_Wp{j}.sim.NN
        u05{j}(i) = d_sol_array{j}(i).site.u_Inf;
        RMSE05{j}(i) = d_sol_array{j}(i).score.RMSE_cline;
        maxError05{j}(i) = d_sol_array{j}(i).score.maxError;
        RMSE_flow05{j}(i) = d_sol_array{j}(i).score.RMSE_flow;
    end
    Wp05{j} = d_Wp{j}; sol_array05{j} = d_sol_array{j}; sys05{j} = d_sys{j};
    scriptOptions05{j} = d_scriptOptions{j}; strucObs05{j} = d_strucObs{j};
end
% plotdWFObs( Wp05{1},sol_array05{1},sys05{1},scriptOptions05{1},strucObs05{1} );
% % plotdCOV( Wp05{1},sol_array05{1},sys05{1},scriptOptions05{1},strucObs05{1} );
% plotdWFObs( Wp05{2},sol_array05{2},sys05{2},scriptOptions05{2},strucObs05{2} );

%%

d_x = cell(2,1);
d_x{1} = d_Wp{1}.state.d_x;        d_x{2} = d_Wp{2}.state.d_x;
dx = unique( union(d_x{1},d_x{2}) );
n = length(dx);

P{1} = strucObs05{1}.Pk;
P{2} = strucObs05{2}.Pk;
x{1} = sol_array05{1}(end).x;
x{2} = sol_array05{2}(end).x;
% [x,P] = d_fuze2(x{1},x{2},P{1},P{2},d_x{1},d_x{2},'CIN','OPTIMAL',0.5,10); 
% % [x,P] = fuze(x,P,d_x,2,n,0,0,[1:n]','C',1,'CONSTANT'); % OPTIMAL % CONSTANT
[x{1},P{1}] = d_fuze(x,P,d_x,2,1,0,0,[1:n]','C',1,'CONSTANT');
[x{2},P{2}] = d_fuze(x,P,d_x,2,2,0,0,[1:n]','C',1,'CONSTANT');
for i = 1:2
    strucObs05{i}.Pk = P{i};
    sol_array05{i}(end).x = x{i};
    Nx = Wp05{i}.mesh.Nx;
    Ny = Wp05{i}.mesh.Ny;
    u = 5*ones(Nx,Ny);
    u(3:end-1,2:end-1)= reshape(x{i}(1:(Nx-3)*(Ny-2)),Ny-2,Nx-3)';
    [u,~,~] = Updateboundaries(Nx,Ny,u,u,u);
    sol_array05{i}(end).u = u;
end
plotdWFObs( Wp05{1},sol_array05{1},sys05{1},scriptOptions05{1},strucObs05{1} );
% plotdCOV( Wp05{1},sol_array05{1},sys05{1},scriptOptions05{1},strucObs05{1} );
plotdWFObs( Wp05{2},sol_array05{2},sys05{2},scriptOptions05{2},strucObs05{2} );