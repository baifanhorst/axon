main_init;
addpath('C:\Users\bai_f\Documents\research\Axon_growth\DelayedFeedback\mine_dde23');
%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 
rhobulkK = 0.5;
rhobulkD = 0.5;
pE = 6;
wE = 1;
pI = 6;
wI = 5;
vK = 1;
vD = 1;




% N=1-->file_id=1 
% N=2-->file_id=0.11 
% N=3-->file_id=10 (this can be different for different parameters)
N = 10;
file_id=3;
% index of the continuation parameter
contpar=ind_wE;
% name of the continuation parameter
namepar='wE';

% min max and step for extending 1par branch
min_bound = 1;
max_bound = 20;
max_step = 0.01;
% number of attempts when extending 1par branch
bf_attempts = 1000;

% Parameters needed for 2par bifurcation.
contpar_2D_1 = ind_wE;
contpar_2D_2 = ind_wI;
max_bound_2D_1 = 20;
max_bound_2D_2 = 20;
min_bound_2D_1 = 1;
min_bound_2D_2 = 1;
max_step_2D_1 = 0.01;
max_step_2D_2 = 0.01;

namepar_2D_1 = 'wE';
namepar_2D_2 = 'wI';

% This is used when continuing the 2par diagram.
% The sign shows the direction.
step_direction = 0.01;

% Figure name (name of the window)
figure_name = sprintf("wE wI N=%g", N);
% Title name
title_name = sprintf("wE wI N=%g", N);
% Number of attempts when continuing 2par branch. 
num_attempts = 2000;
% x and y range for the figure
x_range = [0,20];
y_range = [0,20];

%% Setup the first fixed point
main_setup_first_fixed_pt;

%% Initialize the branch
main_initialize_branch;

%% Continuation of the branch
main_cont_branch;

%% Plot bifurcation diagram
% No need to plot again. This can be done in 1-par code, such as
% main_rhobulkK_all.m

%% Continuation of Hopf bifurcation in two parameters
main_cont_branch_2par;



% %% Setup the first fixed point
% stst.kind = 'stst';
% stst.parameter = [rhobulkK, rhobulkD, pE, wE, pI, wI, vK, vD, N];
% % This values are approximations obtained from dde23.
% JK = vK*J(rhobulkK);
% JD = vD*J(rhobulkD);
% tauK = delay(rhobulkK, N, vK);
% tauD = delay(rhobulkD, N, vD);
% delays = [tauK, tauD];
% tspan = [0, 200];
% sol = dde23(@rhs_dde23, delays, @history, tspan);
% 
% E1_st=sol.y(1,end);
% E2_st=sol.y(2,end);
% I1_st=sol.y(3,end);
% I2_st=sol.y(4,end);
% 
% 
% stst.x = [E1_st; E2_st; I1_st; I2_st];
% % Correct the approximations above.
% method = df_mthod(funcs, 'stst');
% % If there is an error message caused by running the following line,
% % you should check set_funcs to see whether the parameter positions are
% % correct or not.
% [stst, success] = p_correc(funcs, stst, [], [], method.point);
% E1_st = stst.x(1);
% E2_st = stst.x(2);
% I1_st = stst.x(3);
% I2_st = stst.x(4);
% 
% contpar=ind_wE;
% namepar='wE';
% 
% %% Initialize the branch
% stst_branch0 = SetupStst(funcs,'x',[E1_st; E2_st; I1_st; I2_st],'parameter',stst.parameter,...
%     'contpar',contpar,'max_step',[contpar,0.1],'min_bound',...
%     [contpar 1],'max_bound',[contpar 20],...
%     'newheuristics_tests',0);
% 
% %% Continuation of the branch
% figure('Name',append('Generating branch for ',namepar),'NumberTitle','off');
% clf;
% ax1=gca;
% title(ax1,sprintf(append('Generating branch for ',namepar)));
% [stst_branch0] = br_contn(funcs,stst_branch0,500);
% [stst_branch_wbifs,stst_testfuncs]=LocateSpecialPoints(funcs,stst_branch0);
% nunst_stst=GetStability(stst_branch_wbifs);
% 
% %% Plot bifurcation diagram
% par_stst=getpar(stst_branch_wbifs,contpar);
% E1_stst=getx(stst_branch_wbifs,1);
% fig_bf = figure('Name',append('Bifurcation diagram for ',namepar),'NumberTitle','off');
% clf;
% ax1=gca;
% hold on;
% plot(ax1,par_stst(nunst_stst==0),E1_stst(nunst_stst==0),'g.','DisplayName','stable');
% plot(ax1,par_stst(nunst_stst>0),E1_stst(nunst_stst>0),'b.','DisplayName','unstable');
% plot(ax1,bgetpar(stst_branch_wbifs,contpar,'hopf'),bgetx(stst_branch_wbifs,1,'hopf'),'ks','DisplayName','hopf');
% 
% legend(ax1,'location','best');
% xlabel(ax1,namepar);
% ylabel(ax1,'E1');
% title(ax1,sprintf(append('Bifurcation diagram for ',namepar)));
% 
% %% Continuation of Hopf bifurcation in two parameters
% fprintf('----- Hopf branch -----\n');
% 
% ind_par1 = ind_wE;
% ind_par2 = ind_wI;
% namepar1 = 'wE';
% namepar2 = 'wI';
% 
% parameter_bd={'max_bound',[ind_par1, 20; ind_par2, 20],...
%     'min_bound',[ind_par1, 1; ind_par2, 1],...
%     'max_step',[0,0.1; ind_par1,0.1; ind_par2,0.1]};
% 
% % Point numbers of the two hopf bifurcation points.
% num_hopf = br_getflags(stst_branch_wbifs,'hopf');
% % It seems that 'SetupHopf' sets hopf_branch0 to be a 2-par branch, so that
% % 'br_contn' knows to continue it this way, rather than as a steady-state
% % branch.
% [hopf_branch0,suc] = SetupHopf(funcs, stst_branch_wbifs, num_hopf(1),...
%     'contpar', [ind_par1,ind_par2],...
%     'dir', ind_par1, 'step', 0.1, parameter_bd{:});
% 
% figure('Name', sprintf('N=%g', N), 'NumberTitle','off');clf;
% ax2=gca;
% title(ax2, sprintf('N=%g', N));
% % xlim(ax2,[0,1]);
% % ylim(ax2,[0,1]);
% 
% hopf_branch0=br_contn(funcs,hopf_branch0,700);
% 
% xlabel(ax2, namepar1);
% ylabel(ax2, namepar2);



