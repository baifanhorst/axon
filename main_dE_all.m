main_init;
addpath('C:\Users\bai_f\Documents\research\Axon_growth\DelayedFeedback\mine_dde23');
%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 


rhobulkK = 0.5;
rhobulkD = 0.5;
pE = 1;
wE = 5;
pI = 6;
wI = 5;
vK = 1;
vD = 1;

contpar=ind_pE;
namepar='pE';

% N=1-->file_id=1 
% N=2-->file_id=0.11 
% N=3-->file_id=10 
N = 10;
file_id=3; % N=1-->file_id=1 

% no_y specifies which signal will be set as y-axis in the bifurcation diagram. 
%E1, E2, I1, I2 -> 1,2,3,4
no_y = 3;
% label_y specify the ylabel in the bifurcation diagram.
label_y = 'I1';


% min max and step for extending 1par branch
min_bound = 1;
max_bound = 50;
max_step = 0.1;
% number of attempts when extending 1par branch
bf_attempts = 800;

% time instants for extracting frequencies
tspan = [0, 1000];
t_start = 500;
t_end = 1000;
dt = 0.1;

%% Setup the first fixed point
main_setup_first_fixed_pt;
%% Initialize the branch
main_initialize_branch;
%% Continuation of the branch
main_cont_branch;
%% Plot bifurcation diagram
main_plot_bf;
%% Extract frequency and amplitude of oscillation
main_extract_freq_amp;
%% output
save_data_1par(namepar, file_id, N, par_stst, y_stst, label_y, nunst_stst, par_list, freq_list, T_list, amp_list, max_list, min_list, par_bf_hopf, y_bf_hopf);


%save_data_1par(namepar, file_id, N, par_stst, E1_stst, nunst_stst, par_list, freq_list, amp_list, max_list, min_list);

% %% 2-par bifurcation
% fprintf('----- Hopf branch -----\n');
% 
% ind_par1 = ind_pE;
% ind_par2 = ind_pI;
% namepar1 = 'pE';
% namepar2 = 'pI';
% 
% parameter_bd={'max_bound',[ind_par1, 50; ind_par2, 50],...
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
% fig_2par = figure('Name', sprintf('N=%g', N), 'NumberTitle','off');clf;
% ax2=gca;
% title(ax2, sprintf('N=%g', N));
% % xlim(ax2,[0,1]);
% % ylim(ax2,[0,1]);
% 
% hopf_branch0=br_contn(funcs,hopf_branch0,700);
% 
% xlabel(ax2, namepar1);
% ylabel(ax2, namepar2);
% 
% 
% 
% %% Periodic branch
% hopf=hopf_branch0.point(1);
% [psol,stp]=p_topsol(funcs,hopf,0.01,3,10);
% mpsol=df_mthod(funcs,'psol');
% [psol,s]=p_correc(funcs,psol,contpar,stp,mpsol.point);
% psol_branch=df_brnch(funcs,contpar,'psol');
% psol_branch.point=psol;
% [psol,stp]=p_topsol(funcs,hopf,0.02,3,10);
% [psol,s]=p_correc(funcs,psol,contpar,stp,mpsol.point);
% psol_branch.point(2)=psol;
% 
% figure('Name', 'Periodic branch', 'NumberTitle','off');clf;
% [xm,ym]=df_measr(0,psol_branch);
% ym.field='period';
% ym.col=1;
% ym.row=1;
% psol_branch.method.continuation.plot_measure.x=xm;
% psol_branch.method.continuation.plot_measure.y=ym;
% % Continuation of the periodic branch.
% [psol_branch,s,r,f]=br_contn(funcs,psol_branch,100);
% xlabel(namepar);ylabel('period');
% 













% %% Periodic branch
% % Get the first Hopf point
% % Get its index first.
% for i=1:length(nunst_stst)
%     if nunst_stst(i)==2
%         ind_firstHopf = i;
%         break;
%     end
% end
% 
% 
% hopf=stst_branch_wbifs.point(ind_firstHopf);
% 
% % Convert the point to a Hopf point.
% hopf=p_tohopf(funcs, hopf);
% method=df_mthod(funcs,'hopf');
% [hopf, success] = p_correc(funcs, hopf, contpar, [], method.point);
% 
% 
% 
% % Estimate the amplitude
% 
% pE = hopf.parameter(contpar);
% JK = vK*J(rhobulkK);
% JD = vD*J(rhobulkD);
% tauK = delay(rhobulkK, N, vK);
% tauD = delay(rhobulkD, N, vD);
% delays = [tauK, tauD];
%         
%         
% sol = dde23(@rhs_dde23, delays, @history, tspan);
% [freq, amp, T, amp_max, amp_min] = extract_freq_amp(sol, 1, t_start, t_end, dt);
% 
% [psol,stp]=p_topsol(funcs,hopf,amp,3,ind_firstHopf); %???amp 到底该用什么？
% [psol,s]=p_correc(funcs,psol,contpar,stp,method.point);

% psol_branch=df_brnch(funcs,contpar,'psol');
% psol_branch.point=psol;
% 
% [psol,stp]=p_topsol(funcs,hopf,amp*1.1,3,ind_firstHopf);
% [psol,s]=p_correc(funcs,psol,contpar,stp,method.point);
% psol_branch.point(2)=psol;
% 
% 
% [xm,ym]=df_measr(0,psol_branch);
% ym.field='period';
% % ?????????
% ym.col=1;
% ym.row=1;
