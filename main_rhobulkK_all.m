% Codes for generating data for Fig. 8AB
main_init;
addpath('C:\Users\bai_f\Documents\research\Axon_growth\DelayedFeedback\mine_dde23');
%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 
rhobulkK = 0.01;
rhobulkD = 0.5;
pE = 6;
wE = 5;
pI = 6;
wI = 5;
vK = 1; % vK=2; % This is for different motor speeds. 
vD = 1;

% N=1-->file_id=1 
% N=2-->file_id=0.11 
% N=3-->file_id=10 (this can be different for different parameters)

N = 10;
file_id=4;
contpar=ind_rhobulkK;
namepar='rhobulkK';

% no_y specifies which signal will be set as y-axis in the bifurcation diagram. 
%E1, E2, I1, I2 -> 1,2,3,4
no_y = 3;
% label_y specify the ylabel in the bifurcation diagram.
label_y = 'I1';


% min max and step for extending 1par branch
min_bound = 0.02;
max_bound = 0.98;
max_step = 0.0005;  % max_step = 0.0001; %This is for different motor speed %May decrease the step size, otherwise error.
% number of attempts when extending 1par branch
bf_attempts = 2000; % bf_attempts = 15000; %This is for different motor speed %May increase #steps, otherwise error

% time instants for extracting frequencies
t_start = 500;
t_end = 1000;
tspan = [0, t_end];
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
main_save_data_1par;

%save_data_1par(namepar, file_id, N, par_stst, y_stst, label_y, nunst_stst, par_list, freq_list, T_list, amp_list, max_list, min_list, par_bf_hopf, y_bf_hopf);

%save_data_1par(namepar, file_id, N, par_stst, E1_stst, nunst_stst, par_list, freq_list, amp_list, max_list, min_list);












% %% Setup the first fixed point
% stst.kind = 'stst';
% stst.parameter = [rhobulkK, rhobulkD, pE, wE, pI, wI, vK, vD, N];
% % This values are approximations obtained from dde23.
% E1_st=4.84493267843259;
% E2_st=1.15067126224841;
% I1_st=0.329106801770075;
% I2_st=0.263338159420767;
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
% contpar=ind_rhobulkK;
% namepar='rhobulkK';
% 
% %% Initialize the branch
% stst_branch0 = SetupStst(funcs,'x',[E1_st; E2_st; I1_st; I2_st],'parameter',stst.parameter,...
%     'contpar',contpar,'max_step',[contpar,0.001],'min_bound',...
%     [contpar 0.02],'max_bound',[contpar 0.98],...
%     'newheuristics_tests',0);
% 
% %% Continuation of the branch
% figure('Name',append('Generating branch for ',namepar),'NumberTitle','off');
% clf;
% ax1=gca;
% title(ax1,sprintf(append('Generating branch for ',namepar)));
% [stst_branch0] = br_contn(funcs,stst_branch0,2000);
% % Sometimes there is an error about the following command:
% % 'Dot indexing is not supported for variables of this type.'
% % In this case, we can set the max_step in the above block smaller.
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
% legend(ax1,'location','north');
% xlabel(ax1,namepar);
% ylabel(ax1,'E1');
% title(ax1,sprintf(append('Bifurcation diagram for ',namepar)));
% 
% %% Extract frequency and amplitude of oscillation
% 
% JD = vD*J(rhobulkD);
% 
% tspan = [0, 400];
% t_start = 100;
% t_end = 400;
% dt = 0.1;
% 
% 
% 
% length_par = length(par_stst);
% freq_list = [];
% max_list = [];
% min_list = [];
% amp_list = [];
% par_list = [];
% 
% for i=1:length_par
%     if nunst_stst(i)>0
%         rhobulkK = par_stst(i)
%         JK = vK*J(rhobulkK);
%         tauK = delay(rhobulkK, N, vK);
%         tauD = delay(rhobulkD, N, vD);
%         delays = [tauK, tauD];
%         sol = dde23(@rhs_dde23, delays, @history, tspan);
%         [freq, amp, T, amp_max, amp_min] = extract_freq_amp(sol, 1, t_start, t_end, dt);
%         if amp<1e-2
%             freq = 0;
%         end
%         par_list(end+1) = rhobulkK;
%         freq_list(end+1) = freq;
%         max_list(end+1) = amp_max;
%         min_list(end+1) = amp_min;
%         amp_list(end+1) = amp;
%     end
% end
% 
% figure('Name','freq','NumberTitle','off')
% plot(par_list, freq_list, 'b.');
% xlabel('\rho_{bulk,K}');
% ylabel('Frequency');
% 
% figure('Name','amp','NumberTitle','off')
% plot(par_list, amp_list, 'b.');
% xlabel('\rho_{bulk,K}');
% ylabel('Amplitude');
% 
% 
% figure(fig_bf);
% ax1 = gca;
% hold on;
% plot(ax1,par_list,max_list,'k.', 'DisplayName','max');
% plot(ax1,par_list,min_list,'r.', 'DisplayName','min');
% legend(ax1,'location','north');
% 
% %% output
% save_data_1par(namepar, file_id, N, par_stst, E1_stst, nunst_stst, par_list, freq_list, amp_list, max_list, min_list);
