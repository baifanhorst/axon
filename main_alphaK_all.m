main_init_more; % Here main_init_more rather than main_init should be used.
addpath('C:\Users\bai_f\Documents\research\Axon_growth\DelayedFeedback\mine_dde23');
%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 
contpar=ind_alphaK;
namepar='alphaK';

alphaK = 0.01;
alphaD = 1;
pE = 6;
wE = 5;
pI = 6;
wI = 5;
vK = 1;
vD = 1;


dK = 0.5;
wK = 1;
dD = 1;
wD = 1;


N = 1;
file_id=1;

min_bound = 0.01;
max_bound = 0.5;
max_step = 0.001;

bf_attempts = 1500;

tspan = [0, 1000];
t_start = 500;
t_end = 1000;
dt = 0.1;


%% Setup the first fixed point
main_setup_first_fixed_pt2;
%% Initialize the branch
main_initialize_branch;
%% Continuation of the branch
main_cont_branch;
%% Plot bifurcation diagram
main_plot_bf;
%% Extract frequency and amplitude of oscillation
main_extract_freq_amp2;
%% output
%save_data_1par(namepar, file_id, N, par_stst, E1_stst, nunst_stst, par_list, freq_list, amp_list, max_list, min_list);










% fprintf('----- Extract frequencies -----\n');
% 
% length_par = length(par_stst);
% freq_list = [];
% max_list = [];
% min_list = [];
% amp_list = [];
% par_list = [];
% T_list =[];
% 
% for i=1:length_par
%     fprintf(append(namepar,sprintf('=%g\n',par_stst(i))));
%     if nunst_stst(i)>0
%         %[alphaK, alphaD, pE, wE, pI, wI, vK, vD, N, dK, wK, dD, wD]
%         switch namepar
%             case 'alphaK'
%                 alphaK = par_stst(i);
%             case 'alphaD'
%                 alphaD = par_stst(i);
%             case 'pE'
%                 pE = par_stst(i);
%             case 'pI'
%                 pI = par_stst(i);
%             case 'wE'
%                 wE = par_stst(i);
%             case 'wI'
%                 wI = par_stst(i);
%             case 'vK'
%                 vK = par_stst(i);
%             case 'vD'
%                 vD = par_stst(i);     
%             case 'N'
%                 N = par_stst(i);
%             case 'dK'
%                 dK = par_stst(i);
%             case 'wK'
%                 wK = par_stst(i);
%             case 'dD'
%                 dD = par_stst(i);
%             case 'wD'
%                 wD = par_stst(i);
%         end
%         
%         gammaK = alpha_to_gamma3(alphaK, dK, wK, vK);
%         rhobulkK = density_bulk(alphaK, gammaK);
%         JK = vK*J(rhobulkK);
%         tauK = delay(rhobulkK, N, vK);
%         
%         gammaD = alpha_to_gamma3(alphaD, dD, wD, vD);
%         rhobulkD = density_bulk(alphaD, gammaD);
%         JD = vD*J(rhobulkD);
%         tauD = delay(rhobulkD, N, vD);
%         
%         
%         delays = [tauK, tauD];
%         
%      
%         JK = vK*J(rhobulkK);
%         JD = vD*J(rhobulkD);
%         tauK = delay(rhobulkK, N, vK);
%         tauD = delay(rhobulkD, N, vD);
%         delays = [tauK, tauD];
%         
%         sol = dde23(@rhs_dde23, delays, @history, tspan);
%         [freq, amp, T, amp_max, amp_min] = extract_freq_amp(sol, 1, t_start, t_end, dt);
%         if amp<1e-2
%             freq = 0;
%         end
%         par_list(end+1) = par_stst(i);
%         freq_list(end+1) = freq;
%         max_list(end+1) = amp_max;
%         min_list(end+1) = amp_min;
%         amp_list(end+1) = amp;
%         T_list(end+1) = T;
%     end
% end
% 
% figure('Name','freq','NumberTitle','off')
% plot(par_list, freq_list, 'b.');
% xlabel(namepar);
% ylabel('Frequency');
% 
% figure('Name','amp','NumberTitle','off')
% plot(par_list, amp_list, 'b.');
% xlabel(namepar);
% ylabel('Amplitude');
% 
% figure('Name','Period','NumberTitle','off')
% plot(par_list, T_list, 'b.');
% xlabel(namepar);
% ylabel('Period');
% 
% 
% figure(fig_bf);
% ax1 = gca;
% hold on;
% plot(ax1,par_list,max_list,'k.', 'DisplayName','max');
% plot(ax1,par_list,min_list,'r.', 'DisplayName','min');
% legend(ax1,'location','north');