main_init;
addpath('C:\Users\bai_f\Documents\research\Axon_growth\DelayedFeedback\mine_dde23');
%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 
rhobulkK = 0.5;
rhobulkD = 0.5;
pE = 6;
wE = 5;
pI = 6;
wI = 5;
vK = 1;
vD = 1;

% N=1-->file_id=1 
% N=2-->file_id=0.11 
% N=3-->file_id=10 (this can be different for different parameters)
N = 0.01;
file_id=1;
% index of the continuation parameter
contpar=ind_N;
% name of the continuation parameter
namepar='N';

% min max and step for extending 1par branch
min_bound = 0.01;
max_bound = 1;
max_step = 0.01;
% number of attempts when extending 1par branch
bf_attempts = 1000;


%% Setup the first fixed point
main_setup_first_fixed_pt;

%% Initialize the branch
main_initialize_branch;

%% Continuation of the branch
main_cont_branch;

% save("main_2par_N_parameter_part1")










