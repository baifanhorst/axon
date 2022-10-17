% Codes for generating data for Fig. 6

%% This folder contains codes that solve the system. 
addpath('C:\Users\bai_f\Documents\research\Axon_growth\DelayedFeedback\mine_dde23');

%% Initialization for BIFTOOL
main_init;
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
N = 0.01;

% [file_id] used when naming files to store results.
% Each file_id corresponds to a different N.
% Note that if N is the bifurcation parameter, file_id will not be used.
file_id=1; 
contpar=ind_N; % The index of the continuation parameter.
namepar='N'; % The name of the continuation parameter.

% no_y specifies which signal will be set as y-axis in the bifurcation diagram. 
%E1, E2, I1, I2 -> 1,2,3,4
no_y = 3;
% label_y specify the ylabel in the bifurcation diagram.
label_y = 'I1';



% min max and step for extending 1par branch
min_bound = 0.01;
max_bound = 10;
max_step = 0.01;
% number of attempts when extending 1par branch
bf_attempts = 1000;

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




