main_init;
addpath('C:\Users\bai_f\Documents\research\Axon_growth\DelayedFeedback\mine_dde23');
%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 
rhobulkK = 0.02;
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
N = 1;
file_id=1;
% index of the continuation parameter
contpar=ind_rhobulkK;
% name of the continuation parameter
namepar='rhobulkK';

% min max and step for extending 1par branch
min_bound = 0.02;
max_bound = 0.98;
max_step = 0.001;
% number of attempts when extending 1par branch
bf_attempts = 1000;

% Parameters needed for 2par bifurcation.
contpar_2D_1 = ind_rhobulkK;
contpar_2D_2 = ind_rhobulkD;
max_bound_2D_1 = 0.98;
max_bound_2D_2 = 0.98;
min_bound_2D_1 = 0.02;
min_bound_2D_2 = 0.02;
max_step_2D_1 = 5e-3;
max_step_2D_2 = 5e-3;

% This is used when continuing the 2par diagram.
% The sign shows the direction.
step_direction = -0.001;

% Figure name (name of the window)
figure_name = sprintf("rhobulkK rhobulkD N=%g", N);
% Title name
title_name = sprintf('N=%g', N);
% Number of attempts when continuing.
num_attempts = 700;
% x and y range for the figure
x_range = [0,1];
y_range = [0,1];

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





