%% Setup the first fixed point
stst.kind = 'stst';
stst.parameter = [rhobulkK, rhobulkD, pE, wE, pI, wI, vK, vD, N];

%% Use the solution given by dde23 to approximate the steady state.
JK = vK*J(rhobulkK);
JD = vD*J(rhobulkD);
tauK = delay(rhobulkK, N, vK);
tauD = delay(rhobulkD, N, vD);
delays = [tauK, tauD];
sol = dde23(@rhs_dde23, delays, @history, [0, 200]);

E1_st=sol.y(1,end);
E2_st=sol.y(2,end);
I1_st=sol.y(3,end);
I2_st=sol.y(4,end);
stst.x = [E1_st; E2_st; I1_st; I2_st];


%% Correct the approximations above.
method = df_mthod(funcs, 'stst');
% If there is an error message caused by running the following line,
% you should check set_funcs to see whether the parameter positions are
% correct or not.
[stst, success] = p_correc(funcs, stst, [], [], method.point);
E1_st = stst.x(1);
E2_st = stst.x(2);
I1_st = stst.x(3);
I2_st = stst.x(4);
% E1_st, E2_st, I1_st, I2_st will be used in the next script: main_initialize_branch