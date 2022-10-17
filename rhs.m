function v = rhs(xx, par)
% This function gives the right-hand sides of the system.
% [xx] is the state of the system, 
% including E1(t), E2(t), I1(t), I2(t) as the 1st column,
% E1(t-tau1), E2(t-tau1), I1(t-tau1), I2(t-tau1) as the 2nd column, and 
% E1(t-tau2), E2(t-tau2), I1(t-tau2), I2(t-tau2) as the 3rd column.
% [par] is the parameter list. In this case,
% par = [rhobulkK, rhobulkD, pE, wE, pI, wI, vK, vD, N].

% These are fixed parameters.
global dE1 dE2 nE KE betaE
global dI1 dI2 nI KI
% global N
% global wa_K wd_K wa_D wd_D K_K K_D R_K R_D
% global M_K M_D eps
% global dK wK Vtip

% v represents the right-hand sides.
v = zeros(4,1);

% Use convenient names for E1, E2, I1, I2, E1(t-tau1) and I2(t-tau2).
% Thus we can write down the right-hand sides easier (see the end).
E1 = xx(1,1);
E2 = xx(2,1);
I1 = xx(3,1);
I2 = xx(4,1);

E1_delay = xx(1,2); %E1_delay = E1(t-tau1)
I2_delay = xx(4,3); %I2_delay = I2(t-tau2)

% Use convenient parameter names rather than using [par] directly.
rhobulkK = par(1);
rhobulkD = par(2);
pE = par(3);
wE = par(4);
pI = par(5);
wI = par(6);
vK = par(7);
vD = par(8);
%N = par(9);

% Use density_bulk to calculate J.
JK = vK*J(rhobulkK);
JD = vD*J(rhobulkD);

% Calculate the production rates of E1 and I2.
pE1 = pE*(1-betaE*hill(I1, nE, KE));
pI2 = pI*hill(E2, nI, KI);



% Write down the right-hand sides.
v(1) = pE1 - dE1*E1 - JK*wE*E1;
v(2) = -dE2*E2 + JK*wE*E1_delay;
v(3) = -dI1*I1 + JD*wI*I2_delay;
v(4) = pI2 - dI2*I2 - JD*wI*I2;

end