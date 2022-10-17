function v = rhs_more(xx, par)
% par = [alphaK, alphaD, pE, wE, pI, wI, vK, vD, N, dK, wK, dD, wD]
global dE1 dE2 nE KE betaE
global dI1 dI2 nI KI
% global N
% global wa_K wd_K wa_D wd_D K_K K_D R_K R_D
% global M_K M_D eps
% global dK wK Vtip


v = zeros(4,1);

E1 = xx(1,1);
E2 = xx(2,1);
I1 = xx(3,1);
I2 = xx(4,1);

E1_delay = xx(1,2);
I2_delay = xx(4,3);





% Use tau to calculate density_bulk. Then use density_bulk to calculate J.
alphaK = par(1);
alphaD = par(2);
pE = par(3);
wE = par(4);
pI = par(5);
wI = par(6);
vK = par(7);
vD = par(8);
%N = par(9);
dK = par(10);
wK = par(11);
dD = par(12);
wD = par(13);

gammaK = alpha_to_gamma3(alphaK, dK, wK, vK);
rhobulkK = density_bulk(alphaK, gammaK);

gammaD = alpha_to_gamma3(alphaD, dD, wD, vD);
rhobulkD = density_bulk(alphaD, gammaD);

JK = vK*J(rhobulkK);
JD = vD*J(rhobulkD);


pE1 = pE*(1-betaE*hill(I1, nE, KE));
pI2 = pI*hill(E2, nI, KI);





v(1) = pE1 - dE1*E1 - JK*wE*E1;
v(2) = -dE2*E2 + JK*wE*E1_delay;
v(3) = -dI1*I1 + JD*wI*I2_delay;
v(4) = pI2 - dI2*I2 - JD*wI*I2;

end