function v = rhs_dde23(t,y,Z)

global pE dE1 dE2 wE nE KE betaE
global pI dI1 dI2 wI nI KI
global JK JD tauK tauD
global N
global vK vD

E1 = y(1);
E2 = y(2);
I1 = y(3);
I2 = y(4);

E1_delay = Z(1,1);
I2_delay = Z(4,2);

pE1 = pE*(1-betaE*hill(I1, nE, KE));
pI2 = pI*hill(E2, nI, KI);

v = zeros(4,1);

v(1) = pE1 - dE1*E1 - JK*wE*E1;
v(2) = -dE2*E2 + JK*wE*E1_delay;
v(3) = -dI1*I1 + JD*wI*I2_delay;
v(4) = pI2 - dI2*I2 - JD*wI*I2;


end