function tau = sys_tau(delay_no, xx, par)
% This function calculates the delays. If delay_no=1, tau1 will be
% calculated. If delay_no=2, tau2 will be calculated.
% xx is the state of the system. including E1(t), E2(t), I1(t), I2(t) as the 1st column,
% E1(t-tau1), E2(t-tau1), I1(t-tau1), I2(t-tau1) as the 2nd column, and 
% E1(t-tau2), E2(t-tau2), I1(t-tau2), I2(t-tau2) as the 3rd column.
% [par] is the parameter list. In this case,
% par = [rhobulkK, rhobulkD, pE, wE, pI, wI, vK, vD, N].
rhobulkK = par(1);
rhobulkD = par(2);
%pE = par(3);
%wE = par(4);
%pI = par(5);
%wI = par(6);
vK = par(7);
vD = par(8);
N = par(9);

if delay_no==1
    tau = N/vK/(1-rhobulkK);
elseif delay_no==2
    tau = N/vD/(1-rhobulkD);
end

end