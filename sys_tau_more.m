function tau = sys_tau_more(delay_no, xx, par)

% par = [alphaK, alphaD, pE, wE, pI, wI, vK, vD, N, dK, wK, dD, wD]

alphaK = par(1);
alphaD = par(2);
%pE = par(3);
%wE = par(4);
%pI = par(5);
%wI = par(6);
vK = par(7);
vD = par(8);
N = par(9);
dK = par(10);
wK = par(11);
dD = par(12);
wD = par(13);


gammaK = alpha_to_gamma3(alphaK, dK, wK, vK);
rhobulkK = density_bulk(alphaK, gammaK);

gammaD = alpha_to_gamma3(alphaD, dD, wD, vD);
rhobulkD = density_bulk(alphaD, gammaD);



if delay_no==1
    tau = N/vK/(1-rhobulkK);
elseif delay_no==2
    tau = N/vD/(1-rhobulkD);
end

end