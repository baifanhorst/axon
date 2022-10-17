% The hill function
% Used to describe excitation or inhibition

function y=hill(x,n,K)

y = x^n/(x^n+K^n);

end