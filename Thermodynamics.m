beta = Lambda^(+0.5); temp(it+1) = (Lambda^(-(it-1)/2))/beta;

z = sum(exp(-beta*E(1:L)));

s(it+1)    = beta*(sum(E(1:L).*exp(-beta*E(1:L)))/z) + log(z);
tchi(it+1) = 0.25*sum(Sz(1:L).*Sz(1:L).*exp(-beta*E(1:L)))/z - 0.5*sum(Sz(1:L).*exp(-beta*E(1:L)))/z;
Cv(it+1)   = sum(E(1:L).*E(1:L).*exp(-beta*E(1:L)))/z - sum(E(1:L).*exp(-beta*E(1:L)))/z;

if(OPMAT)
    for i = 1:NOP
        op_exp(it+1,i) = sum(diag(op_n(1:L,1:L,i)).*exp(-beta*E(1:L)))/z;
    end
end