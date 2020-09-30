% Definition of Operator(s) whose expectation values are to be calculated.
% Op = <S1.S2>
for i = 0:(L-1)
    z = dec2dual(i,NS);
    
    p = 1;
    q = 2;
    
    op_n(i+1,i+1) = op_n(i+1,i+1) + 0.25*(2*z(p)-1)*(2*z(q)-1);
    
    if (z(p) + z(q) == 1)
        i_prime = i + (1-2*z(p))*2^(p-1) + (1-2*z(q))*2^(q-1);
        op_n(i_prime+1,i+1) = op_n(i_prime+1,i+1) + 0.5;
    end
end