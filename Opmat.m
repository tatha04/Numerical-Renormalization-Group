op(1:L*F_NS,1:L*F_NS,1:NOP) = 0.0;
op(1:L,1:L,1:NOP) = op_n(1:L,1:L,1:NOP);
op(L+1:2*L,L+1:2*L,1:NOP) = op_n(1:L,1:L,1:NOP);
op(2*L+1:3*L,2*L+1:3*L,1:NOP) = op_n(1:L,1:L,1:NOP);
op(3*L+1:4*L,3*L+1:4*L,1:NOP) = op_n(1:L,1:L,1:NOP);

