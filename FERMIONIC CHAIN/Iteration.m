tic;

% Set up the Hamiltonian Matrix elements for current It.
H(1:L*F_NS,1:L*F_NS)   =  0.0;
H(1:L,L+1:2*L)         =  t*fdag_u(1:L,1:L);
H(2*L+1:3*L,3*L+1:4*L) = -t*fdag_u(1:L,1:L);
H(1:L,2*L+1:3*L)       =  t*fdag_d(1:L,1:L);
H(L+1:2*L,3*L+1:4*L)   =  t*fdag_d(1:L,1:L);

% Symmetrize the Hamiltonian 
H = H+H'; 
% Set up main diagonal entries
for i=1:L
    H(i,i) = E(i)*Lambda^0.5; H(L+i,L+i) = E(i)*Lambda^0.5;
    H(2*L+i,2*L+i) = E(i)*Lambda^0.5; H(3*L+i,3*L+i) = E(i)*Lambda^0.5;
end

% Set up the matrices for fdag_{u,d} for current It.
fdag_u(1:L*F_NS,1:L*F_NS) = 0.0; fdag_d(1:L*F_NS,1:L*F_NS) = 0.0;
fdag_u(L+1:2*L,1:L)   = eye(L); fdag_u(3*L+1:4*L,2*L+1:3*L) =  eye(L);
fdag_d(2*L+1:3*L,1:L) = eye(L); fdag_d(3*L+1:4*L,L+1:2*L)   = -eye(L);

% Assign the Kodes for all the basis states in terms of Kodes for prev. It.
Kode(1:L)       = Kodep(1:L) + delta_enkode(dSz(1),dQ(1));
Kode(L+1:2*L)   = Kodep(1:L) + delta_enkode(dSz(2),dQ(2));
Kode(2*L+1:3*L) = Kodep(1:L) + delta_enkode(dSz(3),dQ(3));
Kode(3*L+1:4*L) = Kodep(1:L) + delta_enkode(dSz(4),dQ(4));

L=L*F_NS;  % Increase No. of basis states by F_NS

t1 = toc; tic;

% -------------------------------------------------------------------------
% Block Diagonalize the Hamiltonian.

% xK : Finds all the unique Kode values.
% nK : No. of unique Kode values.
% iK : Counter that counts the blocks/states computed.
xK = unique(Kode(1:L)); nK = length(xK); iK = 0;
E(1:L) = 0.0; c(1:L,1:L) = 0.0; HK(1:L,1:L) = 0.0;

% Loop over all blocks
for blk = 1:nK
    % Find the position of all the basis states for the Blk.
    % b = (b1,b2,..) stores the Pos. of basis states in the original H.
    % LK : Dim. of the block k.
    b = find(Kode(1:L) == xK(blk)); LK = length(b);
    
    % Form Hamilton Matrix HQ for the current block k
    for i = 1:LK
        for j = 1:LK
            HK(i,j) = H(b(i),b(j));
        end
    end
    
    % Find eigenvalues and eigenvectors of the Blk. k
    [cK(1:LK,1:LK), EK(1:LK)] = eig(HK(1:LK,1:LK),'vector');
    
    % Save the Eig. Vec. and Eig. Energies in the original Matrices.
    for i = 1:LK
        for j = 1:LK
            c(b(i),b(j)) = cK(i,j);
        end
        E(b(i)) = EK(i);
        Kodep(b(i)) = xK(blk);
    end
    iK = iK + LK;
end

t2 = toc; tic;

%--------------------------------------------------------------------------
% Sort the states in ascending order in Energies and shift the Energies.

% Sort Energies + c + Kode
[E(1:L),c(1:L,1:L),Kodep(1:L)] = Esort(E(1:L),c(1:L,1:L),Kodep(1:L),L);

% Keep only the first KEPT states and ignore the other states
K = L;
if(L>KEPT)
    L=KEPT;
end

% Imposing degeneracy on nearly degenate energy states
if (IMPOSE)
    E(1:L) = Impose(E(1:L),L,SMALL);
end

% Shift energies such that the ground state has an energy of 0
E(1:L) = E(1:L) - E(1);

t3 = toc; tic;

%--------------------------------------------------------------------------

% Rotate the Matrices for the fdag_{u,d} Operators - to be used in next It.
fdag_u(1:L,1:L) = c(1:K,1:L)'*fdag_u(1:K,1:K)*c(1:K,1:L);
fdag_d(1:L,1:L) = c(1:K,1:L)'*fdag_d(1:K,1:K)*c(1:K,1:L);

% -------------------------------------------------------------------------

% Extract Qtm. Nos. Sz and Q for printing
[Sz(1:L),Q(1:L)] = dekode(Kodep(1:L));

% Compute Thermodynamics of the system
if (THERMO)
    Thermodynamics();
end

t4 = toc;

% End of Iteration. Print out energy levels and Thermodynamics
ITend();

% The times calculated are as follows:
% t1: Set up H, fdag_{u,d}, Kodep; t2 : Block diagonalizing H;
% t3: Sort and shift E, t4: Thermodynamics + Opmat + Printing.
fprintf(FNOUT, 'Elapsed Time = %8.5f + %8.5f +%8.5f + %8.5f = %8.5f \n\n', t1, t2, t3, t4, t1+t2+t3+t4 );
