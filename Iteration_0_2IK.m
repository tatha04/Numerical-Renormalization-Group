it = 0;      % Iteration 0

tic;

% This code is a modified version of a problem that deals with a generic
% Imp. cluster consisting of spin d.o.f's.
% NS : No. of spins in the Imp. cluster.
% LS : No. of links in the spin cluster.
% LF = No. of links between f0 site and Imp. cluster.
NS = 2; LS = 3; LF = 2;

site1 = zeros(1,LS);
site2 = zeros(1,LS);
alpha = zeros(1,LS);
J     = zeros(1,LS);

% Def. of Kitaev Cluster
site1(1:LS) = [1 1 1];
site2(1:LS) = [2 2 2];
alpha(1:LS) = [1 2 3];
J(1:LS) = J12;

% Def. of f0- links to the cluster
sitef = zeros(1,LF); JK = zeros(1,LF);
% It is assumed that sitef = [1 2 3 4 ..] etc.
% It needs to be changed according to the resp. model.
for i = 1:LF
    sitef(i) = i;
end
JK(1)   = Jk1; JK(2)   = Jk2;

L = 2^NS;

% Set up the Hamiltonian for the Imp. cluster.
for i = 0:(L-1)
    z = dec2dual(i,NS);
    for j = 1:LS
        p = site1(j);
        q = site2(j);
        i_prime = i + (1-2*z(p))*2^(p-1) + (1-2*z(q))*2^(q-1);
        if (alpha(j) == 1)
            H(i_prime+1,i+1) = H(i_prime+1,i+1) + 0.25*J(j);
        end
        if (alpha(j) == 2)
            H(i_prime+1,i+1) = H(i_prime+1,i+1) - 0.25*J(j)*(2*z(p)-1)*(2*z(q)-1);
        end
        if (alpha(j) == 3)
            H(i+1,i+1) = H(i+1,i+1) + 0.25*J(j)*(2*z(p)-1)*(2*z(q)-1);
        end
    end
end

% Copy the Imp. Hamiltonian into Diagonal Blocks.
H(L+1:2*L,L+1:2*L)     = H(1:L,1:L);
H(2*L+1:3*L,2*L+1:3*L) = H(1:L,1:L);
H(3*L+1:4*L,3*L+1:4*L) = H(1:L,1:L);

% Set up the Hamiltonian for Imp.-f0 links
for i = 0:(L*F_NS-1)
    k = mod(i,2^NS); % k = 0,1: spin config. of Imp. site(s).
    z = dec2dual(k,NS); % z : spin rep. of Imp. site(s)
    q = (i-k)/2^NS; % q = 0,1,2,3 : config. of f0 site.
    
    for j = 1:LF
        p = sitef(j); % p = Imp. cluster site index.
        
        % z component interactions
        if(q == 1 || q == 2) % q = 1 : ^, q = 2 : v
            H(i+1,i+1) = H(i+1,i+1) + 0.25*JK(j)*(2*z(p)-1)*(3-2*q);
        end
        % x,y components interactions.
        if((2*z(p)-1) + (3-2*q) == 0)
            i_prime = (3-q)*2^NS + k + (1-2*z(p))*2^(p-1);
            H(i_prime+1,i+1) = H(i_prime+1,i+1) + 0.5*JK(j);
        end
    end
end

% Set up Qtm. Nos. 2Sz and Q for Imp. cluster.
% For Kitaev Imps. Sz is not conserved. So Sz is assumed to be 0 for all
% basis states.
for i = 1:L
    z = dec2dual(i-1,NS);
    Sz(i) = sum(2*z-1); Q(i) = 0;
end
Kodep(1:L) = enkode(Sz(1:L),Q(1:L));

% Set up the matrices for fdag_{u,d} for It. 0.
fdag_u(1:L*F_NS,1:L*F_NS) = 0.0; fdag_d(1:L*F_NS,1:L*F_NS) = 0.0;
fdag_u(L+1:2*L,1:L)   = eye(L); fdag_u(3*L+1:4*L,2*L+1:3*L) =  eye(L);
fdag_d(2*L+1:3*L,1:L) = eye(L); fdag_d(3*L+1:4*L,L+1:2*L)   = -eye(L);

% Assign the Kodes for all the basis states in terms of Kodes for prev. It.
Kode(1:L)       = Kodep(1:L) + delta_enkode(dSz(1),dQ(1));
Kode(L+1:2*L)   = Kodep(1:L) + delta_enkode(dSz(2),dQ(2));
Kode(2*L+1:3*L) = Kodep(1:L) + delta_enkode(dSz(3),dQ(3));
Kode(3*L+1:4*L) = Kodep(1:L) + delta_enkode(dSz(4),dQ(4));

if (OPMAT)
    Opmat0_2IK(); Opmat();
end

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

E(1:L) = E(1:L)*Lambda^(-0.5);

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

% Compute Operator Matrix elements
if (OPMAT)
    % Rotate Matrix using Unitary transf. and save it for next It.
    op_n(1:L,1:L) = c(1:K,1:L)'*op(1:K,1:K)*c(1:K,1:L);
end

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
