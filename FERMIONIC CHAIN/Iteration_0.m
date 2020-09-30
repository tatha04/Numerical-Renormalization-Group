it = 0;      % Iteration 0

L = 4;       % Possible No. of Config. of the f0 site.

% Set up the matrices for fdag_{u,d} for current It.
fdag_u(1:L,1:L) = [0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 1 0];
fdag_d(1:L,1:L) = [0 0 0 0; 0 0 0 0; 1 0 0 0; 0 -1 0 0];

% Find out the (trivial) eigenstates and eigenenergies of the system.
[c(1:L,1:L),E(1:L)] = eig(H(1:L,1:L),'vector');

% Rotate the Matrices for the fdag_{u,d} Operators - to be used in next It.
fdag_u(1:L,1:L) = c(1:L,1:L)'*fdag_u(1:L,1:L)*c(1:L,1:L);
fdag_d(1:L,1:L) = c(1:L,1:L)'*fdag_d(1:L,1:L)*c(1:L,1:L);

% Define Q and 2Sz qtm. nos for f0 site. 
Q(1:L)  = [-1 0 0 1];
Sz(1:L) = [0 1 -1 0];

% Assign Kodes to the eigenstates.
Kodep(1:L) = enkode(Sz(1:L),Q(1:L));

% Compute thermodynamics of the system
if(THERMO)
    Thermodynamics();
end

% End of Iteration. Print out energy levels and Thermodynamics
ITend();