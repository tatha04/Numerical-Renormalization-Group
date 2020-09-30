% Declaration of Run variables ============================================
%
% RUN VARIABLES------------------------------------------------------------
% F_NS : D.o.f of each site of the fermionic chain.
% For a One-channel chain, the states are {|0>, |1>, |-1>, |2>}.
%
% MAXDIM = F_NS*KEPT : Max. dim. of Hamiltonian matrix.
% H(MAXDIM,MAXDIM) : Hamiltonian Matrix.
% c(MAXDIM,MAXDIM) : Coeff. of eigenvectors after diagonalization.
% E(MAXDIM) : Array that stores Energies before and after sorting.
% MAXBLKDIM = KEPT: Max, dim. of a Hamiltonian block.
% HK(MAXBLKDIM,MAXBLKDIM) : Hamiltonain Matrix. of a Blk.
% cK(MAXBLKDIM,MAXBLKDIM) : Coeff. of eigenvectors after Blk. diag.
% EK(MAXBLKDIM) : Array that stores Energies after Blk. diag.
% fdag_u(MAXDIM,MAXDIM) : Matrix rep. of fdag_u at each It.
% fdag_d(MAXDIM,MAXDIM) : Matrix rep. of fdag_u at each It.
%
% QTM NOS.-----------------------------------------------------------------
% Kode = enkode(Sz,Q) stores the qtm. numbers 2Sz and Q together.
% They can be retrieved using [Sz,Q] = dekode(Kode)
% Kode(MAXDIM) : Stores the Kode of of each state in the current iteration.
% Kodep(MAXDIM): Stores the Kode of of each state in the prev. iteration.
% Sz(MAXDIM) : Stores Qtm. No. 2Sz of the eigenstates.
% Q(MAXDIM)  : Stores charge Qtm. No. Q = n^-nv of the eigenstates.
% Sz_good(T/F): Sz is a good Qtm. No.
% Q_good(T/F) : Q is a good Qtm. No.
% dSz/dQ : Change in Sz Qtm. No. for the new basis states.
% (by taking a cross product with the new f0-site).
%
% THERMODYNAMICS-----------------------------------------------------------
% temp(ITMAX+1): Array that stores temp. of each It.
% s(ITMAX+1)   : Array that stores Entropy <S> at each It.
% savg(ITMAX+1): Array that stores 3-pt. avg. of <S>.
% tChi(ITMAX+1): Array that stores <T\Chi> at each It.
% tchiavg(ITMAX+1) : Array that stores 3-pt. avg. of <TChi>.
% Cv(ITMAX+1 )  : Array that stores Specific heat <Cv> at each It.
% Cvavg(ITMAX+1): Array that stores 3-pt. avg. of <Cv>.
% s0, tchi0, Cv0 : Array that stores resp. thermodynamics from band file
% simp, tchiimp, Cvimp : Imp. contribution to resp. thermodynamic qts.
%
% =========================================================================

F_NS  = 4; MAXDIM = F_NS*KEPT; MAXBLKDIM = KEPT;     
H = zeros(MAXDIM); c = zeros(MAXDIM); E = zeros(MAXDIM,1);
HK = zeros(MAXBLKDIM); cK = zeros(MAXBLKDIM); EK = zeros(MAXBLKDIM,1);
fdag_u = zeros(MAXDIM); fdag_d = zeros(MAXDIM);

Kode = zeros(MAXDIM,1); Kodep = zeros(MAXDIM,1);
Sz = zeros(MAXDIM,1)  ; Q = zeros(MAXDIM,1); 

if (THERMO)
    temp = zeros(ITMAX+1,1);
    s    = zeros(ITMAX+1,1); savg    = zeros(ITMAX+1,1);
    tchi = zeros(ITMAX+1,1); tchiavg = zeros(ITMAX+1,1);
    Cv   = zeros(ITMAX+1,1); Cvavg   = zeros(ITMAX+1,1);
end

if (PLOT_E || EO_PRINT)
    EvenE = zeros(int16(ITMAX/2)+1,ENMAX_PL); EvenIt= zeros(int16(ITMAX/2)+1,1);
    OddE  = zeros(int16(ITMAX/2),ENMAX_PL); OddIt = zeros(int16(ITMAX/2),1);  
end

Sz_good = true; Q_good = true;

dSz = [0 0 0 0]; dQ = [0 0 0 0];
if ( Sz_good == true )
    dSz = [0 1 -1 0];
end
if (Q_good == true )
    dQ = [-1 0 0 1];

end

% =========================================================================
