% Input Parameters ========================================================
%
% MODEL : Type of Model
% FC (Fermionic chain)
%
% ITMAX : Max. No. of It. it = 0,1,2,...,ITMAX.
% KEPT : Max. No. of states kept after each It.
% Lambda : Discretization Parameter, t_n ~ \Lambda^{-n/2}
%
% IMPOSE(T/F) : Impose degeneracy of nearly degenerate states.
% SMALL  : Tolerance for IMPOSE.
%
% PLOT(T/F) : Plot Energies and Thermodynamics at the end of the run. 
% 
% THERMO(T/F) : Calculate Thermodynamics
%
% OPMAT(T/F) : Calculate Matrix. Elements of operator(s) \hat Op.
% NOP : No. of Operators to calculate Opmat.
%
% PLOT_E(T/F) : Plot Odd and Even Energies at the end of the run.
% PLOT_T(T/F) : Plot Thermodynaic Qts. at the end of the run.
% EO_PRINT(T/F) : Print Even/Odd Energies at the end of each It.
%
% ENMAX_PL : Number of states to plot for even/odd iterations.
% ENMAX_PR : Number of states to print for even/odd iterations.
% ENMAX_PR <= ENMAX_PL
% =========================================================================

MODEL = 'FC';

KEPT = 600; ITMAX = 30; Lambda = 3.0; 

IMPOSE = true; SMALL = 1e-10;

THERMO = true;

PLOT_T = false; PLOT_E = false; EO_PRINT = false;

ENMAX_PL = 16; ENMAX_PR = 16;

% =========================================================================

