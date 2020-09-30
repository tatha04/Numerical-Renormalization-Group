% Input Parameters ========================================================
%
% MODEL : Type of Model
% KO (KONDO), AN (ANDERSON), 2IK (2-IMP. KONDO)
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
% FIND_Tk : Calculate T_k/T_c. Use <Simp> to determine T_k/T_c.
% Simp_UL / Simp_LL : Simp_(UpperLower)_Limit.
% Simp(T_k) = (Simp_UL + Simp_UL)/2
%
% Kondo Model:
%   Jk : Kondo coupling between site-f0 and Imp. Simp.
% Anderson Model:
%   Gamma : Hybridization between site-f0 and Imp. site d
%   ed = \epsilon_d : Imp. occupation energy.
%   U : Coulomb repulsion between ^ and v spins at Imp.
% 2-Imp Kondo:
%   Jk1 : Kondo coupling between site-f0 and Imp. Spin S1.
%   Jk2 : Kondo coupling between site-f0 and Imp. Spin S2.
%   J12 : Exchange coupling between spins 1 and 2.
%
% himp = Imp. magnetic field (TO BE IMPLEMENTED)
%
% PLOT_E(T/F) : Plot Odd and Even Energies at the end of the run.
% PLOT_T(T/F) : Plot Thermodynaic Qts. at the end of the run.
% EO_PRINT(T/F) : Print Even/Odd Energies at the end of each It.
%
% ENMAX_PL : Number of states to plot for even/odd iterations.
% ENMAX_PR : Number of states to print for even/odd iterations.
% ENMAX_PR <= ENMAX_PL
% =========================================================================

MODEL = 'KO';

% Kondo:
%Jk = 0.0;

% Anderson Model
%Gamma = 0.1; ed = -0.2; U = 0.4;

% 2 Imp. Kondo
%Jk1 = 0.25; Jk2 = 0.25; J12 = 0.01;

KEPT = 800; ITMAX = 60; Lambda = 3.0; 

IMPOSE = true; SMALL = 1e-6; 

THERMO = true;

OPMAT = false; NOP = 1;

PLOT_T = false; PLOT_E = false; EO_PRINT = false;

ENMAX_PL = 32; ENMAX_PR = 16;

FIND_Tk = true; Simp_UL = log(4); Simp_LL = log(2);

% =========================================================================

