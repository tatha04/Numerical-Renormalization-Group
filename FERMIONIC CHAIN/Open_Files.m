% Open output files =======================================================

% out.dat : Print out Energies
% even.dat/odd.dat : Print out Energies for even/odd It.
% Therm.dat : Print out Thermodynamic qts. <S>, <TChi>, <Cv> vs Temp.
% Thermavg.dat : Print out 3-pt. avgs <Savg.>, <TChiavg>, <Cvavg> vs Temp.
% OpTherm.dat: Print out the expectation values of the Imp. Operators.

%--------------------------------------------------------------------------

FNOUT   = fopen('out.dat' ,'w+');
if (EO_PRINT)
    FNEVEN  = fopen('even.dat','w+');  FNODD = fopen('odd.dat' ,'w+');
end
if(THERMO)
    FNTHERM = fopen('Therm.dat','w+'); FNTHERMAVG = fopen('Thermavg.dat','w+');
end

% -------------------------------------------------------------------------

% Print out header information ============================================

header_line = '1-Channel Fermionic Chain \n';

fprintf(FNOUT,header_line); clear header_line;

fprintf(FNOUT, 'Lambda = %5.3f, KEPT = %4i, ITMAX = %3i \n \n', Lambda,KEPT,ITMAX);

% Therm.dat: <S>, <TChi>, <Cv>
if (THERMO)
    fprintf(FNTHERM,   '# %3s %10s %12s %16s %12s \n', 'It','Temp.','<S>','<TChi>','Cv');
    fprintf(FNTHERM,   '#---------------------------------------------------------------\n');
    fprintf(FNTHERMAVG,'# %10s %12s %16s %12s \n', 'Temp.','S_avg','TChi_avg','Cv_avg');
    fprintf(FNTHERMAVG,'#- --------------------------------------------------------\n');
end

% even.dat + odd.dat
if (EO_PRINT)
    fprintf(FNODD, '#it'); fprintf(FNEVEN, '#it');
    for i = 1:ENMAX_PR
        fprintf(FNODD, '    E%2i ',i); fprintf(FNEVEN, '   E%2i  ',i);
    end
    fprintf(FNODD, '\n'); fprintf(FNEVEN, '\n');
end

%--------------------------------------------------------------------------
%==========================================================================