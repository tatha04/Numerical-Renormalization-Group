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
if (OPMAT)
    FNOPMAT = fopen('OpTherm.dat', 'w+');
end

% -------------------------------------------------------------------------

% OPTIONAL ----------------------------------------------------------------
% Choose filename according to Jk (or other parameter) values.

% file_out   = num2str(Jk,'out_Jk=%5.3f.dat');
% FNOUT   = fopen(file_out ,'w+');
% 
% if (EO_PRINT)
%     file_even  = num2str(Jk,'even_Jk=%5.3f.dat');
%     file_odd   = num2str(Jk,'odd_Jk=%5.3f.dat');
%     FNEVEN  = fopen(file_even,'w+');  FNODD = fopen(file_odd ,'w+');
% end
% 
% if(THERMO)
%     file_therm     = num2str(Jk,'Therm_Jk=%5.3f.dat');
%     file_thermavg  = num2str(Jk,'Thermavg_Jk=%5.3f.dat');
%     FNTHERM = fopen(file_Therm,'w+'); FNTHERMAVG = fopen(file_thermavg,'w+');
% end
% 
% if (OPMAT)
%     file_optherm = num2str(Jk,'OpTherm_Jk=%5.3f.dat');
%     FNOPMAT = fopen(file_optherm, 'w+');
% end
%--------------------------------------------------------------------------

% Print out header information ============================================

% out.dat:
if (strcmp(MODEL,'KO'))
    header_line = '1-Imp. 1-Channel Kondo Model \n';
elseif (strcmp(MODEL,'AN'))
    header_line = '1-Imp. 1-Channel Anderson Model \n';
elseif (strcmp(MODEL,'2IK'))
    header_line = '2-Imp. 1-Channel Kondo Model \n';
end

fprintf(FNOUT,header_line); clear header_line;

if (strcmp(MODEL,'KO'))
    fprintf(FNOUT,'Jk = %8.4f \n', Jk);
elseif (strcmp(MODEL,'AN'))
    fprintf(FNOUT,'Ed = %8.4f, U = %8.4f, Gamma = %8.4f \n ', ed, U, Gamma);
elseif (strcmp(MODEL,'2IK'))
    fprintf(FNOUT,'J1 = %8.4f, J2 = %8.4f, J12 = %8.4f \n', Jk1, Jk2, J12);
end

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

if (THERMO)
    % Read Thermodynamics of fermionic chain from band file <therm.data>
    
    delimiterIn = ' '; headerlinesIn = 2;
    A = importdata('therm.data', delimiterIn,headerlinesIn);
    [band_itmax, ncol]=size(A.data);
    for i = 1:band_itmax
        s0(i+1)= A.data(i,2); tchi0(i+1)= A.data(i,3); Cv0(i+1)= A.data(i,4);
        %temp0(it) = A.data(i-1,2);
    end
    
    % Check if ITMAX of band file = ITMAX of current run.
    % If not, make sure (manually) that band file has converged and fill the
    % rest of the values with the last entry.
    if (ITMAX > band_itmax)
        for it = band_itmax+1:ITMAX
            s0(it)  = s0(band_itmax); tchi0(it)  = tchi0(band_itmax); Cv0(it)  = Cv0(band_itmax);
        end
    end
    clear A;
end

%==========================================================================