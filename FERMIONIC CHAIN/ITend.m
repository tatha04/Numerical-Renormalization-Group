% Print out the energies in output file

fprintf(FNOUT, '%s %i \n','Iteration =', it);
fprintf(FNOUT,'---------------\n');
fprintf(FNOUT, ' %3s %2s %5s %3s %10s \n', 'It.', '#', '2Sz', 'Q', 'E');
fprintf(FNOUT,'--------------------------------------\n');
for i=1:min(KEPT,L)
    fprintf(FNOUT, '%3i %3i %4i %4i %20.16f \n', it, i, Sz(i), Q(i), E(i));
end
fprintf(FNOUT, '\n');

%  Save even/odd energies for plotting

if(PLOT_E || EO_PRINT)
    if( mod(it,2) == 0 )
        EvenIt((it/2)+1) = it;
        for i = 1:min(L,ENMAX_PL)
            EvenE((it/2)+1,i) = E(i);
        end
    else
        OddIt(int16(it/2)) = it;
        for i = 1:min(L,ENMAX_PL)
            OddE(int16(it/2),i) = E(i);
        end
    end
end

% Print Even/Odd Energies

if(EO_PRINT)
    if( mod(it,2) == 0 )
        NP = FNEVEN;
    else NP = FNODD;
    end
    fprintf(NP,'%3i', it);
    for i = 1:min(L,ENMAX_PR)
        fprintf(NP, '%8.4f', E(i));
    end
    fprintf(NP,'\n');
end

% Print out Thermodynamic qts.

if(THERMO)
    fprintf(FNTHERM, '%4i %14.6e %14.8f %14.8f %14.8f \n', it, temp(it+1), s(it+1), tchi(it+1), Cv(it+1));
end
