function NRG()

% Load Input parameters.
Input_Variables();

% Declare Run Variables.
Parms();

% Open output files
Open_Files();

% Set timer
t0 = cputime;

% Iteration 0
if (strcmp(MODEL,'KO'))
    Iteration_0_Kondo();
elseif (strcmp(MODEL,'AN'))
    Iteration_0_Anderson();
elseif (strcmp(MODEL,'2IK'))
    Iteration_0_2IK();
end

% =========================================================================

% NRG Iterations
for it=1:ITMAX
    
    % Set up hopping coeffecients
    n = it-1;
    t = (1+Lambda^-1)*(1-Lambda^(-n-1))/(2*sqrt(1-Lambda^(-2*n-1))*sqrt(1-Lambda^(-2*n-3)));
          
    % Set up the Hamiltonian and Diagonalize.
    Iteration();
       
end

if (THERMO)
    % Print out average thermodynamics
    for it=1:ITMAX-1
        savg(it+1)   = 0.25*(s(it)+2*s(it+1)+s(it+2)); simp(it+1) = savg(it+1) - s0(it+1);
        tchiavg(it+1)= 0.25*(tchi(it)+2*tchi(it+1)+tchi(it+2)); tchiimp(it+1) = tchiavg(it+1) - tchi0(it+1);
        Cvavg(it+1)  = 0.25*(Cv(it)+2*Cv(it+1)+Cv(it+2)); Cvimp(it+1) = Cvavg(it+1) - Cv0(it+1);
        fprintf(FNTHERMAVG, '%14.6e %14.8f %14.8f %14.8f \n', temp(it+1), simp(it+1), tchiimp(it+1), Cvimp(it+1));
    end
    
    % Find T_k/T_c
    if(FIND_Tk)
        Tk = 0.0;
        found = false; i = 2; simp_cr = 0.5*(Simp_UL+Simp_LL);
        while (found == false && i<ITMAX-1 )
            if(simp(i) > simp_cr && simp(i+1) < simp_cr)
                found = true;
                Tk = exp(log(temp(i+1)) + (log(temp(i))-log(temp(i+1)))*(simp_cr-simp(i+1))/(simp(i)-simp(i+1)));
            end
            i = i+1;
        end
        fprintf(FNOUT,'Kondo/Crossover Temp. T_K = %12.4e \n', Tk);
    end
end

% Print the total time for the run.
tf = cputime;
fprintf(FNOUT, 'Total time = %8.4f \n', tf-t0);

% Close all output files
Close_Files

% Plot Energies and Thermodynamics
Plot_Figures

end
