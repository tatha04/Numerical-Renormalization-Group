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
Iteration_0();

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
        savg(it+1)   = 0.25*(s(it)+2*s(it+1)+s(it+2));
        tchiavg(it+1)= 0.25*(tchi(it)+2*tchi(it+1)+tchi(it+2));
        Cvavg(it+1)   = 0.25*(Cv(it)+2*Cv(it+1)+Cv(it+2));
        fprintf(FNTHERMAVG, '%14.6e %14.8f %14.8f %14.8f \n', temp(it+1), savg(it+1), tchiavg(it+1), Cvavg(it+1));
    end
end

% Print the total time for the run.
tf = cputime;
fprintf(FNOUT, 'Total time = %8.4f \n', tf-t0);

% Close all output files
Close_Files();

% Plot Energies and Thermodynamics
Plot_Figures();

end
