USAGE:
------

1. The program can be run using Matlab. The main function is contained in
	the file NRG.m. Thus to run this, change Matlab's working directory 
	to the directory where all the files are stored.
	In Matlab's command window, type NRG().
	If running from a terminal, Matlab can be run without graphics using:
	matlab -nodisplay -nojvm
	In that case, please set the variables PLOT_T and PLOT_E = false.

2. Edit the file Input_Variables.m and select the appropriate model.
	eg, FC: Fermionic chain, KO: Kondo Model, etc.
	Set the appropriate model parameter(s).
	eg, for the Kondo model, set the value of Jk;
	for the 2Imp. Kondo model, set the values of Jk1, Jk2, J12;
        for the Anderson model, set the values of ed, U, Gamma; etc,..
	and comment out the rest of the input lines.

3. For Imp. Thermodynamics, one needs to subtract off the contribution towards
	the thermodynamic quantity due to the band alone.
        Run NRG for the FERMIONIC CHAIN with THERMO = true.
        The necessary files are contained in the directory './FERMIONIC CHAIN'
	Ensure that ITMAX is high enough such that the <xavg>'s have converged.
	Make a copy or a symbolic link to file Thermavg.dat and save it as therm.data.

4. By default the code assumes that Sz and Q are good Qtm. Nos. In case this
	is not true, make necessary changes in file Parms.m so that the values
	of the relevant Qtm. No is set to zero for all basis states.
	

Notes and Advanced diagnostics:
------------------------------

3. Make sure that the variables KEPT (Total number of Kept states),
	MAXDIM =4*KEPT (Max. Dimension of the Hamiltonian) for a standard
	fermionic chain, and MAXBLKDIM (Max. Dimension of a Block for a
	given set of Quantum Nos.) are consistent.

4. Logical Variables Sz_good and Q_good determines whether Sz and Q are good
	Qtm. Nos. If yes, the the Hamiltonian is broken up into sub-blocks
	labelled by the Qtm. Nos.

5. Function [kode = enkode(Sz,Q)] and [Sz,Q = dekode(kode)] encodes and decodes
	respectively the Qtm. Nos. into a single integer code and vice versa.
	Function delta_enkode(dSz,dQ) calculates the change in the kode of the
	basis state when taking a cross product with a new fermionic site
	with Qtm. Nos dSz and dQ resp.
	The kodes of the eigenstates are then saved into the array Kode.

6. Let the eigen states for It. N be |i> where i = 1,..,L.
	Then the basis states for It. N+1 are as follows (for 1-channel):
	|j> = |i> x |0> , j = 1,...,L;     |j> = |i> x |^> , j = L+1,...,2L;
	|j> = |i> x |v> , j = 2L+1,...,3L; |j> = |i> x |2> , j = 3L+1,...,4L;

7. Iteration 0: The file Iteration 0 is different for different Imp. models.
	The basis states for It. 0 are as follows:
	{[0,..,0], ..,[1,,,,1]} x {|0>, |^>, |v>, |2>},
	where 0:v and 1:^ in the Imp. (spin) part of the Imp. Hamiltonian
	representation [first {}], and |0>, |^>, |v>, |2> denote null, up,
        down and doubly occupied states, respectively, of the f_0 site.

8. For Imp. that are fermionic in nature, the basis states are as follows:
	{|0>, |^>, |v>, |2>} x {|0>, |^>, |v>, |2>}.
	Situations with more than 1 Fermionic Imp. needs to be implemented.


========================================================================

OPTIONAL:

Imp. Thermodynamics: If the parameter THERMO == TRUE, then the Impurity
	contribution to the following quantities are calculated.
	Entropy (S), Magnetic Susceptibility (T*\Chi), and Specific Heat (Cv).
	These are calculated at the end of each Iterations and the values
	are stored in the arrays s, tchi and Cv resp.
	At the end of the NRG run, the 3-point averages are calculated
	and saved in the arrays savg, tchiavg and Cvavg resp.
	Finally the Imp. contribution (simp, tchiimp, and Cvimp) is 
	calculated by subtracting off the contribution due to bath alone.

Find_Tk: Calculate the Kondo Temp. or optionally any crossover Temp.
	defined as the Temp. at which the entropy crosses half-way
	between Simp_UL and Simp_LL.

Opmat: Calculate Operator matrix elements for 1 or more Operators.
	At the end of each iteration, the thermodynamic expectation
	value of the Operator(s) is printed.

=======================================================================

Definition of Output and Input Files:
-------------------------------------

out.dat: Print out the Renormalized energy levels at the end of each NRG
	iteration. Also prints out (i) the time taken for each Iteration,
	and (ii) total time taken for the run.
	The Qtm. Nos 2Sz and Q are also printed beside each energy.
	In case Sz or Q are not good quantum nos, the corresponding values
	are set to 0.

even/odd.dat: Prints out first few energy levels for even/odd Iterations.
	They can then be plotted using standard plotting tools.

Therm.dat: Print out the total values of the Thermodynamic quantities against
	Iteration # as well as Temp.

Theravg.dat: Prints out the 3-point average of the Imp. Thermodynamic
	contributions.

OpTherm.dat: Prints out the expectation value of the Imp. Operator(s).

therm.data: This files contain the Thermodynamics of the Fermionic
	chain (FC) by itself. The program Open_Files() reads this file and subtract
	off the contributions from the bath.

======================================================================

Definition of Matlab functions and sub-routines:
------------------------------------------------

NRG.m : Contains the main function NRG().

Input Variables.m : Declares the set of parameters for the NRG run.

Parms.m : (pre-)Declares the various variables/Parameters  used for the run.

Open_Files.m: Subroutine that opens all the necessary files for the NRG run.

Iteration0: Set up and Diagonalize the Hamiltonian for Iteration 0 of the run.

Iteration.m: Set up and Diagonalize the Hamiltonian for Iterations > 0. 

ITend.m: End of Iteration. Print out NRG energies (measured from the ground
	state), Thermodynamics and <Op>.

Plot_Figures.m: Plot (using Matlab's plotting tool) the Energy flow diagrams
	for even and odd iterations. Also plots the various thermodynamic
	quantities calculated against temperature.

Close_Files.m: Close all the files that were opened during the run.

Impose.m: Imposes an artifical degeneracy over closely lying energy levels
    (\delta E < SMALL).

Esort.m: Sorts the eigenstates with respect to Energy. In addition to energy,
	also sorts the matrix c(containing eigenvectors) and Kode(containing
	Qtm. Nos.).

===========================================================================

Sample output files and plots for the Fermionic Chain, Kondo and Andersoncan
        be found at the './Sample_Outputs'

Last edited: 30.09.2020.

