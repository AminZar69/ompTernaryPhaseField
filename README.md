# ompTernaryPhaseField
"ompTernaryPhaseField" is a multiphase flow solver based on the lattice Boltzmann method which has been speeded up using the OpenMP API. The code models
the 2D case of two immicible droplets on a flat surface surrounded by a third ambient immiscible phase. Using this code and adding suitable boundary conditions, more complex problems such as ternary flow in porous media can be simulated. 
### Requirements
The gFortran compiler is the only requirement for running this package. The "-fopenmp" flag has been added for compilation and linking in the Makefile. The output files are generated in the ascii VTK format readable by the paraview which is an open-source visualisation package found via <https://www.paraview.org/download/>.

### Build and execution by passing the desired number of threads (n)
    make clean
	make NUM_THREADS = n run

### Build and execution when NUM_THREADS is hard-coded
    make clean
    make 
	cd bin
	./Ternary_Contact_Angle_On_Flat_Surface





   