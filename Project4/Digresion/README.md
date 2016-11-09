##A small digresion
The files here are implemented using the CUDA library.
The program makes a visualization of how an initial configuration at the ground state
changes per 2000 MC-cycle for a given temperature and number of spins.
The number of spins must be a multiple of 16. During the process of creating this visualization, a number of 256 spins
proved to be decent.

* To run the visualization, run ising.exe found in ising/Debug.
  The executable will not work if moved to another folder or all the files in in this folder are not cloned
  For now, Microsoft Visual Studio 2013 Runtime Libraries must be installed to run the .exe file,otherwise the system will not find         MSVCP120D.dll This is under investigation how this can be fixed.
