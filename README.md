# TASEP
Solver for one dimensional TASEP Systems

This Software is designed and maintained by Kilian Pioch at the University of Bayreuth. 

It's main purpose is the solving of different types of one dimensional TASEP systems with independent jump processes.

Currently three distinct models can be simulated. 

1. Master Equation
2. Mean Field Models of arbitrary order (this includes the Ribosome Flow Model as it is the Mean Field Model of order 1)
3. Exact Steady State Calculation, using the Matrix Algebra developed by Derrida, Evans, Blythe et al.

Additional information on the deriavtion of the higher order mean field models can be found at https://arxiv.org/abs/2402.10234.



## Dependencies

Rudimentary Parallelization is implemented using OpenMP, see https://www.openmp.org/ and the installation instructions therein.\
The Steady State calculations require GnuMP , see https://gmplib.org/ and the installation instructions therein

## Installation

Currently the Software is confirmed to reliably compile on linux systems with make and gcc.

Compilation is handled by make. 
Additional parameters can be passed to make for different compilation flavours: 

1. make warn=1:   Additional Warnings
2. make debug=1:  Include Debug Symbols when compiling
3. make asan=1:   Compile for use with Adress Sanitizer (use in combination with debug=1)
4. make no-gmp=1: Don't compile executables reliant on GNU-MP


# Usage
The following will explain the usage and functionality of the provided programs.

In general all programs are called analogously. They mainly differ in their underlying models. 

### masterEquation

   This program uses RK4 to solve the [Master Equation](https://en.wikipedia.org/wiki/Master_equation)
   of the TASEP-System which is specified in the command line options. 
   
   The master equation describes an ODE of the temporal evolution of the probability of all full length lattice states. The trajectories of these states are then used to calculate the individual expected values for the simulation time. 
   This top down approach is the most accurate to-date.
   Note that due to the exponential growth of the dimensionality of the ODE, Lattices of size >31 will require upwards of 128GiB of RAM.

   A typical invocation my look like this: 

   `bin/masterEquation -i <jumprate file> -o <output file> -n <lattice size> -t <start time> -T <end time> -s <step number> -f --thread-num 4`

   The command line options are used as follows: 
   
   - `-i <jumprate file>` this switch enables reading of n+1 jump rates from some file <jumprate file>. Entries in the file have to be white space separated.

   - `-o <output file>` this switch points the program to the output file. Additional files and folders with similar names may be created to save different values. If this switch is not present a file called "generic_output_vec" will be created in the directory from which the program is run.

   - `-n <lattice size>` this switch specifies the lattice size

   - `-t <start time>` this switch specifies start time. Note: The actual value here does not have any effect on the simulation. 

   - `-T <end time>` this switch specifies end time. Note: The actual value here dies not have any effect on the simulation and that the program may exit earlier than the time specified here. During the solving process the solver keeps track of the rate of change of the ODE system, if the absolute change of all states drops below $`10^{-10}* h`$ for 100 consecutive iterations, the solver exits.

   - `-s <step number>` this switch specifies the number of timesteps. In combination with `-t` and `-T` this switch defines the temporal resolution of the simulation. 
   Reccomendended values are in the range [(T-t) * 10,(T-t) * 1000].

   - `-f` this switch enables a much more memory intensive, but much quicker compilation

   - `--thread-num 4` this switch tells openMP to use at most 4 threads

In Addition to these, there are several other options avaliable. 


Jump Rates:
   - `--alpha`, `--beta`, `--gamma` these switches are used to set rates for entry flow (alpha), exit flow (beta) and lattice flow (gamma). They can act as a substitute for `-i <input_file>`

   - `--save-select [a|e|p|c|m]` this switch can control which values will be saved to file during or after the calculation (depeding on `-a <num>`). Several can be selected simultaneously.
   Possible options are:
      `-a` All of the below\
      `-e` Expected values ordered from entry to exit node\
      `-p` Probabilities of all possible states ordered from the state represented by binary 0 to state $`2^N - 1`$, with N=chain length\
      `-c` Covariance Matrices\
      `-m` Multiplicative expected values (depreciated)\

   The solver is also able to solve systems with time dependend jump rates.  
   It is important to note here that in this case the default behaviour of the solver is to run until `-T <num>` is reached. Because of entrainment it is not expected for the trajectories to settle at some steady state

   Necessary switches for time variant flow rates

   - `--time-var` this enables time variant simulation
   with only this switch selected, the solver will assume that every jump rate is of the form $`sin(1 * t)`$

   Optional switches for time variant flow rates

   - `--function-type <filename>` program tries to read from given filename. The entries in the file will set the function types for the hop rates from left to right. The contents of the file may look like this "sscls" this would give:
      - $`\alpha(t) = sin(t)`$
      - $`h_0(t) = sin(t)`$
      - $`h_1(t) = cos(t)`$
      - $`h_2(t) = t`$
      - $`\beta(t) = sin(t)`$

   - `--amplitude <filename>` program tries to read from given filename. The content of the file will set the amplitude in the functions for the jump rates. For a given file with the content $`4\quad -1.3\quad 1\quad 0\quad 1.1`$ this would give:
      - $`h_0(t) = \alpha(t) = 4sin(t)`$
      - $`h_1(t) = -1.3sin(t)`$
      - $`h_2(t) = 1cos(t)`$
      - $`h_3(t) = 0t`$
      - $`h_4(t) = \beta(t) = 1.1sin(t)`$

   - `--period <filename>` program tries to read from given filename. The content of the file will set the period in the functions for the jump rates. For a given file with the content $`1\quad 2\quad 4\quad 2\quad 3`$ this would give:
      - $`h_0(t) = \alpha(t) = 4sin(1t)`$
      - $`h_1(t) = -1.3sin(2t)`$
      - $`h_2(t) = 1cos(4t)`$
      - $`h_3(t) = 0t`$ (ignored in this case)
      - $`h_4(t) = \beta(t) = 1.1sin(3t)`$
   
   - `--offset <filename>` program tries to read from given filename. The content of the file will set the offset in the functions for the jump rates. For a given file with the content $`1\quad -0.5\quad 0\quad 0.3\quad 1`$ this would give:
      - $`h_0(t) = \alpha(t) = 4sin(1t+1)`$
      - $`h_1(t) = -1.3sin(2t-0.5)`$
      - $`h_2(t) = 1cos(4t)`$
      - $`h_3(t) = 0t`$ (ignored in this case)
      - $`h_4(t) = \beta(t) = 1.1sin(3(t+1))`$

   - `--constant <filename>` program tries to read from given filename. The content of the file will set the offset in the functions for the jump rates. For a given file with the content $`-1\quad -0.5\quad 0\quad 0.5\quad 1`$ this would give:
      - $`h_0(t) = \alpha(t) = 4sin(1t+1)-1`$
      - $`h_1(t) = -1.3sin(2t-0.5)-0.5`$
      - $`h_2(t) = 1cos(4t)+0`$
      - $`h_3(t) = 0t+0.5`$
      - $`h_4(t) = \beta(t) = 1.1sin(3(t+1))+1`$

### mf_full

   This program utilises a system description obtained by expanding the ODEs for individual expected values. 
   The construction of the ODE can be found at the arxiv link given at the beginning of this document.

   A typical execution my look like this: 

   `bin/mf_full -i <input> -o <output> -n <size> -c <order> -t 0 -T <end_time> -s <steps> -v`

   The command line options are used as follows: 
   
   - `-i <input>` this switch points the program to the file "file_a" from which it will attempt to read n+1 jump rates.
   Entries in the file have to be white space separated.

   - `-o <output>` this switch points the program to the file "file_b" into which it will attempt to write the results of the comutation. 
   Additional files and folders with similar names may be created to save different values. If this switch is not present a file called "generic_output_vec" will be created in the directory from which the program is run.

   - `-n <size>` this switch specifies the chain length

   - `-t 0` this switch specifies start time. Note: The actual value here does not have any effect on the simulation. 

   - `-T <end_time>` this switch specifies end time. Note: The actual value here dies not have any effect on the simulation and that the program may exit earlier than the time specified here. During the solving process the solver keeps track of the rate of change of the ODE system, if the absolute change of all states drops below $`10^{-5}*\delta h`$ for 100 iterations, the solver exits.
   !!This should be the default behaviour!!

   - `-s <steps>` this switch specifies the number of timesteps. In combination with `-t` and `-T` this switch defines the temporal resolution of the simulation. 
   Reccomendended values are in the range [(T-t) * 10,(T-t) * 1000].


   Several other command line options can be found by runnign the program using the `-h` switch.



### steadystates_gmp

   Similarily to the previous program, this will calculate the steady states of chains with homogenious interior jump rates and variable entry/exit hop-rates.

   The main difference is the utilisation of the GNU multi precision library (GMP), which significantly improves accuracy of the results. 

   A typical invocation of the program may look as follows:

   `bin/steadystates_gmp -a <alpha> -b <beta> -n 10`

   This will print the results to stdout. 
   If desired, the results can be printed to a file specified by the switch:

   `-o --output <filename>`

   One functionality of GMP is the variable bit-width of the datatypes used. 
   If you know beforehand, that the simulation that will be run will need more than 64 bits, this can be specified by:

   `-p --precision <int value>` This will set the starting bit width to the value specified. 
   GMP will, however widen the used datatypes by itself, if numerical inaccuracies are detected internally. This is not a feature of `steadystates_gmp` but GMP itself.

   If the results should be printed with more or less accuracy, the printing accuracy can be specified by:

   `--pp --print-precision <int value>` This will set the number of digits printed to the specified output stream
