# particle-precipitation

DOCUMENTATION for particle-precipitation.c

Author: Torben Frey

The file "particle-precipitation.c" is developed for ANSYS Fluent v2021r2 on Linux and Windows.

"particle-precipitation.c" is a compiled User Defined Function (UDF) and generates a particle file
for a Euler-Lagrangian two-phase model, the discrete phase model (DPM) in ANSYS Fluent.

Input:
Converged (steady) single-phase solution with at least 1 reaction product
particle-precipitation.c with user input

Output: 
Array of precipitating particles (ANSYS Fluent injection file) of the form:
(( x-location y-location z-location u-velocity v-velocity w-velocity diameter temperature mass-flow) name/number )



INSTRUCTIONS
1) How to use particle-precipitation.c with ANSYS Fluent

1.1) Read case & data file (Steady, single phase, species-transport-model, at least 1 reaction)
    steady: not required, but recommended to start with minimalistic example
    single phase: continuous Eulerian phase, not yet tested with multiphase models
    species-transport-model: one species (the product of a reaction) triggers the precipitation

1.2) set up UDF: particle-precipitation.c 
    [line 42] enter fluid domain ID (find in ANSYS Fluent TUI: "define/boundary-conditions/list-zones")
    [line 54] specify molecular weight of precipitating species [kmol/kg]
    [line 55] give a formation rate threshold for precipitation [kmol/m3/s]
    [line 56] give the density of the precipitated particle in [kg/m3]
    [line 57] give a radius for clustering. All cells inside this radius will be consolidated into one injection
    [line 58] give a minimum number of cells per cluster. Smaller clusters will be ignored
    [line 63] enter the name of the injection file, the UDF will overwrite existing files

1.3) Save particle-precipitation.c into working directory with *.cas and *.dat files

1.4) In the ANSYS Fluent GUI - User-defined tab, add 1 UDM slot

1.5) Compile and load UDF:
	TUI: define/user-defined/compiled-functions compile "libudf" yes "particle-precipitation.c" "" ""
	TUI: define/user-defined/compiled-functions load "libudf"
    --> if UDM error message occurs: unload and reload UDF, troubleshoot with log file

1.6) Execute UDF on demand (this will write the injection file):
	TUI: define/user-defined/execute-on-demand "particle-precipitation::libudf"

Depending on the size of your domain and the UDF-specific inputs, the injection file generation might take a while.


2) How to use the injection file *.inj in a two phase simulation

2.1) Read case & data file (Steady, DPM-model)

2.2) Create injection from file (change name of file accordingly)
	TUI: define/models/dpm/injections/create-injection UDF-injection no yes file no "precipitating-particle-file.inj" no no no no

2.3) Set up DPM model specifics 

2.4) Run simulation (Tip: simplify models, e.g. compute only conti and momentum eqs)


3) Troubleshooting 

3.1) ANSYS Fluent does not load compiled UDF
    - check if UDF is in working directory with *.cas and *.dat files
    - use Buili-in compiler only with Windows OS
    - delete folder libudf and recompile UDF
    - check log file for "errors" during compilation (warnings are ok), especially syntax errors

3.2) "UDM has not yet been reserved for library" or "You need to define 1 extra UDM in GUI and then reload current library"
    - in the ANSYS Fluent GUI User-defined tab, set 1 UDM slot and reload UDF

3.3) ANSYS Fluent crashes when executing UDF on demand
    - check if *.dat file is loaded. UDF needs data in cells to work.
