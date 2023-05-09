/**********************************************************************
   UDFs for kinetic switch, mixture viscosity, DPM injections
***********************************************************************

Preamble:
Author: Torben Frey (+49 40 42878-4124, torben.frey@tuhh.de)

General Info for DPM injections:
- This UDF clusters particles that precipitate from the Eulerian phase into an injection file. 
- A particle is generated if the (mass) concentration of a component exceeds a solubility limit.
- The injection file can be used in the DPM model as discrete Lagrangian particles.
- It is recommended to use the species-transport (and volumetric reaction) model with this UDF.

You input is required according to your model:
1) [ln 42] Give fluid domain ID (define/boundary-conditions/list-zones)
2) [ln 51] Optional: change the name of the UDF
2) [ln 54] Specify molecular weight of precipitating species [kmol/kg]
3) [ln 55] Give a formation rate threshold for precipitation [kmol/m3/s]
4) [ln 56] Give the density of the precipitated particle in [kg/m3]
5) [ln 57] Give a radius for clustering. All cells inside this radius will be consolidated into one injection.
6) [ln 58] Give a minimum number of cells per cluster. Smaller clusters will be ignored.
7) [ln 63] Optional: Change the name of the injection file to be written.

How to use this UDF:
1) Read case & data file (Steady, single phase, species-transport-model, at least 1 reaction)
2) Compile and load UDF:
	TUI: define/user-defined/compiled-functions compile "libudf" yes "particle-precipitation.c" "" ""
	TUI: define/user-defined/compiled-functions load "libudf"
    --> if UDM error message occurs: unload and reload UDF
3) Execute UDF on demand (this will write the injection file):
	TUI: define/user-defined/execute-on-demand "particle-precipitation::libudf"
4) Read case & data file (Steady, DPM-model)
5) Create injection from file (change name according to ln 63)
	define/models/dpm/injections/create-injection UDF-injection no yes file no "precipitating-particle-file.inj" no no no no
6) Set DPM model specifics as needed

Have fun! */

#include "udf.h"

/* find fluid_id of interior (here 170), TUI: define/boundary-conditions/list-zones */
# define FLUID_ID 170 

/* static real are defined for all following UDFs */
static real totalX[ND_ND], totalU[ND_ND], totalD, totalM; static int udm_offset = UDM_UNRESERVED;

/**********************************************************************
   UDF writes injection file for steady particle tracking
***********************************************************************/

DEFINE_ON_DEMAND(particle-precipitation) {

/* define custom parameters for clustering */
		real part_mol_wt = 100; /* molecular weight [kmol/kg] of precipitating species */ 
		real prec_thresh = 1e-6; /* threshold of formation rate [kmol/m3/s] for clustering */ 
		real part_densit = 1000; /* specify particle density [kg/m3] (important for mass flow and particle size) */
		real radi_thresh = 0.001; /* sphere radius [m] around concentration maximum for clustering */
		real num_cluster = 10; /* minimum number of cells in cluster */

/* Specify injection filename */
		#if !RP_NODE	
			FILE *fp1 = NULL;
			char filename[]="precipitating-particle-file.inj";
		#endif

/* Do not change any code below */

		/* Different variables are needed on different nodes */
		#if !RP_HOST
			Domain *domain=Get_Domain(1);
			Thread *thread=Lookup_Thread(domain,FLUID_ID);
			cell_t c; 
		#else
			int i;
		#endif


		int size, pe;  /* data passing variables */
		real *arrayx, *arrayy, *arrayz, *arrayu, *arrayv, *arrayw, *arraym, x[ND_ND];
		real max_loc[ND_ND], max_species;
		int max_index, j, l;

		#if !RP_NODE /* HOST */
			if ((fp1 = fopen(filename, "w"))==NULL)
				Message("Warning: Unable to open %s for writing!\n",filename);
			else
				Message("Writing %s as injection file . . .\n",filename);
		#endif

		/* UDF Now does 2 different things depending on NODE or HOST */
		#if RP_NODE
			/* Each Node loads up its data passing array */
			size = THREAD_N_ELEMENTS_INT(thread);
			arrayx = (real *)malloc(size * sizeof(real)); arrayy = (real *)malloc(size * sizeof(real)); arrayz = (real *)malloc(size * sizeof(real));
			arrayu = (real *)malloc(size * sizeof(real)); arrayv = (real *)malloc(size * sizeof(real)); arrayw = (real *)malloc(size * sizeof(real));
			arraym = (real *)malloc(size * sizeof(real));

			begin_c_loop_int(c,thread) {	/* create field variable arrays for injection data*/
					
				if (C_UDMI(c,thread,0) > prec_thresh) { 
					C_CENTROID(x,c,thread);
					arrayx[c] = x[0]; /* x-position */
					arrayy[c] = x[1]; /* y-position */
					arrayz[c] = x[2]; /* z-position */
					arrayu[c] = C_U(c,thread); /* x-velocity */
					arrayv[c] = C_V(c,thread); /* y-velocity */
					arrayw[c] = C_W(c,thread); /* z-velocity */
					/* mass flow = formation-rate * Volume / molar mass */
					arraym[c] = C_UDMI(c,thread,0) * C_VOLUME(c,thread) * part_mol_wt;
				} 
			} end_c_loop_int(c,thread)
		
			/* Set pe to destination node */
			/* If: on node_0 send data to host, Else: send to node_0 because compute nodes connect to node_0 & node_0 to host */
			pe = (I_AM_NODE_ZERO_P) ? node_host : node_zero; 

			PRF_CSEND_INT(pe, &size, 1, myid); /* send arrays on compute nodes to node_0 */
			PRF_CSEND_REAL(pe, arrayx, size, myid); PRF_CSEND_REAL(pe, arrayy, size, myid); PRF_CSEND_REAL(pe, arrayz, size, myid);
			PRF_CSEND_REAL(pe, arrayu, size, myid); PRF_CSEND_REAL(pe, arrayv, size, myid); PRF_CSEND_REAL(pe, arrayw, size, myid);
			PRF_CSEND_REAL(pe, arraym, size, myid);
			free(arrayx); free(arrayy); free(arrayz);
			free(arrayu); free(arrayv); free(arrayw);
			free(arraym);

			/* node_0 now collect data sent by other compute nodes and sends it straight on to the host */
			if (I_AM_NODE_ZERO_P) 
				compute_node_loop_not_zero (pe) {
					PRF_CRECV_INT(pe, &size, 1, pe); /* receive arrays from compute nodes */
					arrayx = (real *)malloc(size * sizeof(real)); arrayy = (real *)malloc(size * sizeof(real)); arrayz = (real *)malloc(size * sizeof(real));
					arrayu = (real *)malloc(size * sizeof(real)); arrayv = (real *)malloc(size * sizeof(real)); arrayw = (real *)malloc(size * sizeof(real));
					arraym = (real *)malloc(size * sizeof(real));
					PRF_CRECV_REAL(pe, arrayx, size, pe); PRF_CRECV_REAL(pe, arrayy, size, pe); PRF_CRECV_REAL(pe, arrayz, size, pe);
					PRF_CRECV_REAL(pe, arrayu, size, pe); PRF_CRECV_REAL(pe, arrayv, size, pe); PRF_CRECV_REAL(pe, arrayw, size, pe);
					PRF_CRECV_REAL(pe, arraym, size, pe); 
					/* send arrays on node_0 to host */
					PRF_CSEND_INT(node_host, &size, 1, myid); 
					PRF_CSEND_REAL(node_host, arrayx, size, myid); PRF_CSEND_REAL(node_host, arrayy, size, myid); PRF_CSEND_REAL(node_host, arrayz, size, myid);
					PRF_CSEND_REAL(node_host, arrayu, size, myid); PRF_CSEND_REAL(node_host, arrayv, size, myid); PRF_CSEND_REAL(node_host, arrayw, size, myid);
					PRF_CSEND_REAL(node_host, arraym, size, myid);
					free((char *)arrayx); free((char *)arrayy); free((char *)arrayz);
					free((char *)arrayu); free((char *)arrayv); free((char *)arrayw);
					free((char *)arraym);
				}
		#endif /* RP_NODE */

		#if RP_HOST
			l = 0; /* l is number of injection */
			compute_node_loop (pe) { /* only acts as a counter in this loop */
				/* Receive data sent by each node and write it out to the file */
				PRF_CRECV_INT(node_zero, &size, 1, node_zero); 
				arrayx = (real *)malloc(size * sizeof(real)); arrayy = (real *)malloc(size * sizeof(real)); arrayz = (real *)malloc(size * sizeof(real)); 
				arrayu = (real *)malloc(size * sizeof(real)); arrayv = (real *)malloc(size * sizeof(real)); arrayw = (real *)malloc(size * sizeof(real));
				arraym = (real *)malloc(size * sizeof(real));
				/* receive data from node_0 */
				PRF_CRECV_REAL(node_zero, arrayx, size, node_zero); PRF_CRECV_REAL(node_zero, arrayy, size, node_zero); PRF_CRECV_REAL(node_zero, arrayz, size, node_zero);
				PRF_CRECV_REAL(node_zero, arrayu, size, node_zero); PRF_CRECV_REAL(node_zero, arrayv, size, node_zero); PRF_CRECV_REAL(node_zero, arrayw, size, node_zero);
				PRF_CRECV_REAL(node_zero, arraym, size, node_zero);
				/* cluster and write injection arrays */
				/* note: clustering runs in the compute_node_loop. It is NOT a global clustering but only on each node */
				/* find first maximum of species mass flow */
				max_species = 0.0; max_index = 0; 
				for (i=0; i<size; i++) {
					if (arraym[i] > max_species) {
						max_loc[0] = arrayx[i]; max_loc[1] = arrayy[i]; max_loc[2] = arrayz[i];
						max_species = arraym[i]; max_index = i;
					}
				}
				while (max_species != 0) {
					/* create (one) cell cluster with spherical support domain (r) from maximum concentration */
					totalX[0] = 0.0; totalX[1] = 0.0; totalX[2] = 0.0;
					totalU[0] = 0.0; totalU[1] = 0.0; totalU[2] = 0.0;
					totalM = 0.0; j = 0; /* j is the number of clustered cells */
					for (i=0; i<size; i++) {
						if (sqrt(pow(arrayx[i]-max_loc[0],2)+pow(arrayy[i]-max_loc[1],2)+pow(arrayz[i]-max_loc[2],2)) < radi_thresh) {
							totalU[0] = totalU[0] + arrayu[i]; */
							totalU[1] = totalU[1] + arrayv[i]; */
							totalU[2] = totalU[2] + arrayw[i]; */
							totalM = totalM + arraym[i]; arraym[i] = 0.0;
							j = j + 1;
						}
					} 
					totalX[0] = max_loc[0]; totalX[1] = max_loc[1]; totalX[2] = max_loc[2];
					totalU[0] = totalU[0]/j; totalU[1] = totalU[1]/j; totalU[2] = totalU[2]/j;
					totalD = pow(6*totalM/part_densit/M_PI,0.33333333);
					l = l + 1;
					if (j >= num_cluster) {
						/* (( x y z u v w diameter temperature mass-flow) name ) */
						fprintf(fp1,"((%e %e %e %e %e %e %e %i %e) injection_%i)\n",totalX[0],totalX[1],totalX[2],totalU[0],totalU[1],totalU[2],totalD,80,totalM,l);
						Message("Injection %i written at (%.2e, %.2e, %.2e) with m = %.2e kg/s and d = %.2e m (%i cells clustered).\n",l,totalX[0],totalX[1],totalX[2],totalM,totalD,j);
					}
					else {
						Message("Injection %i not written - only %i singular cells.\n",l,j);
					}	
					/* find next maximum of species mass flow */
					max_species = 0.0; max_index = 0;
					for (i=0; i<size; i++) {
						if (arraym[i] > max_species) {
							max_loc[0] = arrayx[i]; max_loc[1] = arrayy[i]; max_loc[2] = arrayz[i];
							max_species = arraym[i]; max_index = i;
						}
					}
				}
				free(arrayx); free(arrayy); free(arrayz);
				free(arrayu); free(arrayv); free(arrayw);
				free(arraym);
				Message("Clustering on Node %i done! \n",pe+1);
			}
			fclose(fp1); /* Close the file */
		#endif
}

/**********************************************************************
   UDF reserves 1 UDM for libudf and renames the UDM to enhance 
   postprocessing and particle injection UDF
***********************************************************************/

#define NUM_UDM 1

DEFINE_EXECUTE_ON_LOADING(set_alphalyon_formation_rate_as_UDMI, libname) {
	if(udm_offset == UDM_UNRESERVED) {
		udm_offset = Reserve_User_Memory_Vars(NUM_UDM);
	}
	if(udm_offset == UDM_UNRESERVED) {
		Message("You need to define %d extra UDM in GUI and then reload current library %s \n", NUM_UDM,libname);
	}
	else {
		Message("%d UDM has been reserved by the current library %s \n", NUM_UDM,libname);
		Set_User_Memory_Name(udm_offset,"reaction-1-formation-rate");
	}
	Message("UDM Offset for current loaded Loaded Library = %d \n", udm_offset);
	
	/**********************************************************************
	This section initializes C_UDMI to zero in all cells
	***********************************************************************/
	Domain *d;
	Thread *thread;
	cell_t cell;
	d = Get_Domain(1);
	
	if(udm_offset != UDM_UNRESERVED) {
		thread_loop_c(thread,d) {
			begin_c_loop(cell,thread) {
				C_UDMI(cell,thread,0) = 0;
			} end_c_loop(cell,thread)
		}
	}
	else {
		Message("UDM has not yet been reserved for library \n");
	}
}

/*********************************************************************
Generic arrhenius reaction rate for reaction "reaction-1"
**********************************************************************/

DEFINE_VR_RATE(arrhenius-reaction,cell,thread,r,mw,yi,rr,rr_t) 
	{
    real ci, prod;
	int i;
	if (!strcmp(r->name, "reaction-1"))
	{
		/* Calculate Arrhenius reaction rate for reaction-1 */
		prod = 1.;
		for(i = 0; i < r->n_reactants; i++)
		{
			/* ci = density * concentration * mol weight [kmol/m3] */
			ci = C_R(cell,thread) * yi[r->reactant[i]] / mw[r->reactant[i]];
			/* prod = c1^exp1 * c2^exp2 * ... * cn^expn */
			prod *= pow(ci, r->exp_reactant[i]);
		}
		/* rr = A * T^b * e^(-E/RT) * c1^exp1 * c2^exp2 * ... * cn^expn */
		/* Fluent will choose the smaller of rr and rr_t as reaction rate */
		*rr = r->A * exp(- r->E / (UNIVERSAL_GAS_CONSTANT * C_T(cell,thread))) * pow(C_T(cell,thread), r->b) * prod;
		*rr_t = *rr;
    
        /* set reaction rate in user defined memory slot 0: C_UDMI(cell,thread,0) */
		if(udm_offset != UDM_UNRESERVED){
			C_UDMI(cell,thread,0) = r->A * exp(- r->E / (UNIVERSAL_GAS_CONSTANT * C_T(cell,thread))) * pow(C_T(cell,thread), r->b) * prod;
		}
		else {
			Message("UDM has not yet been reserved for library 1\n");
		}
	}
	else
	{
		Message("Unknown Reaction\n");
	}
}