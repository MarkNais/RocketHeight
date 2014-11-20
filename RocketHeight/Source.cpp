/***********************************************************************************
* Mech 7171 Engineering Programing
* LAB #4
*
* Purpose: The of this lab is to solve the heat conduction equation in 2D
*          using finite difference (F-D) methods.  The solution will be compared
*          to an analytic solution and solution convergence will be tracked.
*          This lab will exercise a facets of C taught in Mech 7171.
*
* Author: Corbin Turner, Cyrus Ang, Mark Naismith
* Student ID: A00780890, A00781218, A00819714
*
* Declaration:
* I, Corbin Turner, Cyrus Ang, Mark Naismith, declare that the following
* program was written by us.
*
* Date Created: Dec.6/2013
* Revised: Mark Naismith, 2/25/2014 - Edited for Fluid Flow Project
* Revised: Mark Naismith, 9/20/2014 - Edited for Koorosh's Fin Lab 01
***********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#pragma warning(disable:4996)

#define PI                       (4.0*atan(1.0))
#define MAX_ITER                 100000            // maximum iterations for F-D
#define FULL_THERMAL_MEAN        0.99959            // Temperature is fully developed here
#define DATA_FILENAME			 "data.txt"
#define MAX_LINE_LENGTH			 1024
#define MASS_ROCKET				 0.035
#define FULL_ENGINE				 0.080
#define EMPTY_ENGINE			 0.040

//  structure for finite-difference solution
typedef struct
{
	double y;      //Fully developed flow velocity
	double v;
	double a;
	double t;
	double Fr;
	double Fd;
}
PLATEPOINT;


typedef struct
{
	int scase;		//Case counter (which simulation is it?)
	int Ny;			//Node count in y
	int tfinal;	//When to stop the simulation - The Theta that determines when the simulation ends.
	double dt;
	int iter;		//Tracking the iteration count
	int Rocket;

	//Interior limits of the simulation, ex: of a 9 by 5 grid, only the inside 7 by 3 grid is calculated.
	int NyInter;
	int run;		//Run the simulatuon? 1=run 0=dont run

	PLATEPOINT *r;		//The pointer to the first  array of plate points. (t) "current"

}
PROGRAMDATA;

//------------------------- FUNCTION PROTOTYPES -----------------------------------------

// Get data from text file
PROGRAMDATA GetProgramData(FILE*);

//Initialize some of the variables inside of the pd struct.
PROGRAMDATA pdInit(PROGRAMDATA);

// Function to contain the specific order of regions to be simulated
void num_simulation(PROGRAMDATA);

// functions for each simulation
void simulate(PROGRAMDATA);

// allocating and free allocation fucntions
PROGRAMDATA allocate(PROGRAMDATA);
void freepp(PROGRAMDATA);

// reading and writing functions
FILE *filewrite(PROGRAMDATA);

//rounding function
int nint(double);
//------------------------- END OF FUNCTION PROTOTYPES ----------------------------------


int main()

{
	int i;  //loop counter


	//header
	printf("MECH 7171 Engineering programming                %c\n", 179);
	printf("Lab #4                                           %c\n", 179);
	printf("Mark Naismith, A00819714, Set 5B                 %c\n", 179);
	printf("Cyrus Ang,     A00781218, Set 5B                 %c\n", 179);
	printf("Corbin Turner, A00780890, Set 5B                 %c\n", 179);
	printf("December 6, 2013                                 %c\n", 179);
	printf("This program reads data from a data.txt file and %c\n", 179);
	printf("stores the data in a structure within the program%c\n", 179);
	printf("that is transfered between modular functions.    %c\n", 179);
	printf("The tasks that this program completes are the    %c\n", 179);
	printf("analytical solution to 2D conductive             %c\n", 179);
	printf("heat transfer probrlems as well as               %c\n", 179);
	printf("the numerical solutions. The end result is a     %c\n", 179);
	printf("given heat transfer problem                      %c\n", 179);
	printf("that is solved numerically.                      %c\n", 179);

	//this loop "closes the box" of the header. Run program for details.
	for (i = 0; i<49; i++)printf("%c", 196);
	printf("%c\n", 217);

	//creation of the structure for simulations. Holds all relevant data
	PROGRAMDATA pd;

	//File stream pointer for the input file. This will be passed to the GetProgramData
	//function for every simulation
	FILE *f;

	//Opens stream to the input file. Read only
	f = fopen(DATA_FILENAME, "r");

	//Check to see if the file exists, if not, exits function, main() complains
	if (f == NULL)
	{
		printf("Input file \"%s\" doesn't exist! Exiting...", DATA_FILENAME);
		getchar();
		return 5;
	}

	i = 0;
	//The proceding blocks are all very similar. Only the first will be explained.
	do{
		i++;
		pd = GetProgramData(f);
		if (pd.run == 1)//Checks to see if user wants to run 
		{
			printf("\nRunning Simulation %d :\n", i);
			simulate(pd);
		}
	} while (pd.run != -1);

	// closes data file once all simulations are complete 
	//(could be closed before Given problem)
	fclose(f);

	//End Program
	printf("Press Enter to end the program....");
	fflush(stdin);//needed if the last file needed to be overwriten.
	getchar();
	return 0;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                          FUNCTIONS
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/********************************************************
* GetProgramData *
* *
* Purpose: Retrieves only the appropriate amount of data from the input file.
* Note: it does not error check for bad lines or an incomplete file.
* *
* Parameters *
* i - loop counter for retrieving data ,
* j - static int to determine case for filewrite function
* *pdp[] - array of pointers to addresses of varibales in PROGRAMDATA structure
* *    ^good idea Dave^
* *
* Returns *
* pd to PROGRAMDATA function  *
********************************************************/
PROGRAMDATA GetProgramData(FILE *f)
{
	//create a new program data structure to overwrite the old one.
	PROGRAMDATA pd;

	char buff[MAX_LINE_LENGTH];   //Receiver of lines from the input file
	int i;                    //loop counter  
	static int j = 1;           //case counter (how many times 
	//has this function been called?)
	//An array of pointers where the read data is to be stored.
	int *pdpi[] = { &pd.run, &pd.Rocket, &pd.tfinal };
	double *pdp[] = { &pd.dt};

	//assume the simulations should not run, Minor protection from read error.
	pd.run = 0;

	// retrieves 3 more lines of 'int' data, and points them to their destination
	for (i = 0; i<3; i++)
	{
		fgets(buff, MAX_LINE_LENGTH, f);
		if (feof(f) != 0){
			pd.run = -1;
			return pd;
		}
		*pdpi[i] = atoi(buff);
	}//End of for

	// retrieves 2 more lines of 'double' data, and points them to their destination
	for (i = 0; i<1; i++)
	{
		fgets(buff, MAX_LINE_LENGTH, f);
		if (feof(f) != 0){
			pd.run = -1;
			return pd;
		}
		*pdp[i] = atof(buff);
	}//End of for

	// if the end of the file is reached, stop trying to run simulations


	pd.scase = j; //sets the case for filewrite later
	j++;        //increments j
	fgets(buff, MAX_LINE_LENGTH, f); // skips a line down in data file

	// since we have not used fclose, the file is still open and the line that this
	// function is looking at remains the same
	return pd;     //return pd to main()   
}//End of GetProgramData

/********************************************************
* allocate *
* *
* Purpose: allocate an array to be used for simulations
* *
* Parameters *
* w,h - loop counters for nodes horizontal, vertical
* W,H - horizontal and vertical nodes calculated
*       from input data *
* *
* Returns *
* Initialized allocated array into PLATEPOINT structure nested in PROGRAMDATA *
********************************************************/
PROGRAMDATA allocate(PROGRAMDATA pd)
{
	int i;
	pd.Ny = (double)pd.tfinal / pd.dt + 1;
	pd.NyInter = pd.Ny - 1;

	if (pd.Ny<3)
	{
		printf("Sorry, you must use a positive integer of 3 or larger");
		printf("for Node Count\n");
		printf("\nPress Enter to end this program...");
		getchar();
		exit(0);
	}

	// allocate memory for a horizontal array of pointers
	pd.r = (PLATEPOINT*)malloc(pd.Ny*sizeof(PLATEPOINT));

	// check allocation for array of pointers
	if (pd.r == NULL)
	{
		printf("Cannot allocate pd.r, exiting program...\n");
		getchar();
		exit(0);
	}

//	//introduce allocation for the new freepp additions.
//	for (i = 0; i<8; i++){
//		*pda[i] = (double*)malloc(pd.Ny*sizeof(double));
//		if (*pda[i] == NULL)
//		{
//			printf("Cannot allocate matrix-array number: %d, exiting program...\n", i);
//			getchar();
//			exit(0);
//		}
//	}
	return pd;
}

/********************************************************
* simulate1 *
* *
* Purpose: Run const. temp coarse simulation
* *
********************************************************/
void simulate(PROGRAMDATA pd)
{
	pd = allocate(pd);
	pd = pdInit(pd);
	num_simulation(pd);
	freepp(pd);
}

/********************************************************
* pdInit *
* *
* Purpose: Set some of the pd's variables that do not come from the data.txt file
* *
* Parameters *
* int i - for loop index
* *
* Returns *
* None (void) *
********************************************************/
PROGRAMDATA pdInit(PROGRAMDATA pd)
{
	pd.iter = 0;
	pd.r[0].y = 0.0;
	pd.r[0].v = 0.0;
	pd.r[0].a = 0.0;
	pd.r[0].t = 0.0;
	pd.r[0].Fd = 0.0;
	pd.r[0].Fd = 0.0;
	return pd;
}

/********************************************************
* num_simulation *
* *
* Purpose: Simulate numerical solution for Thermal Entry
* *
* Parameters *
* FILE *f - file pointer
* i,j,iter - loop counters for nodes horizontal, vertical,
*          and iteration number *
* W,H - horizontal and vertical nodes calculated
*       from input data *
* *
* Returns *
* None (void) *
********************************************************/
void num_simulation(PROGRAMDATA pd)
{
	int i;
	FILE *f;

	f = filewrite(pd);

	// error checking
	if (f == NULL)
	{
		printf("File can be read but not writen to.\n");
		getchar();
		exit(1);
	}

	fprintf(f, "iter,time,drag,position\n");
	fprintf(f, "0,0,0,0\n");

	pd.iter = 0;

	//This loop will look through every core temperature node
	/*The loop will continue as long as the Maximum residual of the
	current iteration is larger than the defined maximum residual*/
	do
	{
		pd.iter++;

		

	} while ((pd.r[pd.iter].t*pd.dt)<pd.tfinal && pd.iter <= MAX_ITER);

	fclose(f);
}

/********************************************************
* filewrite *
* *
* Purpose: Creates file names, checks for file existance, and streams files for writing.
* *
* Parameters *
* f - File stream, also the return
* resp - yes or no response
* filen - The filename to be composed procedurally
* j - static int acting as a switch. See below.
* ocase - copy of pd.scase, and comparator to the last time filewrite was called.
* *
* Returns *
* Initialized allocated array into PLATEPOINT structure nested in PROGRAMDATA *
********************************************************/
FILE *filewrite(PROGRAMDATA pd)
{

	FILE *f;    //Filestream pointer - the return
	char resp, filen[MAX_LINE_LENGTH]; //response char and filename composer

	//This static int will be used to first differentiate
	//between the RMS file and the plate simulation results file
	static int j = 0;

	//This static int is used to check to see if a new 
	//simulation has started since the last time the function was called.
	static int ocase = pd.scase;

	//If the simulation has changed, reset j to 0
	if (ocase != pd.scase) j = 0;
	ocase = pd.scase;//set ocase to scase, whether or not scase is differemt or not.

	//Write to first part of the file name.
	sprintf(filen, "Thermal Entry Simulation %d.csv", ocase);

	//attempts to open the file in read mode for checking
	f = fopen(filen, "r");


	if (f != NULL) //If the file exists...
	{

		printf("File \"%s\" exists. Ok to overwrite? (y/n): ", filen);

		fflush(stdin);
		scanf("%c", &resp);        // get response from user

		// if response is "no", tell user the file will not be overwritten
		if (resp == 'n' || resp == 'N')
		{
			printf("The existing file will not be overwritten.\n");
			fclose(f);     //close the file stream
			j++;
			return NULL;   //exit the function with no stream.
		}
	}//end overwrite check

	printf("Opening \"%s\" for writing\n\n", filen);

	//opens the file for writing
	f = fopen(filen, "w");

	j++;//increment j

	// if f returns a NULL send the user a message and return a 
	//NULL to the appropriate function that called filewrite 
	if (f == NULL)
	{
		printf("File can be read but not writen to!\n");
		getchar();
		return NULL;
	}

	return f;
}
/********************************************************
* nint *
* *
* Purpose: floors resulting double to nearest int
* *
* Returns *
* int for W and H calculations*
********************************************************/
int nint(double d)
{
	//Rounds to the nearest int
	return (int)floor(d + 0.5);
}

/********************************************************
* freepp *
* *
* Purpose: free allocated array
********************************************************/
void freepp(PROGRAMDATA pd)
{
	int i;
	free(pd.r);
	/*for (i = 0; i<4; i++){
		free(pd.tri_hold[i]);
	}*/
}