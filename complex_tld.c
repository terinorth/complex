#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "myutil.h"

/*declare functions to be defined later*/
void readparams(void); // read in parameter file
void initialise(void); // set initial values
void evolve(void); // perform evolution between 0 and Stoptime
void mutate(int new); // create mutations, called from evolve()
void copy(int b1, int b2, int d); // mixes b1 and b2's DNA and allocates to d's space, called from evolve()
void clone(int x, int y); // called from copy()
void recomb(int b1, int b2, int d); // a recombination function which compares the mutations in two individuals, copies over the intersect and 
									// assigns probability of 0.5 for copying of any non-intersecting mutations
void purge(void); // a function to remove all fixed mutations from the population and adjust the baseline
					// trait and fitness value for each person. Called from evolve() and main(). 
					// also prints information to file about the frequency and effects of non-extinct mutations before purging. 
void printstuff(void); // output for post-simulation analysis
void free_mem(void); // a function to free all memory dynamically allocated at the end of the simulation 

/*mathematical functions from <math.h>*/
double nearbyint(double arg);

/*define a structure for containing the individual identifiers of the parents in each reproductive evolution. Used in evolve() and pickfit()*/
struct parents{
	int indiv1; // identifier of first individual - between 0 and Psize-1
	int indiv2; // identifier of second individual - between 0 and Psize-1 (can be the same as indiv1)
};

/*define pickfit() which needs struct parents to be defined first*/
struct parents* pickfit(void); // chooses mum and dad from population, called from evolve()

/*define the structure for containing information on mutational effects on fitness and the trait*/
struct gtype{ // linked list of mutational effects on (1)fitness and (2)trait
	double mutrait; // the mutational effects on the trait
	double mufit; // the mutational effects on fitness
	int age; // the 'age' of the mutation out of all mutations in the population - a number indicating the (relative) timing of this mutation
	struct gtype *fp; // the forward pointer 
	struct gtype *rp; // the reverse pointer
};

/*declare the add_mut() function which requires struct gtype to be defined first*/
void add_mut(struct gtype *n, struct gtype **last_mut, struct gtype **start_mut); // called from mutate()

/*individual specific variables*/
struct phen{
	double bday;  // birthday
	double trait; // trait of individual 
	double fit;   // fitness of individual
	struct gtype *last; // the latest ('youngest') mutation currently in the linked list for individual 'id' - this always has a forward pointer to NULL
						// (as the last mutation) and a reverse pointer to the next most recent mutation
	struct gtype *start; // the  'oldest' mutation currently in the linked list for individual 'id' - this always has a forward pointer to the next oldest
						// mutation and a reverse pointer to NULL
	int nmut; // number of mutations for this individual
};


/*declare the recalc_stats(), recomb1() and map() functions which require struct phen to be defined first*/
struct gtype* recomb1(struct phen ti, struct gtype *pi, struct gtype *curr, struct gtype *new); // called from recomb() -> copies the mutational information e.g. trait/fitness effects over to temporary individual
void map(struct phen new, int old); // called from recomb() -> maps the recombined information from temporary individual over to the new birth
struct phen recalc_stats(struct phen new); // called from recomb() -> recalculates statistics for new birth e.g. number of mutations


/*define a structure to contain the site frequency info and trait/fitness effects info for each mutation*/
struct mut_stats{
	int freq; // this is the total number of individuals in the population with this mutation (between 0 and Psize)
	double mut_eff_fit; // the effect of the mutation on fitness
	double mut_eff_trait; // the effect of the mutation on the trait
};


/*declare the ms() function, remove_this_mut() function and the equil() function which need struct mut_stats to be defined first*/
struct mut_stats* ms(void); // this is called from purge() - it is the specific function which assesses the frequency and trait/fitness effects
							// of each mutation in the population
void remove_this_mut(int age_mut, double mut_fiteff, double mut_traiteff); // called from purge() after ms() assessment - removes fixed mutations from each individual's linked list and updates
																			// baseline trait and fitness
void equil(struct mut_stats *rmi); // called from purge() before remove_this_mut(). This prints files to externally plot the current site frequency distribution and fitness/trait effects


/*define a structure to contain the reduced site frequency info and trait/fitness effects info for the population - only concerns mutations with non-zero/non-fixed frequency
 * This has an additional element which logs the age of each mutation, as the index of the resulting malloc'd array no longer relates to this
 */
struct mut_stats_r{
	int mutage; // the age of the mutation 
	int freq; // this is the total number of individuals in the population with this mutation (between 1 and Psize-1)
	double mut_eff_fit; // the effect of the mutation on fitness
	double mut_eff_trait; // the effect of the mutation on the trait
};

/*declare the reduce_muts() function which requires mut_stats_r to be defined first*/
struct mut_stats_r* reduce_muts(struct mut_stats *frmi); //this is called from printstuff() at the end of the simulaiton. It reduces the mutational info array down to only mutations with
														// non-zero/non-fixed frequency in the population. This is to make the trait_effects file as compact as possible



/*global variables*/
struct phen *Indiv; // declare a pointer, Indiv, which points to a structure, phen
int Psize; // define a variable for the pop size - this stays constant throughout the simulation (every death is 
			// matched with a birth). Psize is a parameter in the simulation
double Murate; // mutation rate - a parameter in the simulation
double Stoptime, Gentime; // time length of simulation, time length of generation; user-defined parameters in the simulation. 
							// Stoptime should always be a multiple of Gentime*Getstats_interval (see below) to ensure that
							// data is outputted correctly at the end of the simulation
double Betashape; // the shape paramater for equation (2) of Eyre-Walker (2010); parameter in the simulation
double MeanS; // the mean strength of selection in equation (2) of Eyre-Walker (2010); parameter in the simulation
int Ne; // the effective population size: this should always equal Psize and is given in the param file
double Tau; // parameter in equation (1) of Eyre-Walker - determines the strength of the relationship between the trait and fitness
int Mutcounter = 0; // Number of mutations in the entire population - update every time mutate() actions a new mutation
double Tt = 0.0; // Universal clock
double Purge_check = 0.0; // a record to check when it is time to purge all fixed mutations from the population
double Purge_interval; // this determines how often we implement the fixed mutation purging (number of generations when gentime=1) - defined in parameter file
double Basetrait = 100.0; // baseline trait value
double Basefit = 1.0;  // baseline fitness value
int Filename_count = 0; // this appends the filenames outputted by equil()
int Total_active_muts = 0; // this is a counter to log the number of mutations with non-zero frequency in the population - updated in equil() and used in get_stats()
int Trait_is_fit_switch = 0; // this is a switch (0 or 1) read in from the parameter file, to flag when trait=fitness (results in delta=1 and epsilon=0 in mutate())
double Ep_mu; // defined in param file - this is the mean of the epsilon normal rv in eqn(1) of Eyre-Walker(2010) - only used if Trait_is_fit_switch=0
double Ep_std; // defined in param file - this is the sd of the epsilon normal rv in eqn(1) of Eyre-Walker(2010) - only used if Trait_is_fit_switch=0
gsl_rng *R; // random number generator from gnu gsl --http://www.gnu.org/software/gsl/-- 

/*master function*/
int main(void) // 
{
	const gsl_rng_type *t;
	gsl_rng_env_setup();
	t=gsl_rng_default;
	R=gsl_rng_alloc(t);

	readparams(); // see below
	initialise(); // see below
	evolve(); // see below
	
	printf("%d\n\n\n",Mutcounter); // total number of mutations at the end of the simulation
	printf("%f\n\n\n",Tt); // check that this is close to Stoptime
		
	Filename_count=12345; // We change the Filename_count so we know that "12345" always refers to the end of the simulation
	purge();
	printstuff(); // see below 
	
	printf("%.14e\n\n\n",Basetrait); // the baseline trait at the end of the simulation
	FILE *basetrait_file;
	basetrait_file = fopen("basetrait_end.txt","w");
	fprintf(basetrait_file,"%.14e",Basetrait); // the baseline trait at the end of the simulation
	fclose(basetrait_file);
	
	/*free all memory dynamically allocated*/
	gsl_rng_free(R);
	free_mem();
}

/* NOTE:
 * We implement purge(), which calls on equil(), one final time after evolution is complete.
 * We label the data written during this phase as "12345"
 * It should mostly be the case that the plots from this data are exactly the same as the last plot created by equil() in purge() during evolution:
 * (Stoptime<=Purge_check<=Tt)
 * There may be a rare occasion when purge() doesn't initiate at the last iteration due to the finite representation of doubles:
 * (Stoptime<=Tt<=Purge_check) -->Purge_check has drifted higher than Stoptime for a Tt to be generated in between
 * We would always know when this happens because the number of plots created by purge()[equil()] will be less than Stoptime/Purge_check
 * It may also be the case that purge() initiates before the final value of Tt is generated in the simulation:
 * Purge_check<=Tt<=Stoptime...leading to...Purge_check<=Stoptime<=Tt --> Purge_check has drifted lower than Stoptime so Tt has been generated again 
 * This would lead to the last plot created by Purge_check being out of date at the end of the simulation
 * To ensure we always use the correct plots and data (trait_effects.txt) to describe the end of the simulation, we re-run purge() after evolution is complete
 * and implement ms() again in printstuff to account for the fact that a final purging of fixed alleles has taken place.
 * When checking equilbirium convergence, always substitue the "12345" file in the place of the last equil() data file, or include the "12345" file
 * if the last equil() data file is missing. This can be done manually when reviewing the plots created by complex_equil.R. complex_equil_ii.R
 * has been coded to perform this substitution. 
 * WARNING: in complex_equil_ii.r we rely on the fact that stoptime will always be exactly divisible by purge_interval. 
 * Every time we purge the fixed mutations, we also assess the sfd and other statistics within the current population using ms()/equil()
 */
 
 
 
/*define the called functions*/

void readparams(void)
{
	FILE *pfile;
	pfile = fopen("complex_params_tld.txt","r");
	fscanf(pfile,"%d %lf %lf %lf %lf %lf %d %lf %d %lf %lf %lf",&Psize,&Stoptime,&Gentime,&Murate, &Betashape, &MeanS, &Ne, &Tau, &Trait_is_fit_switch, &Ep_mu, &Ep_std, &Purge_interval); 
	fclose(pfile);
}

void initialise(void)
{
	int j;
	Indiv = (struct phen *)malloc(Psize*sizeof(struct phen));
	for(j=0;j<Psize;++j){ // increment across every individual in the population
		Indiv[j].bday = Tt; // every individual has their own 'clock'. We begin by setting everyone's birthday to Tt=0.0
		Indiv[j].trait = Basetrait; // some arbitrary value
		Indiv[j].fit = Basefit; // some arbitrary value. At this stage, fitness just manifests itself through differential fertility
		Indiv[j].nmut = 0; // no mutations initially
		Indiv[j].last=NULL; // no mutations to begin with, hence no youngest mutation
		Indiv[j].start=NULL; // no mutations to begin with, hence no oldest mutation
	}
}


void evolve(void)
{
	double tick, product1;
	struct parents *iborns=NULL;
	int iborn1,iborn2,idie,imut,j;
	
	tick = Gentime/(Psize/2.0);// [MB comment] why this? Gentime is standard Wright-Fisher generation time (e.g. in years)
	                               // In the Moran model the same amount of inbreeding is experienced as one Wright-Fisher
	                               // generation after N/2 birth-death events. So the waiting time between each birth-death
	                               // event should have expectation equal to Gentime/(Psize/2.0). 
	while(Tt < Stoptime){ // If we set Gentime to one, Stoptime will be the number of generations.
		Tt += gsl_ran_exponential(R,tick); // an exponential rv with expectation tick
		iborns=pickfit(); // choose two individuals to give birth based on fitness, could be the same individual
		iborn1 = iborns->indiv1; 
		iborn2 = iborns->indiv2; 
		idie = gsl_rng_uniform_int(R,Psize); // chosen uniformly randomly, and can include an individual that gives birth
		copy(iborn1,iborn2,idie); // copies with recombination 
		imut = gsl_ran_poisson(R,Murate); // A poisson rv with mean = Murate. We mutate at birth in line with convention.
		for(j=0;j<imut;++j)mutate(idie); // 08/11/12 modified to Poisson. We do this to follow the assumption of the Moran model.
										// with Murate=0.1, we expect exp(-0.1)=0.905 of the calls to give zero for imut
										//  - so then the loop is never entered because 0<0 is false
		if (Tt >= (Purge_check + (Purge_interval*Gentime))){
			purge(); 
			product1=Purge_interval*Gentime;
			Purge_check=Purge_check + nearbyint(product1); // update this for the next Purge() to take place in Purge_interval generations time
			
		}							
		free(iborns); // malloc'd in pickfit - need to free
		iborns=NULL; // to avoid spurious assignments after freeing 
	}
	

}

/*define the functions called by evolve */

struct parents* pickfit(void)
/*Pick two individuals to give birth weighted by fitness*/

{
	double r1, r2; // random numbers to determine which individuals are chosen (weighted by fitness) 
	double f1=0; // parameter for gsl_ran_flat
	double f2=1; // parameter for gsl_ran_flat
	struct parents *ret_parents; // a structure to contain the identifiers of the parents - returned to evolve()
	int par1, par2; // for recursion: the parents
	
	ret_parents = (struct parents *)malloc(sizeof(struct parents)); // allocate some memory to keep the parent identifiers
	
	while(1){
		r1=gsl_ran_flat(R,f1,f2); // a number of type double between 0 and 1 (maximum fitness) - actual range outputted by function is [0,1)
		par1=gsl_rng_uniform_int(R,Psize); // chooses an individual at random between 0 and (total number individuals-1)
		if(Indiv[par1].fit>r1){
			ret_parents->indiv1=par1;
			break;
		}
	}

	while(1){
		r2=gsl_ran_flat(R,f1,f2);
		par2=gsl_rng_uniform_int(R,Psize);
		if(Indiv[par2].fit>r2){
			ret_parents->indiv2=par2;
			break;
		}
	}
	
	return(ret_parents);
	
}


void mutate(int id)
/*generate mutational effects on the trait and fitness, which are linked by equation (1) in Eyre-Walker(2010)*/
/*id is the individual to mutate (the individual who dies and is replaced)*/

{
	double scale,S,epsilon;
	int delta;
	struct gtype *newg;
	double new_mutrait; // the new mutational effect on the trait
	double new_mufit; // the new mutational effect on fitness
	int new_age, new_S; // the 'age' of the new mutation out of all mutations in the population and a special case of "S" (see below)

	newg = (struct gtype *)malloc(sizeof(struct gtype));
	
	scale=MeanS/Betashape;
	S=gsl_ran_gamma(R,Betashape,scale); // generate a random gamma variable as the strength of selection; where S=4Ns, s=effect of mutation on fitness, N=effective 
										// population size - c/f Eyre-Walker (2010)  [note we use 2Ns as we have a haploid model]

	if (Trait_is_fit_switch==1){
		delta=1;
		epsilon=0;
		Tau=1; // just in case this isn't specified in the parameter file
	}
	else if (Trait_is_fit_switch==0){
		delta=gsl_rng_uniform_int(R,2)*2-1; // generate a rv with equal probability of being either -1 or 1
		epsilon=(gsl_ran_gaussian(R,Ep_std))+Ep_mu; // generate a normal rv with mean Ep_mu and standard deviation Ep_std
	}
	else printerr("mutate:trait_is_fit_switch problem");

	++Indiv[id].nmut; // prospectively update the mutation count for this individual
	++Mutcounter; // propectively update the population mutation count - this will be the 'age' of the new mutation
	
	new_mufit = -S/(2*Ne); //negative as we assume all mutations are deleterious. 2Ne not 4Ne as our model is haploid
	if (new_mufit<-1){
		new_mufit=-1; // restrict fitness effects as minimum -1 (any individual with this mutation will have fitness=0)	
		new_S=2*Ne;
		new_mutrait = delta*(pow(new_S,Tau))*(1+epsilon); // because s has been restricted to -1, this means that S=2Ne and we pass this on to the trait effect calculation
	}
	else {
	new_mutrait = delta*(pow(S,Tau))*(1+epsilon); // equation (1) of Eyre-Walker (2010)
	}
	
	new_age=Mutcounter; // the age of this new mutation is the current value of Mutcounter
	
	/*assign the appropriate values to newg*/
	newg->mutrait = new_mutrait;
	newg->mufit = new_mufit;
	newg->age = new_age;
	
	
	/*
	 * change the pointers so that newg is the last entry in the linked list - this happens because any mutation happening 'now' will always be 
	 * the youngest mutation for individual id by definition
	 */
	
	add_mut(newg, &Indiv[id].last, &Indiv[id].start);
	
	/*now update the individual's trait value and fitness value subject to additional mutations*/
	
	Indiv[id].trait = Indiv[id].trait + new_mutrait;
	
	Indiv[id].fit= Indiv[id].fit * (1+new_mufit);
		
}

void add_mut(struct gtype *n, struct gtype **last_mut, struct gtype **start_mut)
{
	if(!*last_mut){
		if (!*start_mut){
			*last_mut = n; // this is the first mutation to hit this individual, so redefine it as the 'youngest' mutation
			*start_mut = n; // similarly, define as the oldest mutation when this is the first mutation (this is needed for the docopy function)
			n->fp = NULL;
			n->rp = NULL;
		}
		else{
			printerr("add_mut: shouldn't be here");
		}
	}

	else {
		(*last_mut)->fp = n;
		n->fp = NULL;
		n->rp = *last_mut;
		*last_mut=n; // update the last mutation to be the current mutation - all previous mutations will be linked via reverse pointers
	}
		
}

void copy(int b1, int b2, int d)
{
	if(b1 == b2 && b1 == d){ // nothing happens (parents are the same, one parent dies - clonal)
		Indiv[b1].bday=Tt; // re-set the birthday to the current value of the universal clock
	}
	else if(b1 != b2){ // commonest case (parents are different)
		if(b1 != d && b2 != d){ // commonest case (neither parent dies)
			recomb(b1,b2,d); // normal recombination replacing the 'space' of individual d
		}
		else if(b1 == d)recomb(b1,b2,b1); // parent 1 dies - normal recombination, overwrite b1
		else if(b2 == d)recomb(b2,b1,b2); // parent 2 dies - normal recombination, overwrite b2
		else printerr("copy-recomb: shouldn't be here 1");
	}
	else if(b1 == b2)clone(d,b1); // parents are the same (hemaphrodite), another indiv dies - copy all of b1 into d
	else printerr("copy-clone: shouldn't be here 2");

}


void clone(int x, int y){ // clone all y mutations over to x

struct gtype *copyg, *curr_g, *next_g, *copyfp1, *copyfp2; 

copyg=NULL;
curr_g=NULL;
next_g=NULL;
copyfp1=NULL;
copyfp2=NULL;

Indiv[x].bday = Tt; // re-set the birthday of the individual
Indiv[x].trait=Indiv[y].trait; // same mutations so same trait value
Indiv[x].fit=Indiv[y].fit; // same mutations so same fitness value
Indiv[x].nmut=Indiv[y].nmut; // same number of mutations
		

/*first free the memory for x's previous mutations*/

if (Indiv[x].start!=NULL){
	copyfp1=Indiv[x].start;
}

while (copyfp1!=NULL){
	copyfp2=copyfp1->fp;
	free(copyfp1);
	copyfp1=copyfp2;
}

Indiv[x].start=NULL;
Indiv[x].last=NULL;

/*now make copies of y's mutations for x.
 *the reason we have to make copies rather than just link up is because each
 * mutation has a trait effect, a fitness effect, an age and importantly, a forward and
 * reverse pointer. While the trait/fitness effects and age of the mutation never change, the
 * forward and reverse pointers are individual specific and reflect how recombination
 * has shaped the mutational profile of a particular individual*/


/*intialise*/

if (Indiv[y].start!=NULL){
	copyg = (struct gtype *)malloc(sizeof(struct gtype));
	copyg->mutrait=Indiv[y].start->mutrait;
	copyg->mufit=Indiv[y].start->mufit;
	copyg->age=Indiv[y].start->age;
	copyg->rp=NULL;
	copyg->fp=NULL;
	curr_g=copyg; // curr_g now points to the same target as copyg
	Indiv[x].start=copyg; // Indiv[x].start now points to the same target as copyg
	next_g=Indiv[y].start->fp;
}

/*make copies of all mutations in the list*/

while (next_g!=NULL){
	copyg = (struct gtype *)malloc(sizeof(struct gtype)); // copyg now points to a new piece of memory, newly allocated
	copyg->mutrait=next_g->mutrait;
	copyg->mufit=next_g->mufit;
	copyg->age=next_g->age;
	copyg->rp=curr_g;
	copyg->fp=NULL;
	curr_g->fp=copyg;
	curr_g=copyg;
	next_g=next_g->fp;
}

Indiv[x].last=curr_g;
		
}


void recomb(int b1, int b2, int d){
/*
 * follow through the linked list of mutations for individual b1 and b2, copying over all mutations which occur in both individuals, and assigning p=0.5
 * for copying over of any remaining mutations in one individual only.
 * The mutations are tagged with (and ordered by) age of first occurence in the whole population, so if a mutation is in both individuals, this will be recognized by the 
 * unique age of the mutation. 
 */

/*begin at the start pointers for the oldest mutations in b1 and b2*/

struct phen temp_birth; // define a temporary individual - all info will be mapped into d at end
struct gtype *new;

new=NULL; // to avoid spurious assignments


temp_birth.bday = Tt; // the birthday is the current value of Tt
temp_birth.trait = Basetrait;
temp_birth.fit = Basefit; 
temp_birth.nmut = 0; // no mutations initially
temp_birth.last=NULL; // no mutations to begin with, hence no youngest mutation
temp_birth.start=NULL; // no mutations to begin with, hence no oldest mutation


struct gtype *b1_current, *b2_current, *d_current, *d_new; // these are the values of the pointers to mutations we are currently considering for persons b1 and b2, the 
															// latest mutation added to d, and the newly created d which will be added next

b1_current=Indiv[b1].start; // begin at the start of the linked lists
b2_current=Indiv[b2].start; 
d_current=NULL;
d_new=NULL;

while (b1_current!=NULL && b2_current!=NULL){ // check both b1 and b2 have mutations left to consider

if ((b1_current->age)==(b2_current->age)){ // if the ages of the mutations match ...
		if (temp_birth.start==NULL){
			temp_birth.start=recomb1(temp_birth, b1_current, d_current, d_new);
			temp_birth.start->fp=NULL; // define the forward and reverse pointers for d's mutation - currently these will point to nothing as this is the only mutation
			temp_birth.start->rp=NULL; 
			d_current=temp_birth.start; // d_current now points to d's latest (and only in this case) mutation
		}
	
		else if (temp_birth.start!=NULL){
			new=recomb1(temp_birth, b1_current, d_current, d_new); // do the same for the last but one mutation - link up with the newest mutation
			d_current->fp=new;
			new->rp=d_current; // we actually do this in the recomb1 function as well (remove?)
			d_current=new; // re-define the latest mutation added to d
				
		
		}
		b1_current=b1_current->fp; // move to the next mutation in b1 via b1's forward pointer saved in its current mutation
		b2_current=b2_current->fp; // repeat for b2

}


else if ((b1_current->age)>(b2_current->age)){ // b2's mutation is older (smaller 'age' value), hence we copy over b2 mutation with probability=0.5
	if (gsl_rng_uniform_int(R,2)==0){
		if (temp_birth.start==NULL){
			temp_birth.start=recomb1(temp_birth, b2_current, d_current, d_new);
			temp_birth.start->fp=NULL; // define the forward and reverse pointers for d's mutation - currently these will point to nothing as this is the only mutation
			temp_birth.start->rp=NULL; 
			d_current=temp_birth.start; // d_current now points to d's latest (and only in this case) mutation
		}
	
		else if (temp_birth.start!=NULL){
			new=recomb1(temp_birth, b2_current, d_current, d_new); // do the same for the last but one mutation - link up with the newest mutation
			d_current->fp=new;
			new->rp=d_current;
			d_current=new; // re-define the latest mutation added to d
				
		
		}
	}
	b2_current=b2_current->fp; // move to the next mutation for b2; irrespective of whether the current mutation was carried over to d
	
}


else if ((b1_current->age)<(b2_current->age)){
	if (gsl_rng_uniform_int(R,2)==0){
		if (temp_birth.start==NULL){
			temp_birth.start=recomb1(temp_birth, b1_current, d_current, d_new);
			temp_birth.start->fp=NULL; // define the forward and reverse pointers for d's mutation - currently these will point to nothing as this is the only mutation
			temp_birth.start->rp=NULL; 
			d_current=temp_birth.start; // d_current now points to d's latest (and only in this case) mutation
		}
	
		else if (temp_birth.start!=NULL){
			new=recomb1(temp_birth, b1_current, d_current, d_new); // do the same for the last but one mutation - link up with the newest mutation
			d_current->fp=new;
			new->rp=d_current;
			d_current=new; // re-define the latest mutation added to d
				
		
		}
		
	}
	b1_current=b1_current->fp; // move to the next mutation for b1; 
}
}



if (b1_current==NULL && b2_current!=NULL){ // no more mutations left for b1

	while (b2_current!=NULL){
		if (gsl_rng_uniform_int(R,2)==0){
			if (temp_birth.start==NULL){
			temp_birth.start=recomb1(temp_birth, b2_current, d_current, d_new);
			temp_birth.start->fp=NULL; // define the forward and reverse pointers for d's mutation - currently these will point to nothing as this is the only mutation
			temp_birth.start->rp=NULL; 
			d_current=temp_birth.start; // d_current now points to d's latest (and only in this case) mutation
			}
	
			else if (temp_birth.start!=NULL){
			new=recomb1(temp_birth, b2_current, d_current, d_new); // do the same for the last but one mutation - link up with the newest mutation
			d_current->fp=new;
			new->rp=d_current;
			d_current=new; // re-define the latest mutation added to d
				
		
			}

		}
		b2_current=b2_current->fp; // move to the next mutation for b2;
		
			
	}
}
	

if (b2_current==NULL && b1_current!=NULL){ // no more mutations left for b2
	while (b1_current!=NULL){
		if (gsl_rng_uniform_int(R,2)==0){
			if (temp_birth.start==NULL){
			temp_birth.start=recomb1(temp_birth, b1_current, d_current, d_new);
			temp_birth.start->fp=NULL; // define the forward and reverse pointers for d's mutation - currently these will point to nothing as this is the only mutation
			temp_birth.start->rp=NULL; 
			d_current=temp_birth.start; // d_current now points to d's latest (and only in this case) mutation
			}
	
			else if (temp_birth.start!=NULL){
			new=recomb1(temp_birth, b1_current, d_current, d_new); // do the same for the last but one mutation - link up with the newest mutation
			d_current->fp=new;
			new->rp=d_current;
			d_current=new; // re-define the latest mutation added to d
				
		
			}
		}
		b1_current=b1_current->fp; // move to the next mutation for b1
	}
}

	
temp_birth.last=d_current; // define the last mutation for d
temp_birth=recalc_stats(temp_birth); // re-calculate the fitness, trait value and number of mutations for the new individual

map(temp_birth,d); // a function which substitutes all values in d for corresponding values in temp_birth


}


struct phen recalc_stats(struct phen new){
/*function to parse through the linked list of mutations and re-calculate the new fitness and trait values for d, in addition to the 
 * total number of mutations for d
 */

double newtrait, newfit;
int new_nmut;
struct gtype *curr_mut;

curr_mut=new.start;
newtrait=new.trait; // baseline
newfit=new.fit; // baseline
new_nmut=new.nmut; // baseline


while (curr_mut!=NULL){
	++new_nmut;
	newtrait = newtrait+curr_mut->mutrait;
		
	newfit = newfit * (1+curr_mut->mufit);
	
	curr_mut=curr_mut->fp; // move on to the next mutation
}


new.trait=newtrait; // update the summary trait/fitness/number mutations statistics for the birth
new.fit=newfit;
new.nmut=new_nmut;
return(new);

}


struct gtype* recomb1(struct phen ti, struct gtype *pi, struct gtype *curr, struct gtype *new){
/*this function copies mutations over from parent (pi) to temporary new individual (ti)*/
		if (ti.start==NULL){ // if this is the first mutation for d ... 
			new = (struct gtype *)malloc(sizeof(struct gtype));
			new->mutrait=pi->mutrait; 
			new->mufit=pi->mufit;
			new->age=pi->age;
			return(new);
		}
		else if (curr!=NULL){ // if this isn't the first mutation for d ... 
			new = (struct gtype *)malloc(sizeof(struct gtype));
			new->mutrait=pi->mutrait;
			new->mufit=pi->mufit;
			new->age=pi->age;
			new->fp=NULL; // re-define the forward pointer to be NULL, as this is d's latest mutation
			new->rp=curr; // re-define the reverse pointer to link back with d's other mutations 
			return(new);
	
	
		}

}

void map(struct phen new, int old){

/*free the memory for the old mutations*/

struct gtype *copyfp1, *copyfp2;
copyfp1=NULL;
copyfp2=NULL;

if (Indiv[old].start!=NULL){
	
		copyfp1=Indiv[old].start;
}

while (copyfp1!=NULL){
	copyfp2=copyfp1->fp;
	free(copyfp1);
	copyfp1=copyfp2;
}

Indiv[old].start=NULL;
Indiv[old].last=NULL;

/*now link in the new mutations for the individual who died to signify the birth*/
	
Indiv[old].bday = new.bday;
Indiv[old].trait = new.trait;
Indiv[old].fit = new.fit;
Indiv[old].nmut = new.nmut;
Indiv[old].last= new.last;
Indiv[old].start= new.start;
	
	
	
}

void purge(void){
    /*
     *  
     * When everyone has the mutation (number=Psize), delete this from everyone's 
     * linked list, freeing the memory in the process and updating the baseline fitness and trait values across the whole
	 * population. 
	 * We do this to speed up the script as otherwise too much time is spent copying fixed mutations, which are uninformative
	 * 
	 */

struct mut_stats *returned_mut_stats_purge; // a structure to contain the frequency/trait effects/fitness effects 
											// of each mutation that has ever been generated in the population
int i; // for recursion
int mutation_age; // calculated using index of returned_mut_stats_purge

returned_mut_stats_purge=ms(); // parse through the linked list of mutations for each individual and calculate the 
								// mutation frequencies in the population, additionally recording the trait and fitness effects
								// note: if previous calls to purge() have taken place then fixed mutations will appear to be extinct with freq=0

equil(returned_mut_stats_purge); // print to file the frequencies and effects of all non-extinct and non-fixed mutations. 
								// this is to monitor the evolution of the system and assess equilibrium convergence
								
for (i=0;i<Mutcounter;++i){
	if (returned_mut_stats_purge[i].freq==Psize){
		mutation_age=i+1;
		remove_this_mut(mutation_age, returned_mut_stats_purge[i].mut_eff_fit, returned_mut_stats_purge[i].mut_eff_trait); 
		// purge the population of all fixed mutations
	}
}

free(returned_mut_stats_purge); // free memory which was malloc'd in ms()

}


struct mut_stats* ms(void){
	
int i, j, l; 
struct gtype *this_mut; // a temporary variable to parse through the linked list of mutations for each individual
struct mut_stats *mut_info; // a pointer to where the info on the site freq spectrum/mutational effects info is stored
mut_info = (struct mut_stats *)malloc(Mutcounter*sizeof(struct mut_stats));
for (l=0; l<Mutcounter; ++l){
	mut_info[l].freq=0; // initial values
	mut_info[l].mut_eff_fit=0;
	mut_info[l].mut_eff_trait=0;
}
this_mut=NULL; // to avoid spurious assignments	
	 
for (i=0; i<Psize;++i){
	j=0;
				
	if (Indiv[i].start!=NULL){ // if individual i has mutations
		this_mut=Indiv[i].start;
		while (this_mut){
			while ((this_mut->age)>(j+1)){
				++j;
			}
			if ((j+1)==(this_mut->age)){
				++mut_info[j].freq;
				mut_info[j].mut_eff_fit=this_mut->mufit;
				mut_info[j].mut_eff_trait=this_mut->mutrait;
			}
			this_mut=this_mut->fp;
					
		}
		
	}
				
}
	
return(mut_info);


}

void remove_this_mut(int age_mut, double mut_fiteff, double mut_traiteff){
/* Seek out the identified fixed mutation (argument to function) from everyone's linked list and delete it.
 * Update the baseline trait and fitness effects so that any new individuals have their trait and fitness adjusted accordingly
 * (they will all have these mutations but we only represent this in the baseline trait/fitness values).
 */

int i; // for recursion
int b; // a binary switch for QC purposes - check to ensure that we locate and remove the mutation for every individual
struct gtype *curr_mut=NULL; // allows us to parse through each individual's linked list of mutations 
struct gtype *mut_before=NULL; // for re-arranging pointers following deletion
struct gtype *mut_after=NULL; // for re-arranging pointers following deletion

/*update baseline values for future births*/
FILE *basestats_file;
basestats_file = fopen("base_updates.txt","a");
Basetrait=Basetrait+mut_traiteff;
fprintf(basestats_file,"%.14e\n",mut_traiteff);
fclose(basestats_file);
Basefit=Basefit * (1+mut_fiteff);

/*find the fixed mutations and remove them from each individual's linked list*/
for (i=0; i<Psize; ++i){ // for each individual in the population
	b=0; // re-set switch
	curr_mut=Indiv[i].start; // start at the oldest mutation (smallest mutation age)
	if (curr_mut!=NULL){ // this should always be the case as the mutation is fixed
		while(curr_mut){
		if (curr_mut->age==age_mut){
			if (curr_mut==Indiv[i].start){ // The oldest mutation for this individual is the target fixed mutation
				Indiv[i].start=curr_mut->fp; // Re-define the first mutation for this individual as we are deleting from the list
				if (Indiv[i].start!=NULL){ // This was not the only mutation for this individual
					Indiv[i].start->rp=NULL;
				}
				else if (Indiv[i].start==NULL){ // This was the only mutation for this individual
					Indiv[i].last=NULL; 
				}
				free(curr_mut);
				b=1; //QC switch for check below
				curr_mut=NULL; // to break out of the while loop
			}
			else if (curr_mut==Indiv[i].last){ // The youngest mutation for this individual is the target fixed mutation
				Indiv[i].last=curr_mut->rp; // Re-define the last mutation for this individual as we are deleting from the linked list
				Indiv[i].last->fp=NULL; // We know there must be another mutation available because the case of one mutation has been dealt with above
				free(curr_mut);
				b=1; //QC switch for check below
				curr_mut=NULL; // to break out of the while loop - we need to re-define the pointer because it is pointing to nonsense
			}
			else{
				mut_before=curr_mut->rp;
				mut_after=curr_mut->fp;
				mut_before->fp=mut_after;
				mut_after->rp=mut_before;
				free(curr_mut);
				b=1; //QC switch for check below
				curr_mut=NULL; // to break out of the while loop				
			}
		}
		else if ((curr_mut->age)<(age_mut)){
			curr_mut=curr_mut->fp; // check the next mutation
		}
		else if ((curr_mut->age)>(age_mut)){
			printerr("remove_this_mut: unable to locate fixed mutation1"); 
		}
		}
		if (b!=1){
			printerr("remove_this_mut:unable to locate fixed mutation2");
		}
	}
	else{
	printerr("remove_this_mut:individual has no mutations");
	}
}
}


void equil(struct mut_stats *rmi){
/*print to file the frequency and trait/fitness effects of all mutations with non-zero & non-fixed frequency*/

/*for simple plotting in R*/
	int i; // for recursion
	FILE *mutstats_file;
	char filename[50];
	snprintf(filename,50,"mut_stats_%d.txt", Filename_count );
	mutstats_file = fopen(filename,"w");
	Total_active_muts=0; //  reset to 0 before re-calculating 
	for (i=0; i<Mutcounter;++i){
		if (rmi[i].freq!=0 && rmi[i].freq!=Psize){
			++Total_active_muts;
			fprintf(mutstats_file,"%d\t%d\t%.14e\t%.14e\n",i,rmi[i].freq,rmi[i].mut_eff_fit,rmi[i].mut_eff_trait);
		}
	}
	fclose(mutstats_file);
	++Filename_count;

}


void printstuff(void)
{
/* Here we want to print out the final data to analyse. This is a table with each
 * individual on a different row, their fitness value, their trait value, and the effects of their mutations on the trait.
 * The columns for the mutations will be ordered by Mutcounter and if an individual doesn't 
 * have a mutation we set the effect to 0.
 * 
 * Note: fixed mutations (which will appear to be extinct) and mutations with genuine zero frequency are skipped and not printed (as in equil())
 * 
 * A new function, reduce_muts(), populates an array with all non-zero/non-fixed freq mutations and their respective frequencies,
 * fitness and trait effects. This is used compare against each individual's linked list of mutations.
 */


/*final mutation frequency distribution*/

struct mut_stats *final_returned_mut_info;
struct mut_stats_r *reduced_returned_mut_info;

final_returned_mut_info=ms(); 

reduced_returned_mut_info=reduce_muts(final_returned_mut_info); // equil() must always be run before reduced_returned_mut_info() as equil() updates 
																// the global variable Total_active_muts which reduced_returned_mut_info() requires
																// to malloc the correct amount of memory. Equil is called by purge() after evolution has
																// finished, and only catalogues non-fixed/non-extinct mutations so ok to proceed.

/*the individual - trait effects table*/

FILE *trait_file;
trait_file = fopen("trait_effects.txt","w");
		
int i; // for recursion
int j; // for recursion
int k; // for recursion
struct gtype *current_mut; // the current mutation under consideration

current_mut=NULL; // to avoid spurious assignments

for(i=0;i<Psize;++i){
	
fprintf(trait_file,"%d\t%.14e\t%.14e\t",i,Indiv[i].fit,Indiv[i].trait); // begin by printing the individual's index and their total fitness and trait value

if (Indiv[i].start==NULL){
	for (k=0; k<Total_active_muts-1;++k){
		fprintf(trait_file,"%d\t",0); // for individuals with no mutations, print out a 0 effect for each mutation with non-zero frequency that we are considering
	}
	fprintf(trait_file,"%d\n",0); // for the final mutation, include an end of line character
}

else if (Indiv[i].start!=NULL){
current_mut=Indiv[i].start;	

for(j=0;j<Total_active_muts;++j){

if (current_mut!=NULL && j!=(Total_active_muts-1)){
	if (current_mut->age==reduced_returned_mut_info[j].mutage){
		fprintf(trait_file,"%.14e\t", current_mut->mutrait); // print out the mutational effect
		current_mut=current_mut->fp;
	}

	else if (current_mut->age!=reduced_returned_mut_info[j].mutage){
		fprintf(trait_file,"%d\t", 0); // print out a zero mutational effect
	}
				
}

else if (current_mut!=NULL && j==(Total_active_muts-1)){
	if (current_mut->age==reduced_returned_mut_info[j].mutage){
		fprintf(trait_file,"%.14e\n", current_mut->mutrait); // print out the mutational effect
		current_mut=current_mut->fp;
	}

	else if (current_mut->age!=reduced_returned_mut_info[j].mutage){
		printerr("print stuff: shouldn't be here 3");
	}
				
} 	

else if (current_mut==NULL && j!=(Total_active_muts-1)){
	fprintf(trait_file,"%d\t", 0); // print out a zero mutational effect	
}

else if (current_mut==NULL && j==(Total_active_muts-1)){
	fprintf(trait_file,"%d\n", 0); // print out a zero mutational effect	
}


	
}

}
}	

fclose(trait_file);
free(final_returned_mut_info);
free(reduced_returned_mut_info); 

}





struct mut_stats_r* reduce_muts(struct mut_stats *frmi){
/*
 *reduce_muts() creates a structure containing the age, freq, effect on fitness and effect on the trait, of each
 *non-zero/non-fixed frequency mutation in the population. Note: all fixed mutations should be represented as zero frequency
 *as they would have been deleted from the linked lists by the last iteration of purge() after evolution was complete
 *
 */
	int i, j; // for recursion
	int count=0; // sanity check and for recursion
	struct mut_stats_r *reduced_mut_info; // a pointer to where the info on the site freq spectrum/mutational effects info is stored (non-zero freq muts only)
	reduced_mut_info = (struct mut_stats_r *)malloc(Total_active_muts*sizeof(struct mut_stats_r));
	for (j=0; j<(Total_active_muts-1); ++j){
		reduced_mut_info[j].mutage=-1; // dummy variable initially
		reduced_mut_info[j].freq=0; // initialise
		reduced_mut_info[j].mut_eff_fit=0; // initialise
		reduced_mut_info[j].mut_eff_trait=0;	// initialise
	}
	
	for (i=0; i<Mutcounter;++i){
		if (frmi[i].freq!=0 && frmi[i].freq!=Psize){ // we check against Psize but really this is redundant as all fixed muts should have freq=0
			reduced_mut_info[count].mutage=i+1; // the index of final_returned_mut_info should correspond to the age of the mutation-1
			reduced_mut_info[count].freq=frmi[i].freq;
			reduced_mut_info[count].mut_eff_fit=frmi[i].mut_eff_fit;
			reduced_mut_info[count].mut_eff_trait=frmi[i].mut_eff_trait;
			++count;	
		}
		
	}
	if (count!=Total_active_muts){ // count should start at 0 and populate reduced_mut_info to the last element (Total_active_muts -1). Count will be updated by 1
									// after the last non-zero mutation has been logged. This means that count should always be exactly equal to Total_active_muts
		printerr("problem with reduce_muts");

	}
	
	return(reduced_mut_info);
}



void free_mem(void){

int i; // for recursion
struct gtype *copyfp1, *copyfp2; // temporary pointers to free the memory for each mutation

copyfp1=NULL; // to avoid spurious assignments
copyfp2=NULL; 



for (i=0;i<Psize;++i){


if (Indiv[i].start!=NULL){
		copyfp1=Indiv[i].start;
}

while (copyfp1!=NULL){
	copyfp2=copyfp1->fp;
	free(copyfp1);
	copyfp1=copyfp2;
}

Indiv[i].start=NULL;
Indiv[i].last=NULL;

}

free(Indiv);

}
