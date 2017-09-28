#include <iostream>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <float.h>

using namespace std;

/*Rosenbrock Function*/
double func(double* position, int number_of_dimensions) {
	int i;
	double result = 0;
	for(i = 0;i < number_of_dimensions-1;i++) {
		result += (100*pow(position[i+1] - position[i] * position[i],2) + pow(position[i] - 1,2));

	}
	return result;
}

double random_normal(double mean, double standard_deviation) {
	double U = rand()/(double)RAND_MAX;
	double V = rand()/(double)RAND_MAX;
	double normal_variable = sqrt(-2*log(U))*cos(2*M_PI*V);
	return (standard_deviation*normal_variable + mean);

}

double** allocate_memory_to_matrix(int row, int column) {
	int i,j;
	double** matrix = (double**) calloc(row,sizeof(double*));
	for(i =0 ; i < row;i++){
		
	 matrix[i] = (double*)calloc(column,sizeof(double));
	}
	return matrix;

}

void initialize_particle_position(double** particle_position, int swarm_size, int number_of_dimensions) {
	int i,j;
	for(i = 0 ; i < swarm_size;i++)
		for(j =0; j < number_of_dimensions;j++){

			particle_position[i][j] = random_normal(0.0,1.0);
		}

}
void initalize_best_known_position(double** best_known_position, double** particle_position, int swarm_size, int number_of_dimensions) {
	int i,j;
	for(i =0;i < swarm_size;i++) best_known_position[i] = (double*) malloc(sizeof(double)*number_of_dimensions);
	for(i =0;i <swarm_size;i++)
		for(j=0 ; j < number_of_dimensions;j++)
			best_known_position[i][j] = particle_position[i][j];
}

void calculate_score(double** particle_position, double* score, int swarm_size, int number_of_dimensions) {
	int i;
	for(i = 0; i < swarm_size;i++) score[i] = func(particle_position[i],number_of_dimensions);
}

double find_best_position(double* g, double* current_score,double** particle_position, int swarm_size, int number_of_dimensions) {
	int i,best_index;
	double best_score = DBL_MAX;
	for(i =0 ; i < swarm_size;i++) {
		if(current_score[i] < best_score){
			best_index = i;
			best_score = current_score[i];
		}
	}
	for(i = 0 ; i< number_of_dimensions;i++) g[i] = particle_position[best_index][i];
	return best_score;
}
void initialize_best_known_score(double* current_score, double* best_known_score, int swarm_size) {
	int i;
	for(i =0 ; i< swarm_size;i++) best_known_score[i] = current_score[i];
}

void initialize_velocity(double** velocity, int swarm_size, int number_of_dimensions) {
	int i,j;
	for(i = 0 ; i < swarm_size;i++)
		for(j =0; j < number_of_dimensions;j++)
			velocity[i][j] = random_normal(0.0,1.0);
}
void update_position(double** particle_position, double** velocity,double** best_known_position,double*g, double omega, double c1, double c2, int swarm_size, int number_of_dimensions){
	int i,j;
	for(i =0 ; i < swarm_size;i++)
		for(j=0;j < number_of_dimensions;j++)
		{
			velocity[i][j] = omega*velocity[i][j] + c1*random_normal(0.0,1.0)*(best_known_position[i][j] - particle_position[i][j]) + c2* random_normal(0.0,1.0) * (g[j] - particle_position[i][j]);
			particle_position[i][j] = particle_position[i][j] + velocity[i][j];
		}

}

void update_best_known_score(double** particle_position, double** best_known_position, double* current_score, double* best_known_score, int swarm_size, int number_of_dimensions)
{
	int i,j;
	for(i =0 ;i < swarm_size;i++){
		if(current_score[i] < best_known_score[i]) {
			best_known_score[i] = current_score[i];
			for(j =0;j<number_of_dimensions;j++) {
				best_known_position[i][j] = particle_position[i][j];
			}

		}
	}


}

int main() {

	//Seed random number generator
	struct timeval time;
    	gettimeofday(&time,NULL);
    	srand((time.tv_sec* 1000) + (time.tv_usec / 1000));

	int i,j;
	int swarm_size= 1000, number_of_dimensions = 5;
	int max_iteration = 5000;

	double** particle_position;
	double** velocity;
	double** best_known_position;
	double* current_score;
	double* best_known_score;
	double* g;

		

	double omega = 0.0, c1 =0.05, c2 =0.8;
	double best_score;

	
	particle_position = allocate_memory_to_matrix(swarm_size,number_of_dimensions);
	best_known_position = allocate_memory_to_matrix(swarm_size,number_of_dimensions);
	velocity = allocate_memory_to_matrix(swarm_size,number_of_dimensions);
	

	current_score = (double*) malloc(sizeof(double)*swarm_size);
	best_known_score = (double*) malloc(sizeof(double)*swarm_size);
	g = (double*) malloc(sizeof(double)*number_of_dimensions);

	
	initialize_particle_position(particle_position, swarm_size, number_of_dimensions);
	initalize_best_known_position(best_known_position,particle_position, swarm_size, number_of_dimensions);
	calculate_score(particle_position, current_score, swarm_size, number_of_dimensions);
	initialize_best_known_score(current_score, best_known_score, swarm_size);
	find_best_position(g,current_score,best_known_position, swarm_size, number_of_dimensions);
	
	initialize_velocity(velocity, swarm_size, number_of_dimensions);

	for(i =0; i < max_iteration;i++) {
		update_position(particle_position,velocity,best_known_position,g,omega,c1,c2,swarm_size,number_of_dimensions);
		calculate_score(particle_position, current_score, swarm_size, number_of_dimensions);
		update_best_known_score(particle_position, best_known_position, current_score, best_known_score, swarm_size, number_of_dimensions);
		best_score = find_best_position(g,best_known_score,best_known_position,swarm_size, number_of_dimensions);

		cout<<"ITERATION: "<<i+1<<" BEST SCORE: "<<best_score<<"\n";

	}

	cout<<"BEST PARTICLE\n";
	cout<<"[ ";
	for(i =0; i < number_of_dimensions;i++)cout<<g[i]<<" ";
	cout<<"]\n"; 



}

