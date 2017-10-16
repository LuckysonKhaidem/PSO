#include <Python.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <float.h>



double random_normal(double mean, double standard_deviation) {
	double U = rand()/(double)RAND_MAX;
	double V = rand()/(double)RAND_MAX;
	double normal_variable = sqrt(-2*log(U))*cos(2*M_PI*V);
	return (standard_deviation*normal_variable + mean);

}

static PyObject* convert_into_python_list(double* particle_position, int number_of_dimensions) {
	int i;
	PyObject* result = PyList_New(0);
	for(i = 0; i < number_of_dimensions;i++) {
		PyList_Append(result,Py_BuildValue("d", particle_position[i]));
	}
	return result;
}

double** allocate_memory_to_matrix(int row, int column) {
	int i;
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

void calculate_score(PyObject* objective_function,double** particle_position, double* score, int swarm_size, int number_of_dimensions) {
	int i;

	for(i = 0; i < swarm_size;i++){

		PyObject* particle_position_as_pylist = convert_into_python_list(particle_position[i],number_of_dimensions);
		PyObject* arg = PyTuple_New(1);
		
		PyTuple_SetItem(arg,0,particle_position_as_pylist);
		score[i] = PyFloat_AsDouble(PyEval_CallObject(objective_function, arg));
		Py_DECREF(particle_position_as_pylist);
		Py_DECREF(arg);
	}
}

double find_best_position(double* g, double* current_score,double** particle_position,double* worst_score, int swarm_size, int number_of_dimensions) {
	int i,best_index;
	double best_score = DBL_MAX;
	*worst_score = DBL_MIN;
	for(i =0 ; i < swarm_size;i++) {
		if(current_score[i] < best_score){
			best_index = i;
			best_score = current_score[i];
		}
		else if(current_score[i] > *worst_score) *worst_score = current_score[i];
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
void update_position(double** particle_position, double** velocity,double** best_known_position,double*g,double* lb, double* ub, double omega, double c1, double c2, int swarm_size, int number_of_dimensions){
	int i,j;
	for(i =0 ; i < swarm_size;i++)
		for(j=0;j < number_of_dimensions;j++)
		{
			velocity[i][j] = omega*velocity[i][j] + c1*random_normal(0.0,1.0)*(best_known_position[i][j] - particle_position[i][j]) + c2* random_normal(0.0,1.0) * (g[j] - particle_position[i][j]);
			particle_position[i][j] = particle_position[i][j] + velocity[i][j];

			/*Perform boundary check*/
			if(particle_position[i][j] < lb[j]) particle_position[i][j] = lb[j];
			else if (particle_position[i][j] > ub[j]) particle_position[i][j] = ub[j]; 
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

void convert_python_list_into_c_array(PyObject* pylist, double* c_array, int n) {
	int i;
	for(i = 0 ; i < n; i++){
		c_array[i] = PyFloat_AsDouble(PyList_GetItem(pylist,i));
	}
}

static PyObject* run_pso(PyObject* self, PyObject* args, PyObject* kwargs) {

	//Parameter names
	static char* parameter_names[] = {"objective_function","swarm_size","number_of_dimensions","omega","c1","c2","tolerance","lower_bound","upper_bound","max_iteration",NULL};

	//Seed random number generator
	struct timeval time;
    gettimeofday(&time,NULL);
    srand((time.tv_sec* 1000) + (time.tv_usec / 1000));

	int i,j;
	int swarm_size,number_of_dimensions;
	int max_iteration;

	double omega,c1,c2;

	double** particle_position;
	double** velocity;
	double** best_known_position;
	double* current_score;
	double* best_known_score;
	double* g;

	double tolerance;
	PyObject* objective_function;
	PyObject* upper_bound;
	PyObject* lower_bound;

	if(!PyArg_ParseTupleAndKeywords(args,kwargs,"OiiddddOOi",parameter_names,&objective_function,&swarm_size,&number_of_dimensions,&omega,&c1,&c2,&tolerance,&lower_bound,&upper_bound,&max_iteration)) {
		PyErr_SetString(PyExc_TypeError, "Invalid Parameters");
		return NULL;
	}
	if (!PyCallable_Check(objective_function)) {
            PyErr_SetString(PyExc_TypeError, "objective_function must be a callable");
            return NULL;
     }
    if(!PyList_Check(lower_bound) || !PyList_Check(upper_bound)) {
    	PyErr_SetString(PyExc_TypeError, "lower_bound and upper_bound should be a list");
        return NULL;

    }
    if(PyList_Size(lower_bound) != number_of_dimensions || PyList_Size(upper_bound) != number_of_dimensions) {
    	PyErr_SetString(PyExc_TypeError, "The size of upper_bound and lower_bound should be equal to the number of number_of_dimensions");
        return NULL;

    }
	double best_score,worst_score;
	double* ub = (double*)malloc(sizeof(double)*number_of_dimensions);
	double* lb = (double*)malloc(sizeof(double)* number_of_dimensions);

	convert_python_list_into_c_array(lower_bound,lb,number_of_dimensions);
	convert_python_list_into_c_array(upper_bound,ub,number_of_dimensions);
	
	particle_position = allocate_memory_to_matrix(swarm_size,number_of_dimensions);
	best_known_position = allocate_memory_to_matrix(swarm_size,number_of_dimensions);
	velocity = allocate_memory_to_matrix(swarm_size,number_of_dimensions);
	
	current_score = (double*) malloc(sizeof(double)*swarm_size);
	best_known_score = (double*) malloc(sizeof(double)*swarm_size);
	g = (double*) malloc(sizeof(double)*number_of_dimensions);

	
	printf("Initializing...\n");
	printf("Swarm Size: %d\n Number of Dimensions: %d\n Omega: %f\n C1: %f\n C2: %f\n",swarm_size,number_of_dimensions,omega,c1,c2);
	initialize_particle_position(particle_position, swarm_size, number_of_dimensions);

	initalize_best_known_position(best_known_position,particle_position, swarm_size, number_of_dimensions);
	calculate_score(objective_function,particle_position, current_score, swarm_size, number_of_dimensions);
	initialize_best_known_score(current_score, best_known_score, swarm_size);
	find_best_position(g,current_score,best_known_position,&worst_score, swarm_size, number_of_dimensions);
	initialize_velocity(velocity, swarm_size, number_of_dimensions);

	printf("Running PSO..\n");
	for(i =0; i < max_iteration;i++) {

		update_position(particle_position,velocity,best_known_position,g,lb,ub,omega,c1,c2,swarm_size,number_of_dimensions);
		calculate_score(objective_function,particle_position, current_score, swarm_size, number_of_dimensions);
		update_best_known_score(particle_position, best_known_position, current_score, best_known_score, swarm_size, number_of_dimensions);
		best_score = find_best_position(g,best_known_score,best_known_position,&worst_score,swarm_size, number_of_dimensions);

		// cout<<"ITERATION: "<<i+1<<" BEST SCORE: "<<best_score<<"\n";
		printf("ITERATION: %d BEST SCORE: %f\n",i+1,best_score);
		if((worst_score - best_score) <= tolerance) break;

	}

	PyObject* best_particle_position = convert_into_python_list(g, number_of_dimensions);
	free(current_score);
	free(best_known_score);
	free(g);
	return best_particle_position;

}

static PyMethodDef methods[] =  {
	{"optimize",(PyCFunction)run_pso,METH_KEYWORDS | METH_VARARGS,NULL},
	{NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initpso(void) {
	(void) Py_InitModule3("pso",methods,NULL);
}



