#include<complex.h>
#include<ctype.h>
#include<math.h>
#include<pthread.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

// Function declarations
void compute_exact_roots();
void * compute_main(void * args);
void compute_item(int item);
void * write_main(void * args);

// Gloabal variables
int ** roots;
int ** iterations;
double complex * exact_roots;
char * item_done;
int nmb_threads;
int size;
int exponent;
pthread_mutex_t item_done_mutex;

void main(int argc, char* argv[]) {
  // Read arguments
  if (argc != 4 || strncmp(argv[1],"-t",2) != 0 || strncmp(argv[2],"-l",2) != 0 ||
      !isdigit(argv[3][0])) {
    printf("Wrong arguments!");
    exit(0);
  }
  nmb_threads = atoi(&argv[1][2]);
  size = atoi(&argv[2][2]);
  exponent = atoi(&argv[3][0]);
  
  // Initialize global variables
  compute_exact_roots();
  roots = (int**) malloc(sizeof(int*)*size);
  iterations = (int**) malloc(sizeof(int*)*size);
  item_done = (char*) malloc(sizeof(char)*size);

  // Start compute threads
  pthread_t * compute_threads = malloc(sizeof(pthread_t)*nmb_threads);
  int ret;
  for (size_t ix = 0; ix < nmb_threads; ++ix) {
    size_t * args = malloc(sizeof(size_t));
    *args = ix;
    if (ret = pthread_create(&compute_threads[ix], NULL, compute_main, (void*)args)) {
      printf("Error creating thread: %\n", ret);
      exit(1);
    }
  }
  
  pthread_t write_thread;
  if (ret =  pthread_create(&write_thread, NULL, write_main, NULL)) {
    printf("Error creating threaf: %\n",ret);
    exit(1);
  }

  pthread_join(write_thread,NULL);
}

void compute_exact_roots() {
  double angle;
  double complex * a = malloc(sizeof(double complex) * (exponent + 1));
  for (size_t ix = 0; ix < exponent; ++ix) {
    angle = (2 / (double) exponent) * ix  * M_PI;
    a[ix] = cos(angle) + sin(angle)*I;
  }
  a[exponent] = 0;
  exact_roots = a;
}

void * compute_main(void * args) {
  size_t offset = *((size_t*)args);
  free(args);
  
  for (size_t ix = offset; ix < size; ix += nmb_threads ) {
    compute_item(ix);
    printf("Klar med item %d\n",ix);
    pthread_mutex_lock(&item_done_mutex);
    item_done[ix] = 1;
    pthread_mutex_unlock(&item_done_mutex);
  }
}

void compute_item(int item) {
  int * roots_result = malloc(sizeof(int) * size);
  int * iterations_result = malloc(sizeof(int) * size);
  double step = 4 / (double) (size-1);
  double complex yI = (-2 + item*step)*I;
  for (size_t ix = 0; ix < size; ++ix) {
    double complex z = -2 + ix*step + yI;
    int done = 0;
    if (z == 0) {
      roots_result[ix] = exponent;
      iterations_result[ix] = 0;
      done = 1;
    }
    for (size_t jx = 0; done == 0; ++jx) {
      double complex z_d_1 = 1;
      for (size_t kx = 0; kx < exponent-1; ++kx)
	z_d_1 = z_d_1*z;
      z = z - (z*z_d_1-1)/exponent/z_d_1;
      for (size_t kx = 0; kx <= exponent; ++kx)
	if (cabs(z-exact_roots[kx]) < 1e-3) {
	  roots_result[ix] = kx;
	  iterations_result[ix] = jx;
	  done = 1;
	  break;
	}
      if (abs(creal(z)) > 1e10 || abs(cimag(z)) > 1e10) {
	roots_result[ix] = exponent;
	iterations_result[ix] = jx;
	done = 1;
      }
    }
  }
  roots[item] = roots_result;
  iterations[item] = iterations_result;
}

void * write_main(void * args) {
  char * item_done_loc = calloc(size, sizeof(char));
  struct timespec sleep_timespec = {.tv_sec = 0, .tv_nsec = 500L};
  char * intro_string = "P3 \n%d %d \n255 \n";
  char * template_string = "%3d %3d %3d ";
  char tmp_string[30];
  sprintf(tmp_string, "newton_convergence_x%d.ppm", exponent);
  FILE * f_roots = fopen(tmp_string,"w+");
  sprintf(tmp_string, "newton_attractors_x%d.ppm", exponent);
  FILE * f_iterations = fopen(tmp_string,"w+");
  sprintf(tmp_string,intro_string,size,size);
  fwrite(tmp_string,sizeof(char),strlen(tmp_string),f_roots);
  fwrite(tmp_string,sizeof(char),strlen(tmp_string),f_iterations);
  for (size_t ix = 0; ix < size;) {
    pthread_mutex_lock(&item_done_mutex);
    if (item_done[ix] != 0) {
      memcpy(item_done_loc, item_done, size*sizeof(char));
    }
    pthread_mutex_unlock(&item_done_mutex);

    if (item_done_loc[ix] == 0) {
      nanosleep(&sleep_timespec, NULL);
      continue;
    }

    for (; ix < size  && item_done_loc[ix] != 0; ++ix ) {
      int * roots_results = roots[ix];
      int * iterations_results = iterations[ix];
      
      for (size_t jx = 0; jx < size; ++jx) {
	int rcol = roots_results[jx]/(double)exponent*255;
	sprintf(tmp_string,template_string,rcol,rcol,rcol);
	fwrite(tmp_string,sizeof(char),strlen(tmp_string),f_roots);
	int icol = roots_results[jx]/(double)50*255;
	sprintf(tmp_string,template_string,icol,icol,icol);
       	fwrite(tmp_string,sizeof(char),strlen(tmp_string),f_iterations); 
      }
      fwrite("\n",sizeof(char),1,f_roots);
      fwrite("\n",sizeof(char),1,f_iterations);

      free(roots_results);
      free(iterations_results);
    }
  }
}
