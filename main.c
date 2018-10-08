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
int ** attractors;
int ** convergences;
double complex * exact_roots;
char * item_done;
int nmb_threads;
int size;
int expo;
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
  expo = atoi(&argv[3][0]);
  
  // Initialize global variables
  compute_exact_roots();
  attractors = (int**) malloc(sizeof(int*)*size);
  convergences = (int**) malloc(sizeof(int*)*size);
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
    printf("Error creating thread: %\n",ret);
    exit(1);
  }

  pthread_join(write_thread,NULL);
  free(attractors);
  free(convergences);
  free(exact_roots);
  free(item_done);
}

void compute_exact_roots() {
  double angle = M_PI * 2.0 / expo;
  double complex * tmp_exact_roots = malloc(sizeof(double complex) * (expo + 1));
  for (size_t ix = 0; ix < expo; ++ix)
    tmp_exact_roots[ix] = cos(angle*ix)+sin(angle*ix)*I;
  tmp_exact_roots[expo] = 0;
  exact_roots = tmp_exact_roots;
}

void * compute_main(void * args) {
  size_t ix = *((size_t*)args);
  free(args);
  for (; ix < size; ix += nmb_threads) {
    compute_item(ix);
    pthread_mutex_lock(&item_done_mutex);
    item_done[ix] = 1;
    pthread_mutex_unlock(&item_done_mutex);
  }
}

void compute_item(int item) {
  int * attr_result = malloc(sizeof(int) * size);
  int * conv_result = malloc(sizeof(int) * size);
  double step = 4.0 / (size-1);
  double complex y_I = (-2 + item*step)*I;
  double expo_inv = 1.0/expo;
  double one_min_expo_inv = 1 - expo_inv;
  for (size_t ix = 0; ix < size; ++ix) {
    double complex z = -2 + ix*step + y_I;
    int done = 0;
    for (size_t jx = 0; done == 0; ++jx) {
      double abs_z_2 = creal(z)*creal(z) + cimag(z)*cimag(z);
      for (size_t kx = 0; kx < expo; ++kx) {
	//if (cabs(z-exact_roots[kx]) < 1e-3) {
	if (abs_z_2 + 1 -2*(creal(z)*creal(exact_roots[kx])+cimag(z)*cimag(exact_roots[kx])) < 1e-6) {
	  attr_result[ix] = kx;
	  conv_result[ix] = jx;
	  done = 1;
	  break;
	}
      }
      if (abs_z_2 < 1e-6 || abs(creal(z)) > 1e10 || abs(cimag(z)) > 1e10) {
	attr_result[ix] = expo;
	conv_result[ix] = jx;
	break;
      }
      double complex z_d_1 = 1;
      for (size_t kx = 0; kx < expo-1; ++kx)
	  z_d_1 = z_d_1*z;
      z = z*one_min_expo_inv + expo_inv/z_d_1;

    }
  }
  attractors[item] = attr_result;
  convergences[item] = conv_result;
}

void * write_main(void * args) {
  char * item_done_loc = calloc(size, sizeof(char));
  struct timespec sleep_timespec = {.tv_sec = 0, .tv_nsec = 100000};
  char * intro_string = "P3\n%d %d\n255\n";
  char * template_string = "%3d %3d %3d ";
  char * tmp_string = (char*) malloc(sizeof(char)*30);
  sprintf(tmp_string, "newton_attractors_x%d.ppm", expo);
  FILE * f_attr = fopen(tmp_string,"w+");
  sprintf(tmp_string, "newton_convergences_x%d.ppm", expo);
  FILE * f_conv = fopen(tmp_string,"w+");
  sprintf(tmp_string,intro_string,size,size);
  fwrite(tmp_string,sizeof(char),strlen(tmp_string),f_attr);
  fwrite(tmp_string,sizeof(char),strlen(tmp_string),f_conv);
  char ** colors = (char**) malloc(sizeof(char*)*expo+1);
  for (size_t ix = 0; ix <= expo; ++ix) {
    colors[ix] = (char*) malloc(sizeof(char)*12);
    int color1 = ix * 255.0/expo;
    int color2 = 255 - color1;
    int color3 = 127;
    sprintf(colors[ix],template_string,color1,color2,color3);
  }
  char ** grays = (char**) malloc(sizeof(char*)*52);
  for (size_t ix = 0; ix <= 51; ++ix) {
    grays[ix] = (char*) malloc(sizeof(char)*12);
    int gray = 5 * ix;
    sprintf(grays[ix],template_string,gray,gray,gray);
  }
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
      int * attr_results = attractors[ix];
      int * conv_results = convergences[ix];
      
      for (size_t jx = 0; jx < size; ++jx) {
	int attr_col = attr_results[jx];
	fwrite(colors[attr_col],sizeof(char),strlen(colors[attr_col]),f_attr);
	int conv_col = conv_results[jx];
	if (conv_col > 51) {
	  conv_col = 51;
	}
       	fwrite(grays[conv_col],sizeof(char),strlen(grays[conv_col]),f_conv); 
      }
      free(attr_results);
      free(conv_results);
      fwrite("\n",sizeof(char),1,f_attr);
      fwrite("\n",sizeof(char),1,f_conv);
    }
  }
}

