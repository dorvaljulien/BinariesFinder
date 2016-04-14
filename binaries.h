#include "hdf5.h"
char* concat(char *s1,char *s2);
char* snapfilepath(int n, char* path);
int locate(int N, int *array, int target);
double get_a(double r, double v, double mt, double G);
double get_e(double rx, double ry, double rz, double vx, double vy, double vz, double a, double mt, double G);
double inv_distance(double x1, double x2, double y1, double y2, double z1, double z2);
int get_Ek_Ep_Rv( double N, double* x,double* y,double *z,double *V, double *M, double *Ek, double *Ep, double *Rv, double G);
double ebin(double x, double y, double z, double vx, double vy, double vz, double mt, double G);
double ebin_norm(double r, double v, double mt, double G);
int collect_binaries_from_data(int N,
			       int *n,
			       double *m,
			       double *x,
			       double *y,
			       double *z,
			       double *vx,
			       double *vy,
			       double *vz,
			       int *n1,
			       int *n2, 
			       double *a_array, 
			       double *e_array, 
			       double *mt_array, 
			       double *m1_array, 
			       double *m2_array, 
			       double *cmr_array, 
			       double *p_array,
			       double *ebin_array,
			       double *ratio_array,
			       double density_ratio,
			       int Nnb,
			       double G);
int collect_binaries(hid_t grp, 
		     int *n1, 
		     int *n2, 
		     double *a_array, 
		     double *e_array, 
		     double *mt_array, 
		     double *m1_array, 
		     double *m2_array, 
		     double *cmr_array, 
		     double *p_array, 
		     double *ebin_array,
		     double *ratio_array,
		     double density_ratio,
    		     int Nnb,
		     double G);
int not_done(int  ni, int  nj, int *n1, int *n2, int nbin);
int follow_binaries(char *path,
		    int nbinmax,
		    int Nnb,
		    double density_ratio,
		    int *n1,
		    int *n2,
		    double *a_array, 
		    double *e_array, 
		    double *mt_array, 
		    double *m1_array, 
		    double *m2_array, 
		    double *cmr_array, 
		    double *p_array,
		    double G);
int read_snapshot(char *snapfilename, 
		  double *t,
		  int *n,
		  double *m,
		  double *x,
		  double *y,
		  double *z,
		  double *vx,
		  double *vy,
		  double *vz );
int read_fort10(char *snapfilename, 
		  double *m,
		  double *x,
		  double *y,
		  double *z,
		  double *vx,
		  double *vy,
		  double *vz );


