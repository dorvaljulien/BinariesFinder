#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <kdtree.h>
#include <argsort.h>
#include <string.h>
#include <unistd.h>
#include "binaries.h"
#include "hdf5.h"

#define DIM 3
#define PI 3.14159265



/* -----------------------------------------------------------------------
 -------------------------------------------------------------------------
 -------------------------------------------------------------------------
 ------------------------------------------------------------------------- */

void ReadData_double(hid_t place, char *set_name, double *data){ 
    hid_t set_id = H5Dopen(place, set_name);
    H5Dread(set_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(set_id);
}


void ReadData_float(hid_t place, char *set_name, float *data){ 
    hid_t set_id = H5Dopen(place, set_name);
    H5Dread(set_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(set_id);
}

void ReadData_int(hid_t place, char *set_name, int *data){ 
    hid_t set_id = H5Dopen(place, set_name);
    H5Dread(set_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(set_id);
}


void ReadAttr_int(hid_t place, char *name, int *value){
    hid_t attr;
    attr = H5Aopen(place, name, H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, value); H5Aclose(attr);
}

void ReadAttr_float(hid_t place, char *name, float *value){
    hid_t attr;
    attr = H5Aopen(place, name, H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_FLOAT, value); H5Aclose(attr);
}

void ReadAttr_char(hid_t place, char *name, char *value){
    hid_t attr;
    attr = H5Aopen(place, name, H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, value); H5Aclose(attr);
}


void ReadAttr_string(hid_t place, char *name, char **value){
    hid_t attr;
    attr = H5Aopen(place, name, H5P_DEFAULT);
    H5Aread(attr, H5T_C_S1, value); H5Aclose(attr);
}


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
			       double G)
{

    int i,j,k,i1,i2;
    float t;

    hid_t set;
    int nstep;
    char precision;

    //Miscellaneous variables
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8;
    double rx, ry, rz, rvx, rvy, rvz, mt, m1, m2, E;
    double cm_r, cm_x, cm_y, cm_z;
    int *n_neighbours=malloc(Nnb*sizeof(int));
    double *d_neighbours=malloc(Nnb*sizeof(double));
    double a,e,rv,V,V_local,m_local,rho_bin,rho_local,loc_r,ratio;
                                
    int tmpd;
    int nbin=0;


 
    /*To build the KDtree, we need a 1d array of coordinates*/
    double *coords=malloc(3*N*sizeof(double));
    for(i=0;i<N;i++){
        coords[i]=x[i];
        coords[N+i]=y[i];
        coords[2*N+i]=z[i];
    }
    // We build the KD tree:
    Tree tree;
    KDtree(coords,N,3,&tree);

    // For each star, we look at the 10 closest neighbours
    for(i=0 ; i < N-1 ; i++){
        (void)nnearest(tree,i,n_neighbours, d_neighbours, Nnb);
        /* We acquired the Nnb closest neighbours.*/
        /* We now go through them to check wether they are bound to i.*/
        for(j=0;j<Nnb;j++){
            rvx=vx[n_neighbours[j]]-vx[i];
            rvy=vy[n_neighbours[j]]-vy[i];
            rvz=vz[n_neighbours[j]]-vz[i];
            rv=sqrt(rvx*rvx + rvy*rvy + rvz*rvz);
	    mt=m[i]+m[n_neighbours[j]];

	    E=ebin_norm(d_neighbours[j],rv,mt,G);
            

	    if( E < 0){ // If the two stars are bound
                if (not_done(n[i],n[n_neighbours[j]],n1,n2,nbin)){ // If this pair wasn't already treated
                    /* We compute the density formed by the two bound stars*/
                    a=get_a(d_neighbours[j],rv,mt,G);
		    /* We compute the eccentricity */
		    rx=x[n_neighbours[j]]-x[i];
		    ry=y[n_neighbours[j]]-y[i];
		    rz=z[n_neighbours[j]]-z[i];
		    e=get_e(rx,ry,rz,rvx,rvy,rvz,a,mt,G);
                    V=(4./3)*PI*a*a*a*(1+e)*(1+e)*(1+e);  
		    /* we take a(1+e) as the radius of the binary sphere */
                    rho_bin=mt/V;
                    
                    /* We now compute the density formed by the Nnb neighbours */
                    m_local=0;
                    // Casertano,Hut 1985:
                    // We don't take the last neighbour:
                    for(k=0;k<Nnb-1;k++) m_local+=m[n_neighbours[k]]; 
                    loc_r=d_neighbours[Nnb-1];
                    V_local=(4./3)*PI*loc_r*loc_r*loc_r;
                    rho_local=m_local/V_local;
		    ratio=rho_bin/rho_local;

                   // We compare local density and binary density
                    if (ratio > density_ratio){ // We validate the binary:
			/* We record the primary and secondary mass */
			i1 = m[i] > m[n_neighbours[j]] ?  i : n_neighbours[j];
			i2 = m[i] <= m[n_neighbours[j]] ?  i : n_neighbours[j];
			m1 = m[i1]; m2 = m[i2];
			if (m1==m2) printf("Error, identical binary components !\n");
			
			/* We record the distance from the center of mass to the origin */
			cm_x = ( m[i]*x[i] + m[n_neighbours[j]]*x[n_neighbours[j]] )/ mt; 
			cm_y = ( m[i]*y[i] + m[n_neighbours[j]]*y[n_neighbours[j]] )/ mt; 
			cm_z = ( m[i]*z[i] + m[n_neighbours[j]]*z[n_neighbours[j]] )/ mt; 
			cm_r = sqrt(cm_x*cm_x +cm_y*cm_y + cm_z*cm_z);

			n1[nbin] = n[i1];
			n2[nbin] = n[i2];
                        a_array[nbin]=a;
                        e_array[nbin]=e;
			mt_array[nbin]=mt;
			m1_array[nbin]=m1;
			m2_array[nbin]=m2;
			cmr_array[nbin]=cm_r;
                        ebin_array[nbin]=E;
			ratio_array[nbin]=ratio;
			p_array[nbin]=2*M_PI*sqrt( a*a*a / ( G * mt )  );

                        nbin++;                            
                        /* If we already found too much binaries, print an error message.*/
                        if( nbin > 2*N){
                            printf("Error: number of binaries exceeds 2*number of stars.\n");
                            printf("please check the density ratio, or consider that your system may be too dynamically cold for this tool.\n");
                            return 0;    
                        }
                    }
                }
            }
        }
	

    }
    //printf("Found %d binaries in total.\n", nbin);
    
    //fclose(file);

    /* free(m); free(n); */
    /* free(x); free(y); free(z);  */
    /* free(vx); free(vy); free(vz); */

    free(n_neighbours); free(d_neighbours);
                                
    return nbin;


}


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
		     double G)
{
    int i,j,k;
    int N,dim;
    float t;

    hid_t set;
    int nstep;
    char precision;

    //file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    ReadAttr_int(grp,"N",&N);

    // We prepare the containers for the data about to be read in the snapshot file
    double *m = malloc(N * sizeof(double));
    /* double *r_data = malloc(3*N * sizeof(double)); */
    /* double **r = malloc(3 * sizeof(*double)); */
    /* for(i=0;i++;i<3) r[0] = r_data + i*N; */
    /* double *v_data = malloc(3*N * sizeof(double)); */
    /* double **v = malloc(3 * sizeof(*double)); */
    /* for(i=0;i++;i<3) v[0] = v_data + i*N; */
    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *z = malloc(N * sizeof(double));
    double *vx = malloc(N * sizeof(double));
    double *vy = malloc(N * sizeof(double));
    double *vz = malloc(N * sizeof(double));
    int *n = malloc(N * sizeof(int));

                                
    int tmpd;
    int nbin=0;

    //****************************************************************************************

    ReadData_int(grp,"n",n);
    ReadData_double(grp,"m",m);
    ReadData_double(grp,"x",x);
    ReadData_double(grp,"y",y);
    ReadData_double(grp,"z",z);
    ReadData_double(grp,"vx",vx);
    ReadData_double(grp,"vy",vy);
    ReadData_double(grp,"vz",vz);

    nbin= collect_binaries_from_data(N,n,m,x,y,z,vx,vy,vz,
				     n1,n2, a_array, e_array, mt_array, m1_array, m2_array, 
				     p_array, cmr_array, ebin_array, ratio_array,
				     density_ratio, Nnb, G);


    free(m); free(n);
    free(x); free(y); free(z); 
    free(vx); free(vy); free(vz);
    //free(n_neighbours); free(d_neighbours);
                                
    return nbin;
}

int not_done(int  ni, int  nj, int *n1, int *n2, int nbin)
{/* This functions checks if the pair n1=j and n2=i was already found as a binary */
 /* If that is the case, it returns 0, if it's a new pair, it returns 1*/
    int k;
    for(k=0;k<nbin;k++){
        //printf("nbin: %d  k: %d\n",nbin,k);
        if ( (n1[k]==nj && n2[k]==ni) ||  (n2[k]==nj && n1[k]==ni) ){
            return 0;
        }
    }
    return 1;
}

/* ------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------ */



int follow_binaries(char *filename,
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
		    double G)
{
    int N, dim, i , j, n;
    double t;
    int nbin,stored,testattr;
    char snapfilename[300];
    char string_nstep[10];

    hid_t  file,attr,step_grp;
    int nstep;
    char *precision_string;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    ReadAttr_int(file,"nstep",&nstep);
    ReadAttr_string(file,"precision",&precision_string);
    /* attr = H5Aopen(file, "nstep", H5P_DEFAULT); */
    /* ret  = H5Aread(attr, H5T_NATIVE_INT, &nstep); H5Aclose(attr); */
    /* attr = H5Aopen(file, "precision", H5P_DEFAULT); */
    /* ret  = H5Aread(attr, H5T_NATIVE_INT, &precision); H5Aclose(attr); */

    step_grp = H5Gopen(file, "0");
    ReadAttr_int(step_grp,"N",&N);
    int nx=nbinmax;
    int nbintot=0,nbinactive=0;
    int rank;


    int *n1_snap=malloc(2*N*sizeof(int));
    int *n2_snap=malloc(2*N*sizeof(int));
    double *a_snap=malloc(2*N*sizeof(double));
    double *e_snap=malloc(2*N*sizeof(double));
    double *mt_snap=malloc(2*N*sizeof(double));
    double *m1_snap=malloc(2*N*sizeof(double));
    double *m2_snap=malloc(2*N*sizeof(double));
    double *cmr_snap=malloc(2*N*sizeof(double));
    double *p_snap=malloc(2*N*sizeof(double));
    double *ebin_snap=malloc(2*N*sizeof(double));
    double *ratio_snap=malloc(2*N*sizeof(double));
    // We collect the binaries of the first snapshot to initialize the arrays
    nbin = collect_binaries(step_grp, 
			    n1_snap, 
			    n2_snap, 
			    a_snap, 
			    e_snap, 
			    mt_snap, 
			    m1_snap, 
			    m2_snap, 
			    cmr_snap, 
			    p_snap,
			    ebin_snap, 
			    ratio_snap,
			    density_ratio,
			    Nnb,
			    G);
    printf("n_step = 0   |    %d active detected binaries\n",nbin);

    H5Gclose(step_grp);

    // We initialize the arrays
    for( i=0; i<nbin; i++ ){
	n1[i*nstep]=n1_snap[i];
	n2[i*nstep]=n2_snap[i];
	a_array[i*nstep]=a_snap[i];
	e_array[i*nstep]=e_snap[i];
	mt_array[i*nstep]=mt_snap[i];
	m1_array[i*nstep]=m1_snap[i];
	m2_array[i*nstep]=m2_snap[i];
	cmr_array[i*nstep]=cmr_snap[i];
	p_array[i*nstep]=p_snap[i];
    }
    nbintot=nbin;

    /* We go through each snapshot, collect the binaries, search the n1 and n2 
       arrays for the corresponding binaries, then store them where they 
       belong: after their previous occurences, or at a brand new spot when 
       it just appeared*/
    for( n=1; n < nstep; n++ ){
	//printf("Reading snap_%7.7d:  ",n+1);
	sprintf(string_nstep, "%d", n);
	step_grp = H5Gopen(file,string_nstep);
	nbin = collect_binaries(step_grp, 
				n1_snap, 
				n2_snap, 
				a_snap, 
				e_snap, 
				mt_snap, 
				m1_snap, 
				m2_snap, 
				cmr_snap,
				p_snap,
				ebin_snap,
				ratio_snap,
				1,
				Nnb,
				G);
	H5Gclose(step_grp);

	/* for( i=0; i<5; i++ ){ */
	/*     printf("n1: %d  n2: %d   a: %lf\n",n1_snap[i], n2_snap[i], a_snap[i]); */
	/* } */

	for( i=0; i<nbin; i++ ){ // We go through each binary we found in the new snapshot
	    rank=0; // This variable keep track of which timestep we should use as a start. See end comments.
	    stored=0;
	    for( j=0; j<nbintot; j++ ){
		while(n1[j*nstep+rank]==0 && rank<n){
		    rank++;
		}// If rank reached n-1, it means this binary has not been stored yet
		if ( n1[j*nstep+rank]==n1_snap[i] &&  n2[j*nstep+rank]==n2_snap[i]){
		    // This binary is pre-existing, it doesn't need to pass a new density test
		    //printf("Found existing binary ! n1 = %d   n2 = %d\n", n1_snap[i], n2_snap[i]);
		    n1[j*nstep+n]=n1_snap[i];
		    n2[j*nstep+n]=n2_snap[i];
		    e_array[j*nstep+n]=e_snap[i];
		    a_array[j*nstep+n]=a_snap[i];
		    mt_array[j*nstep+n]=mt_snap[i];
		    m1_array[j*nstep+n]=m1_snap[i];
		    m2_array[j*nstep+n]=m2_snap[i];
		    cmr_array[j*nstep+n]=cmr_snap[i];
		    p_array[j*nstep+n]=p_snap[i];
		    stored=1;
		    break;
		}
	    }
	    if(!stored){ // It means we went through all binaries without finding the one we were searching for
		if(ratio_snap[i]>density_ratio){ // We validate this as a valid new binary
		    n1[nbintot*nstep+n]=n1_snap[i];
		    n2[nbintot*nstep+n]=n2_snap[i];
		    e_array[nbintot*nstep+n]=e_snap[i];
		    a_array[nbintot*nstep+n]=a_snap[i];
		    mt_array[nbintot*nstep+n]=mt_snap[i];
		    m1_array[nbintot*nstep+n]=m1_snap[i];
		    m2_array[nbintot*nstep+n]=m2_snap[i];
		    cmr_array[nbintot*nstep+n]=cmr_snap[i];
		    p_array[nbintot*nstep+n]=p_snap[i];
		    nbintot++;
		    if(nbintot>nbinmax){
			printf("\nError: too many binaries ! If you're not running on a lot of snapshots,\nyou should probably increase the density ratio as you're likely \npicking up weakly bound stars. If you are using this tool on many snapshots,\nyou could increase nbinmax and recompile.");
			return 0;
		    }
		}
	    }
	}// End of loop over binaries found in each file
	/* printf("n = %d\n",n); */
	nbinactive=0;
	for( i=0; i<nbintot; i++ ){
	    if ( a_array[i*nstep+n]>0) nbinactive++;
	}
	printf("n_step = %d   |    %d active detected binaries over %d\n",n,nbinactive,nbintot);
    }// End of file loop
    printf("\n");
    H5Fclose(file);

    free(n1_snap);  free(n2_snap);
    free(a_snap); free(e_snap);
    free(mt_snap); free(m1_snap); free(m2_snap);
    free(cmr_snap); free(ebin_snap); free(ratio_snap);

    return nbintot; 
}





/* ------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------ */

char* concat(char *s1,char *s2){
    char *result=malloc(strlen(s1)+strlen(s2)+1);
    strcpy(result,s1);
    strcat(result,s2);
    return result;
}

char* snapfilepath(int n, char* path)
{
    char *filename;
    asprintf(&filename,"snap_%7.7d",n);
    char *result=concat(path,filename);
    return result;
}


int locate(int N, int *array, int target)
{
    int i;
    for(i=0;i<N;i++){
        if(array[i]==target) return i;
    }
    return -1;
}

/* ------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------ */

double get_a(double r, double v, double mt, double G)
{//    """This function takes the relative radius and velocity vector of a system, the total  mass, and returns the semi-major axis.
    return r / ( 2. - r*v*v / (G*mt) );
}    

double get_e(double rx, double ry, double rz, double vx, double vy, double vz, double a, double mt, double G)
{
    double h1=ry*vz-rz*vy;
    double h2=rz*vx-rx*vz;
    double h3=rx*vy-ry*vx;
    double h = sqrt(h1*h1 + h2*h2 + h3*h3);
    double e2 = 1- h*h / ( a * G * mt );
    return sqrt(e2);
}

double inv_distance(double x1, double x2, double y1, double y2, double z1, double z2)
{
    double distance2;
    distance2= (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
    return 1/ sqrt(distance2);
}


int get_Ek_Ep_Rv(double N, double* x,double* y,double *z,double *V, double *M, double *Ek, double *Ep, double *Rv, double G)
{
    double Ektmp=0, Eptmp=0, invRv=0,invd=0, Mt=0;
    int i,j;
    for( i=0; i<N; i++ ){
	Mt+=M[i];
    }

    for( i=0 ; i < N ; i++ )
    {
        Ektmp=Ektmp+0.5*M[i]*V[i]*V[i];
        for( j=0 ; j < N ; j++ )
        {
            if(i!=j)
            {
                invd=inv_distance(x[i],x[j],y[i],y[j],z[i],z[j]);
                invRv=invRv+M[i]*M[j]*invd;
                if(j >= i+1)
                    {
                        Eptmp=Eptmp-G*M[i]*M[j]*invd;
                    }
            }
        }
    }
    *Rv=Mt*Mt/invRv;
    *Ek=Ektmp;
    *Ep=Eptmp;
    return 0;
}


void get_particles_E(int N, double *m, double *x, double* y, double *z, double *V, double *Ek, double *Ep, double G){
    int i=0,j=0,k=0;
    double Ep_tmp, invd;
    for( i=0 ; i < N ; i++ ){
        Ek[i] = 0.5*m[i]*V[i]*V[i];
	Ep_tmp = 0.0;
        for( j=0 ; j < N ; j++ ){
	    k=0;
            if(i!=j) {
                invd=inv_distance(x[i],x[j],y[i],y[j],z[i],z[j]);
                Ep_tmp += - G*m[i]*m[j]*invd;
            }
	}
	Ep[i] =  Ep_tmp;
	/* E[i] = Ek + Ep; */
    }
}

double ebin(double x, double y, double z, double vx, double vy, double vz, double mt, double G)
// This function returns the binding energy from the relative 3d distance, 3d velocity and total mass in nbody units
{
    double norm_r,norm_v,Ebin;
    norm_r = sqrt( x*x + y*y + z*z );
    norm_v = sqrt( vx*vx + vy*vy + vz*vz );
    Ebin= norm_v*norm_v/2. - G * mt /norm_r;
    return Ebin;
}

double ebin_norm(double r, double v, double mt, double G)
// This function returns the binding energy from the relative 3d distance, 3d velocity and total mass in nbody units
{
    double Ebin;
    Ebin= v*v/2. - G*mt/r;
    return Ebin;
}






int read_snapshot(char *snapfilename, 
		  double *t,
		  int *n,
		  double *m,
		  double *x,
		  double *y,
		  double *z,
		  double *vx,
		  double *vy,
		  double *vz )

{
    int i,j,k;
    int N,dim;
    FILE* file = fopen(snapfilename, "r");
    if( file == NULL)
    {
        printf("%s not found\n",snapfilename);
        return 0;
    }
    //We first read the number of particles, the dimension, and the time of the simulation
    fscanf(file, "%d",&N); 
    fscanf(file, "%d",&dim); 
    fscanf(file, "%lf",t);


    //Miscellaneous variables
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8;
    int tmpd;

    //************************************************************************************

    //We first count the number of lines in the snapshot file
    int ch;
    int nlines=0;
    while (EOF != (ch=fgetc(file)))
        if (ch=='\n')
            nlines++;
    rewind(file);

    fscanf(file, "%d",&N); 
    fscanf(file, "%d",&dim); 
    fscanf(file, "%lf",t);
    if (nlines > 2*N){ //We are dealing with an oldstyle snapshot file
            for( i=0 ; i < N ; i++ ){// We read the masses
                fscanf(file,"%le",&tmp1);
                m[i]=tmp1;
            }
            for( i=0 ; i < N ; i++ ){// We read the positions
                fscanf(file,"%le %le %le",&tmp1,&tmp2,&tmp3);
                x[i]=tmp1; y[i]=tmp2; z[i]=tmp3;
            }
            for( i=0 ; i < N ; i++ ){// We read the velocities
                fscanf(file,"%le %le %le",&tmp1,&tmp2,&tmp3);
                vx[i]=tmp1; vy[i]=tmp2; vz[i]=tmp3;
            }
            for( i=0 ; i < N ; i++ ){// We read the identities of the stars
                fscanf(file,"%d",&tmpd);
                n[i]=tmpd;
            }
        }
    else{// We are dealing with a "modern" snapshot file
            for( i=0 ; i < N ; i++ ){// We read the parameters
                fscanf(file," %le %le %le %le %le %le %le %le ",
                       &tmp1,&tmp2,&tmp3,&tmp4,&tmp5,&tmp6,&tmp7,&tmp8);
                n[i]=(int) tmp1;
                m[i]=tmp2;
                x[i]=tmp3; y[i]=tmp4; z[i]=tmp5;
                vx[i]=tmp6; vy[i]=tmp7; vz[i]=tmp8;
            }
        }
    fclose(file);
    return N;
}






int read_fort10(char *snapfilename, 
		  double *m,
		  double *x,
		  double *y,
		  double *z,
		  double *vx,
		  double *vy,
		  double *vz )

{
    int i,j,k;
    int N,dim;
    FILE* file = fopen(snapfilename, "r");
    if( file == NULL)
    {
        printf("%s not found\n",snapfilename);
        return 0;
    }


    //Miscellaneous variables
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8;
    int tmpd;

    //************************************************************************************

    //We first count the number of lines in the snapshot file
    int ch;
    int nlines=0;
    while (EOF != (ch=fgetc(file)))
        if (ch=='\n')
            nlines++;
    rewind(file);
    N = nlines;
    

    for( i=0 ; i < N ; i++ ){// We read the parameters
	fscanf(file," %le %le %le %le %le %le %le ",
	       &tmp2,&tmp3,&tmp4,&tmp5,&tmp6,&tmp7,&tmp8);
	m[i]=tmp2;
	x[i]=tmp3; y[i]=tmp4; z[i]=tmp5;
	vx[i]=tmp6; vy[i]=tmp7; vz[i]=tmp8;
    }
        
    fclose(file);

    return N;

}

    //**********************************************************************************




int main()
{
    int i,BIG=1000000;
    double *ebin_array=malloc(BIG*sizeof(double));
    int *n1=malloc(BIG*sizeof(int));
    int *n2=malloc(BIG*sizeof(int));
    double *a=malloc(BIG*sizeof(double));
    double *e=malloc(BIG*sizeof(double));
    double *mt=malloc(BIG*sizeof(double));
    double *m1=malloc(BIG*sizeof(double));
    double *m2=malloc(BIG*sizeof(double));
    double *cmr=malloc(BIG*sizeof(double));
    double *p=malloc(BIG*sizeof(double));
    double *ratio=malloc(BIG*sizeof(double));
    /* int nbin=collect_binaries("/home/dorval/heavy_data/nbody6/test/test/snap_0000000" */
    /* 			      ,n1,n2,a,e,mt,m1,m2,cmr, ebin_array,ratio,10,5,1.); */
    int nbin=follow_binaries("/home/dorval/heavy_data/nbody6/test/testhdf5/run64.hdf5",
			     100000,5,10,
			     n1,n2,a,e,mt,m1,m2,cmr,p,1.);
    //printf("nbin: %d\n",nbin);
    /* for() */
    FILE *file=fopen("a_test","w");
    for(i=0;i<nbin;i++){
        fprintf(file,"%lf\n",a[i]);
}
    fclose(file);
    int na,nb;

    printf("nbintot: %d\n",nbin);
    for(i=0;i<10;i++){
        printf("i:%d  n1: %d  n2: %d  a: %lf\n",i,n1[i*11],n2[i*11],a[i*11]);
    }


    free(ebin_array); free(n1); free(n2); free(a); free(mt); free(m1); free(m2); free(cmr);
    free(p); free(ratio);

}


