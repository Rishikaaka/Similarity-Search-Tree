#if !defined BKMEANS_H
#define BKMEANS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define min(x,y)        ((x)>(y) ? (y) : (x))

double calc_dist_square(int dim, double *datum1, double *datum2);
#if 0
double calc_dist_square(int dim, double *datum1, double *datum2) { 
   unsigned short i, j;
   double tmp, dist = 0.0;

   for(i=0; i<dim; i++) { tmp = datum1[i]-datum2[i] ;   dist += tmp*tmp ; }
   return dist;   /* Square of Euclidean distance */
}  /****** end of calc_dist_square  ******/
#endif



int two_means(int iterat_limit,
     int dim, int i0, int im, double *data, 
     int *cluster_assign, double *datum, double *center0,
     double *radius_pt, double *center, int start[2], int size[2], double ssd[2]);
/*************************************************************************
 * array sizes:                                                          *
 *       cluster_assign[ndata], datum[dim], center[2*dim]        *
 * Input:                                                                *
 *       center[0] - centroid of the cluster                     *
 *       radius_pt[0] - farthest pt to the centroid                      *
 * Output:                                                               *
 *       radius_pt, center, start, size, ssd     *
 * buffers: cluster_assign[], datum, cluster_center0[]                   *
 * radius_pt[2*dim]: - the radius pt of the cluster[k]                   *
 *************************************************************************/






/******************************************************************************
 kk :               number of clusters, i.e. the K in K-mean.
 cluster_center[kk*dim]: input  -- stores initial kk centers
                         output -- stores kk centers
 cluster_radius[kk]:output -- the radius of each output cluster
 cluster_start[kk]: output -- the index of the 1st in each cluster
 cluster_size[kk]:  output -- the num of datapoints in each cluster
 dataset_size and mem_capacity in unit of dim*sizeof(double),          
                               i.e. in unit of data items               
*******************************************************************************/
int kmeans(int iterat_limit, int kk, 
    int dim, int i0, int im, double *data,
    int *cluster_assign, double *datum, double *cluster_center0, double *radius_pt,
    double *cluster_center, double *cluster_radius, int *cluster_start,
              int *cluster_size, double *cluster_ssd);




int bkmeans(int iterat_limit, int kk, 
    int dim, int ndata, double *data,
    int *cluster_assign, double *datum,
    double *cluster_center, double *cluster_radius, int *cluster_size, double *cluster_ssd);


#endif
