/****** File: bkmeans.c ******/
#include "bkmeans.h"


double calc_dist_square(int dim, double *datum1, double *datum2) {
   int    i;
   double tmp, dist = 0.0;

   for(i=0; i<dim; i++) { tmp = datum1[i]-datum2[i] ;   dist += tmp*tmp ; }
   return dist;   /* Square of Euclidean distance */
}  /****** end of calc_dist_square  ******/




/*************************************************************************
 * array sizes:                                                          *
 *       cluster_assign[ndata], datum[dim], center[2*dim]                *
 * Input:                                                                *
 *       center[0] - centroid of the cluster                             *
 *       radius_pt[0] - farthest pt to the centroid                      *
 * Output:                                                               *
 *       radius_pt, center, start, size, ssd                             *
 * buffers: cluster_assign[], datum, center0[]                           *
 * radius_pt[2*dim]: - the radius pt of the cluster[k]                   *
 *************************************************************************/
int two_means(int iterat_limit,
        int dim, int i0, int im, double *data,        /* line of input dataset */
        int *cluster_assign, double *datum, double *center0, /* line of buffer */
        double *radius_pt, double *center, int start[2], int size[2], double ssd[2])
{
   int    i, j, k, i_max, iterations, iterat_limit, change, offset,
          start0, start1, end0, end1 ;
   double tmp, dist_max, dist0, dist1 ;

   for(j=0;j<dim;j++) {/*** Choosing initial pair of centers ***/
      center[dim+j] = radius_pt[j] ;
      center[j] = 2.0*center[j] - radius_pt[j] ;
   } /*** End of choosing initial pair of centers ***/


   change = 1 ;
   iterations = 0 ;
   while( (iterations<iterat_limit) && (change!=0) ) {
      iterations++ ;
      change = 0 ;
      size[0]=0 ;   size[1]=0 ;
      for(j=0;j<dim;j++) {
         center0[j]    = center[j] ;
         center0[dim+j]= center[dim+j];
         center[j]=0.0;  center[dim+j]=0.0;
      }/* center0 needed for calculating center at AAA */

      for(i=i0; i<im; i++) {/****** Passes 2 to iterat_limit ******/
         memcpy(datum, data+i*dim, dim*sizeof(double)) ;
         dist0 = calc_dist_square(dim, datum, center0);
         dist1 = calc_dist_square(dim, datum, center0+dim);
         k = ((dist0 < dist1) ? 0 : 1 ) ;
         if( cluster_assign[i] != k ) change++ ;
         cluster_assign[i] = k ;
         size[k]++ ;
         for(j=0; j<dim; j++) center[k*dim+j]+=datum[j]; // AAA
      } /****** Now, each datum has been assigned to a cluster ******/

      if(size[0] > 0) 
         for(j=0; j<dim; j++) center[j] /= size[0] ;
      if(size[1] > 0)
         for(j=0; j<dim; j++) center[dim+j] /= size[1] ;

   }/****** End of while(iterations < iterat_limit) ******/


   /****** Data Re-ordering and compute sse[k]: Pass 1+iterat_limit ******/
   start[0]=i0;       start[1] = i0+size[0];
   start0= start[0];  end0 = start0 + size[0] ;
   start1= start[1];  end1 = start1 + size[1] ;
   ssd[0]= 0.0 ;    ssd[1] = 0.0 ;
   dist0 = 0.0 ;    dist1 = 0.0 ; /*** dist0, dist1 hold radii of two clusters ***/

   while(start0 < end0) {/*** end0 and end1 MUST NOT EXCEED im ***/
      while( (start0<end0) && (cluster_assign[start0]==0) ) start0++ ;
      while( (start1<end1) && (cluster_assign[start1]==1) ) start1++ ;
      if( (start0<end0) && (start1<end1) ) {
         memcpy(datum, data+start0*dim, dim*sizeof(double)) ;
         memcpy(data+start0*dim, data+start1*dim, dim*sizeof(double)) ;
         memcpy(data+start1*dim, datum, dim*sizeof(double)) ;
         start0++;   start1++;

         for(i=start[0]; i<start0; i++) {
            tmp = calc_dist_square(dim, data+i*dim, center+k*dim) ;
            if(dist0 < tmp) { 
               dists0 = tmp ;
               memcpy(radius_pt, data+i*dim, dim*sizeof(doouble));
            }
            ssd[0] += tmp;
         }
         for(i=start[1]; i<start1; i++) {
            tmp = calc_dist_square(dim, data+i*dim, center+k*dim) ;
            if(dist1 < tmp) { 
               dists1 = tmp ;
               memcpy(radius_pt+dim, data+i*dim, dim*sizeof(doouble));
            }
            ssd[1] += tmp;
         }
      }
      else if( (start0<end0)&&(start1==end1) ) {
         printf("Error: start0<end0 && start1==end1\n");   return 0;
      }
      else if( (start0==end0)&&(start1<end1) ) {
         printf("Error: start0==end0 && start1<end1\n");   return 0;
      }
      else {;}
   }/*** End of Data Re-ordering ***/

   return 1 ;
} /****************** End of function two_means() ******************/












int bkmeans(int iterat_limit, int kk, 
    int dim, int ndata, double *data,       /* line of input dataset */
    int *cluster_assign, double *datum,     /* line of buffers       */
    double *cluster_center, double *cluster_radius, int *cluster_size, double *cluster_ssd)
{ /*** cache_capacity: input values are in bytes ***/
   int    i0, im, i, j, k, k_max, nclusters, start[2], size[2] ;
   double tmp, dist_max, dist, ssd_initial, ssd[2], 
          *center, *radius_pt, *cluster_center0, *cluster_radius_pt ;

   center  = (double *) calloc((2*dim), sizeof(double)) ;
   radius_pt=(double *) calloc((2*dim), sizeof(double)) ;
   cluster_center0  = (double *) calloc((kk*dim), sizeof(double)) ;
   cluster_radius_pt =(double *) calloc((kk*dim), sizeof(double)) ;

   memset(cluster_radius, 0.0, kk*sizeof(double)) ;
   memset(cluster_ssd,    0.0, kk*sizeof(double)) ;
  
/*** Compute centroid of the whole dataset ***/
   memset(cluster_center, 0.0, dim*sizeof(double)) ;
   for(i=i0; i<im; i++)
      for(j=0; j<dim; j++) cluster_center[j] += data[i*dim+j];
   for(j=0; j<dim; j++) cluster_center[j] /= (double)dataset_size ;
/*** End of computing centroid of the whole dataset ***/ 

/*** Compute farthest pt to centroid. Store it in radius_pt[0] ***/
   ssd_initial = 0.0 ;
   for(i=i0; i<im; i++) {
      tmp = calc_dist_square(dim, data+i*dim, cluster_center) ;
      ssd_initial += tmp ;
      if(cluster_radius[0] < tmp) { 
         cluster_radius[0] = tmp ;
         memcpy(cluster_radius_pt, data+i*dim, dim*sizeof(doouble)) ;
      }
   }/*** End of computing farthest pt to centroid ***/


   i0 = 0 ;   im = ndata ;
   two_means(iterat_limit, 
       dim, i0, im, data,
       cluster_assign, datum, cluster_center0,
       cluster_radius_pt,cluster_center,cluster_start,cluster_size,cluster_ssd) ;
   nclusters = 2 ;

   while( nclusters < kk ) {
      tmp = 0.0 ;
      for(k=0; k<nclusters; k++) {  /* Find cluster with largest ssd */
         if( cluster_ssd[k] > tmp ) { tmp=cluster_ssd[k]; k_max=k; }
      }

      /*** Split a cluster into 2 clusters ***/
      i0 = cluster_start[k_max];
      im = i0 + cluster_size[k_max];
      for(j=0;j<dim;j++) radius_pt[j]= cluster_radius_pt[k_max*dim+j] ;
      for(j=0;j<dim;j++) center[j]   = cluster_center[k_max*dim+j] ;
      two_means(iterat_limit,
                dim, i0, im, data,
                cluster_assign, datum, cluster_center0,
                radius_pt, center, start, size, ssd) ;

      cluster_start[k_max]  = start[0] ;
      cluster_start[nclusters]= start[1] ;
      cluster_size[k_max]  = size[0] ;
      cluster_size[nclusters]= size[1] ;
      cluster_ssd[k_max]  = ssd[0] ;
      cluster_ssd[nclusters]= ssd[1] ;
      memcpy(cluster_radius_pt+k_max*dim,    center,    dim*sizeof(double));
      memcpy(cluster_radius_pt+nclusters*dim,center+dim,dim*sizeof(double));
      memcpy(cluster_center+k_max*dim,     center,     dim*sizeof(double)) ;
      memcpy(cluster_center+nclusters*dim, center+dim, dim*sizeof(double)) ;
      nclusters++;
   } /************ End of while() loop ************/

   for(k=0; k<nclusters; k++) {
      if(cluster_size[k] == 0) { printf("In bkmeans: cluster %d empty!\n"); }
   }

   if(nclusters > 2) {
      iterat_limit = 5;   i0 = 0;   im = ndata;
      nclusters = kmeans(iterat_limit, nclusters,
                  dim, i0, im, data,
                  cluster_assign, datum, cluster_center0, cluster_radius_pt,
                  cluster_center,cluster_radius,cluster_start,cluster_size,cluster_ssd);
   }

   free(center); free(radius_pt);  free(cluster_center0); free(cluster_radius_pt);
   return nclusters ;
} /****************** End of bkmeans() ******************/
/****************** End of File bkmeans.c ******************/
