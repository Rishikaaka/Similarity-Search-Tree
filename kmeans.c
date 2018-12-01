/****** File: kmeans.c ******/
#include "bkmeans.h"


/******************************************************************************
 kk : number of clusters, i.e. the K in K-mean.
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
    double *cluster_center, double *cluster_radius, 
          int *cluster_start, int *cluster_size, double *cluster_ssd)
{  int i, j, k, k_max, iterations, membership, start0, end0, start1, end1,
       *cluster_size0, position, change, *radius_index, k_max_used, nclusters,
       *local_cluster_start, **local_cluster_size ;
   double  tmp, dist_min;

   nclusters = kk ;
   cluster_size0 = (int *) calloc(kk, sizeof(int)) ;
   radius_index  = (int *) calloc(kk, sizeof(int)) ;

  /****** Start of k-means iterations ******/
   change = 1 ;
   iterations = 0 ;
   while( (iterations < iterat_limit) && (change!=0) ) {
      iterations++ ;
      change = 0 ;
      for(k=0; k<kk; k++) {
         cluster_size0[k] = cluster_size[k] ;
         cluster_size[k]=0 ;
         cluster_radius[k]=0.0;
         for(j=0;j<dim;j++) cluster_center0[k*dim+j]=cluster_center[k*dim+j];
         for(j=0;j<dim;j++) cluster_center[k*dim+j] = 0.0 ;
      }/* cluster_center0 needed for calculating cluster_center at BBB */

      for(i=i0; i<im; i++) {
         memcpy(datum, data+i*dim, dim*sizeof(double));
         dist_min = 987654321012345.0 ; /* Find closest center to datum */
         for(k=0; k<kk; k++) {
            if(cluster_size0[k]>0) {
               tmp = calc_dist_square(dim, datum, cluster_center0+k*dim);
               if(tmp < dist_min) { membership=k ;    dist_min=tmp ; }
            }
         }
         if(cluster_assign[i] != membership) change++ ;
         cluster_assign[i] = membership ;
         cluster_size[membership]++ ;
         for(j=0;j<dim;j++) cluster_center[membership*dim+j]+=datum[j];//BBB
         if(dist_min > cluster_radius[membership]) {
            cluster_radius[membership] = dist_min ;
            radius_index[membership] = i ;
            memcpy(radius_pt+dim*membership, datum, dim*sizeof(double));
         }
      }/*** Each data item has been assigned to a cluster ***/

      k_max=0;   tmp=cluster_radius[0];
      for(k=1; k<kk; k++) {
         if(tmp<cluster_radius[k]){ k_max=k; tmp=cluster_radius[k];}
      }
      k_max_used = 0 ;
      for(k=0; k<kk; k++) {
         if(cluster_size[k] > 0) 
            for(j=0;j<dim;j++) cluster_center[k*dim+j] /= (double)cluster_size[k];
         else if(!k_max_used) {
            k_max_used = 1 ;
            i = radius_index[k_max] ;

            cluster_size[k_max]--;
            cluster_size[k]++ ;
            cluster_assign[i] = k ;
            memcpy(cluster_center+k*dim, radius_pt+k_max*dim, dim*sizeof(double));
         }
      }
   }/*** End of while( iterations < iterat_limit ) ***/


/****** Data Re-ordering ******/
   position = 0 ;
   for(k=0;k<kk;k++) {
      cluster_ssd[k] = 0.0 ;        cluster_radius[k]=0.0;
      cluster_start[k]= position;   position += cluster_size[k];
   }
   for(k=0; k<kk-1; k++) {
      start0= cluster_start[k] ;
      end0  = start0 + cluster_size[k];
      start1= cluster_start[k+1];
      while(start0 < end0) {
         while((start0<end0)&&(cluster_assign[start0]==k)) start0++ ;
         if(start0<end0) {
            while((start1<im)&&(cluster_assign[start1]!=k)) start1++ ;
            if(start1==im) {
               printf("\nError: start1==im.\n");   return 0;
            }
            memcpy(datum, data+start0*dim, dim*sizeof(double)) ;
            memcpy(data+start0*dim, data+start1*dim, dim*sizeof(double)) ;
            memcpy(data+start1*dim, datum, dim*sizeof(double)) ;
            cluster_assign[start1] = cluster_assign[start0] ;
            cluster_assign[start0] = k;
            start0++;   start1++;
         }
      }
   }/****** End of loop for(k=0; k<kk-1; k++). End of Data Re-ordering ******/

   for(k = 0; k < kk; k++) {/*** Calculate sse, radius for every chunk ***/
      cluster_ssd[k] = 0.0 ;
      cluster_radius[k]=0.0;
      start0 = cluster_start[k] ;
      end0 = start0 + cluster_size[k] ;
      for(i=start0; i<end0; i++) {
         tmp = calc_dist_square(dim, data+i*dim, cluster_center+k*dim) ;
         if(tmp > cluster_radius[k]) cluster_radius[k] = tmp ;
         cluster_ssd[k] += tmp ;
      }
      if(cluster_size[k]==0) nclusters-- ;
   }
   for(k=0; k<kk; k++) cluster_radius[k] = sqrt(cluster_radius[k]) ;

   free(cluster_size0) ;   free(radius_index) ;
   return nclusters ;
} /****************** End of function kmeans() ******************/

/****************** End of File kmeans.c ******************/
