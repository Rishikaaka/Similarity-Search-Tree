/****** File:  SphereTree4L.c ******/

#include "SphereTree4L.h"


int LPtree_construc(int dim, int ndata, double *data, struct SphereTree4L *ptr_root)  
{
   int  i0,im, i, j, k, kk, cluster_start[KMAX], cluster_size[KMAX], 
        nclusters, nchldrn, node_indx, N1, N2, N3;
   char *cluster_assign;
   double tmp, dist_min, radius, *datum, *buf, *cluster_center, cluster_radius[KMAX] ;

   datum  = (double *)calloc(dim, sizeof(double)) ;
   buf  = (double *)calloc(ndata*dim, sizeof(double)) ;
   cluster_assign = (char *)calloc(ndata, sizeof(char)) ;
   cluster_center =(double *)calloc(KMAX*dim, sizeof(double));

   kk = KMAX ; ////// or use kk = sqrt(sqrt(ndata))
 
   ptr_root->dim = dim;
   ptr_root->NumData = ndata;  
 
   i0=0; im=ndata;
   N1 = bkmeans(kk, dim, i0, im, data,      // dataset info
                cluster_assign, buf, datum, // buffers.  May delete if use a different bkmeans() interface
                cluster_center, cluster_radius, cluster_start, cluster_size) ; // output
   ptr_root->N1 = N1 ;
   ptr_root->L1node_StartChild = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1node_NumChldrn  = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1node_StartDatum = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1node_DataSize   = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1CentersRadii = (double *) calloc(N1*(dim+1), sizeof(double)) ;

   for(k=0; k<(ptr_root->N1); k++) {
      ptr_root->L1node_StartDatum[k] = cluster_start[k] ; 
      ptr_root->L1node_DataSize[k]   = cluster_size[k]  ;  
      for(j=0; j<dim; j++) ptr_root->L1CentersRadii[k*(dim+1)+j] = cluster_center[k*dim+j];
      ptr_root->L1CentersRadii[k*(dim+1)+dim] = cluster_radius[k] ;

//printf("LPT L1 done  N1=%d, radius=%f, size=%d, start=%d\n", ptr_root->N1, 
//  cluster_radius[k], cluster_size[k], cluster_start[k]);
   }


   ptr_root->L2node_StartChild = (int *) calloc(N1*KMAX, sizeof(int)) ;
   ptr_root->L2node_NumChldrn  = (int *) calloc(N1*KMAX, sizeof(int)) ;
   ptr_root->L2node_StartDatum = (int *) calloc(N1*KMAX, sizeof(int)) ;
   ptr_root->L2node_DataSize   = (int *) calloc(N1*KMAX, sizeof(int)) ;
   ptr_root->L2CentersRadii = (double *) calloc(N1*KMAX*(dim+1), sizeof(double)) ;

   node_indx = 0 ;
   for(k=0; k<(ptr_root->N1); k++) { 
      i0= ptr_root->L1node_StartDatum[k];   im= i0 + ptr_root->L1node_DataSize[k] ;

      nclusters=bkmeans(kk, dim, i0, im, data, 
                        cluster_assign, buf, datum,
                        cluster_center,cluster_radius,cluster_start,cluster_size) ; 
      ptr_root->L1node_NumChldrn[k]  = nclusters ;
      ptr_root->L1node_StartChild[k] = node_indx ;

      for(i=0; i<nclusters; i++) { 
         ptr_root->L2node_StartDatum[node_indx] = cluster_start[i] ;
         ptr_root->L2node_DataSize[node_indx]   = cluster_size[i]  ;  
         for(j=0;j<dim;j++) ptr_root->L2CentersRadii[node_indx*(dim+1)+j] = 
                            cluster_center[i*dim+j];
         ptr_root->L2CentersRadii[node_indx*(dim+1)+dim] = cluster_radius[i] ;
         node_indx++ ;

//printf("LPT L2 done  Nclusters=%d, radius=%f, size=%d, start=%d\n", nclusters, 
//       cluster_radius[i], cluster_size[i], cluster_start[i]);
      }

//printf("LPT L2 done  Nclusters=%d, center=(%f, %f), radius=%f, size=%d, start=%d\n", nclusters, 
//cluster_center[k][0],cluster_center[k][1], cluster_radius[k],cluster_size[k],cluster_start[k]);
   }

/*
   ptr_root->N2 = node_indx ;
   N2 = node_indx ;
   ptr_root->L3node_StartChild = (int *) calloc(N2*KMAX, sizeof(int)) ;
   ptr_root->L3node_NumChldrn  = (int *) calloc(N2*KMAX, sizeof(int)) ;
   ptr_root->L3node_StartDatum = (int *) calloc(N2*KMAX, sizeof(int)) ;
   ptr_root->L3node_DataSize   = (int *) calloc(N2*KMAX, sizeof(int)) ;
   ptr_root->L3CentersRadii = (double *)calloc(N2*KMAX*(dim+1), sizeof(double)) ;

   node_indx = 0 ;
   for(k=0; k<ptr_root->N2; k++) { 
      i0 = ptr_root->L2node_StartDatum[k];   im = i0 + ptr_root->L2node_DataSize[k] ;

      nclusters = bkmeans(kk, dim, i0, im, data,
                          cluster_assign, buf, datum,
                          cluster_center, cluster_radius, cluster_start, cluster_size) ; 
      ptr_root->L2node_NumChldrn[k]  = nclusters ;
      ptr_root->L2node_StartChild[k] = node_indx ;
      for(i=0; i<nclusters; i++) { 
         ptr_root->L3node_NumChldrn[node_indx]  = cluster_size[i] ;
         ptr_root->L3node_StartChild[node_indx] = cluster_start[i] ;
         ptr_root->L3node_StartDatum[node_indx] = cluster_start[i] ;
         ptr_root->L3node_DataSize[node_indx]   = cluster_size[i]  ;  
         for(j=0; j<dim; j++) ptr_root->L3CentersRadii[node_indx*(dim+1)+j] = 
                              cluster_center[i*dim+j];
         ptr_root->L3CentersRadii[node_indx*(dim+1)+dim] = cluster_radius[i] ;
         node_indx++ ;

//printf("LPT L3 done  Nclusters=%d, radius=%f, size=%d, start=%d\n", nclusters, 
//  cluster_radius[i], cluster_size[i], cluster_start[i]);
      }
   }
   ptr_root->N3 = node_indx ;
*/

#if 0
   FILE *L1_radius;
   FILE *L1_size;
   FILE *L1_r_s;
   L1_radius = fopen("LP_L1_radius.xvg", "w");
   L1_size   = fopen("LP_L1_size.xvg", "w");
   L1_r_s    = fopen("LP_L1_r_s.xvg", "w");

   for(k=0; k<ptr_root->N1; k++) {
      fprintf(L1_radius, "%d   %f\n", k, ptr_root->L1CentersRadii[k*(DIM+1)+DIM] );
      fprintf(L1_size, "%d    %d\n", k, ptr_root->L1node_DataSize[k]);
      fprintf(L1_r_s, "%d    %f\n", k, ptr_root->L1CentersRadii[k*(DIM+1)+DIM]/ptr_root->L1node_DataSize[k]);

   }

   fclose(L1_radius);
   fclose(L1_size);
   fclose(L1_r_s);

   FILE *L2_radius;
   FILE *L2_size;
   FILE *L2_r_s;
   L2_radius = fopen("LP_L2_radius.xvg", "w");
   L2_size   = fopen("LP_L2_size.xvg", "w");
   L2_r_s    = fopen("LP_L2_r_s.xvg", "w");

   for(k=0; k<ptr_root->N2; k++) {

      fprintf(L2_radius, "%d   %f\n", k, ptr_root->L2CentersRadii[k*(DIM+1)+DIM] );
      fprintf(L2_size, "%d    %d\n", k, ptr_root->L2node_DataSize[k]);
      fprintf(L2_r_s, "%d    %f\n", k, ptr_root->L2CentersRadii[k*(DIM+1)+DIM]/ptr_root->L2node_DataSize[k]);

   }
   fclose(L2_radius);
   fclose(L2_size);
   fclose(L2_r_s);
/*
   FILE *L3_radius;
   FILE *L3_size;
   FILE *L3_r_s;
   L3_radius = fopen("LP_L3_radius.xvg", "w");
   L3_size   = fopen("LP_L3_size.xvg", "w");
   L3_r_s    = fopen("LP_L3_r_s.xvg", "w");

   for(k=0; k<ptr_root->N3; k++) {
      fprintf(L3_radius, "%d   %f\n", k, ptr_root->L3CentersRadii[k*(DIM+1)+DIM] );
      fprintf(L3_size, "%d    %d\n", k, ptr_root->L3node_DataSize[k]);
      fprintf(L3_r_s, "%d    %f\n", k, ptr_root->L3CentersRadii[k*(DIM+1)+DIM]/ptr_root->L3node_DataSize[k]);

   }

   fclose(L3_radius);
   fclose(L3_size);
   fclose(L3_r_s);
*/
# endif 


   /*** free memory: buf, cluster_assign ***/
   free(datum);
   free(buf); 
   free(cluster_assign); 
   free(cluster_center);     

   /*** free_memory: ptr_root->ChldrnCentersRadii, ptr_root->ChildNodeArrays  in main() ***/
} /************ End of function LPtree_construc() ************/












int SStree_construc(int dim, int ndata, double *data, struct SphereTree4L *ptr_root) {

   int  i0,im, i, j, k, kk, child_start, child_end,
        cluster_start[KMAX], cluster_size[KMAX], 
        nclusters, nchldrn, node_indx, N1, N2, N3;
   char *cluster_assign;
   double tmp, dist_min, radius, *datum, *center, *buf,
          *cluster_center, cluster_radius[KMAX] ;

   datum  = (double *)calloc(dim, sizeof(double)) ;
   current_center= (double *)calloc(dim, sizeof(double)) ;
   child_center  = (double *)calloc(dim, sizeof(double)) ;
   buf  = (double *)calloc(ndata*dim, sizeof(double)) ;
   cluster_assign = (char *)calloc(ndata, sizeof(char)) ;
   cluster_center =(double *)calloc(KMAX*dim, sizeof(double));

   kk = KMAX ;////// or use kk = sqrt(sqrt(ndata))
 
   ptr_root->dim = dim;
   ptr_root->NumData = ndata;  
 
   i0=0; im=ndata;
   N1 = bkmeans(kk, dim, i0, im, data,
                cluster_assign, buf, datum,
                cluster_center, cluster_radius, cluster_start, cluster_size) ;

   ptr_root->N1 = N1 ;
   ptr_root->L1node_StartChild = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1node_NumChldrn  = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1node_StartDatum = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1node_DataSize   = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1CentersRadii = (double *) calloc(N1*(dim+1), sizeof(double)) ;

   for(k=0; k<(ptr_root->N1); k++) {
      ptr_root->L1node_StartDatum[k] = cluster_start[k] ; 
      ptr_root->L1node_DataSize[k]   = cluster_size[k]  ;  
      for(j=0; j<dim; j++) ptr_root->L1CentersRadii[k*(dim+1)+j] = cluster_center[k*dim+j];
      ptr_root->L1CentersRadii[k*(dim+1)+dim] = cluster_radius[k] ;

//printf("SST L1 done  N1=%d, radius=%f, size=%d, start=%d\n", ptr_root->N1, 
//  cluster_radius[k], cluster_size[k], cluster_start[k]);

   }


   ptr_root->L2node_StartChild = (int *) calloc(N1*KMAX, sizeof(int)) ;
   ptr_root->L2node_NumChldrn  = (int *) calloc(N1*KMAX, sizeof(int)) ;
   ptr_root->L2node_StartDatum = (int *) calloc(N1*KMAX, sizeof(int)) ;
   ptr_root->L2node_DataSize   = (int *) calloc(N1*KMAX, sizeof(int)) ;
   ptr_root->L2CentersRadii = (double *) calloc(N1*KMAX*(dim+1), sizeof(double)) ;

   node_indx = 0 ;
   for(k=0; k<(ptr_root->N1); k++) { 
      i0= ptr_root->L1node_StartDatum[k];   im= i0 + ptr_root->L1node_DataSize[k] ;

      nclusters = bkmeans(kk, dim, i0, im, data,
                         cluster_assign, buf, datum,
                         cluster_center, cluster_radius, cluster_start, cluster_size) ; 

      ptr_root->L1node_NumChldrn[k]  = nclusters ;
      ptr_root->L1node_StartChild[k] = node_indx ;

      for(i=0; i<nclusters; i++) { 
         ptr_root->L2node_StartDatum[node_indx] = cluster_start[i] ;
         ptr_root->L2node_DataSize[node_indx]   = cluster_size[i]  ;  
         for(j=0;j<dim;j++) ptr_root->L2CentersRadii[node_indx*(dim+1)+j] = 
                            cluster_center[i*dim+j];
         ptr_root->L2CentersRadii[node_indx*(dim+1)+dim] = cluster_radius[i] ;
         node_indx++ ;

//printf("SST L2 done  Nclusters=%d, radius=%f, size=%d, start=%d\n", nclusters, 
//       cluster_radius[i], cluster_size[i], cluster_start[i]);

      }

//printf("SST L2 done  Nclusters=%d, center=(%f, %f), radius=%f, size=%d, start=%d\n", nclusters, 
//cluster_center[k][0], cluster_center[k][1],  cluster_radius[k], cluster_size[k], cluster_start[k]);

   }

/*
   ptr_root->N2 = node_indx ;
   N2 = node_indx ;
   ptr_root->L3node_StartChild = (int *) calloc(N2*KMAX, sizeof(int)) ;
   ptr_root->L3node_NumChldrn  = (int *) calloc(N2*KMAX, sizeof(int)) ;
   ptr_root->L3node_StartDatum = (int *) calloc(N2*KMAX, sizeof(int)) ;
   ptr_root->L3node_DataSize   = (int *) calloc(N2*KMAX, sizeof(int)) ;
   ptr_root->L3CentersRadii = (double *)calloc(N2*KMAX*(dim+1), sizeof(double)) ;

   node_indx = 0 ;
   for(k=0; k<ptr_root->N2; k++) { 
      i0 = ptr_root->L2node_StartDatum[k];   im = i0 + ptr_root->L2node_DataSize[k] ;

      nclusters = bkmeans(kk, dim, i0, im, data,
                          cluster_assign, buf, datum,
                          cluster_center, cluster_radius, cluster_start, cluster_size) ; 

      ptr_root->L2node_NumChldrn[k]  = nclusters ;
      ptr_root->L2node_StartChild[k] = node_indx ;
      for(i=0; i<nclusters; i++) { 
         ptr_root->L3node_NumChldrn[node_indx]  = cluster_size[i] ;
         ptr_root->L3node_StartChild[node_indx] = cluster_start[i] ;
         ptr_root->L3node_StartDatum[node_indx] = cluster_start[i] ;
         ptr_root->L3node_DataSize[node_indx]   = cluster_size[i]  ;  
         for(j=0;j<dim;j++) ptr_root->L3CentersRadii[node_indx*(dim+1)+j] = 
                            cluster_center[i*dim+j];
         ptr_root->L3CentersRadii[node_indx*(dim+1)+dim] = cluster_radius[i] ;
         node_indx++ ;

//printf("LPT L3 done  Nclusters=%d, radius=%f, size=%d, start=%d\n", nclusters, 
//  cluster_radius[i], cluster_size[i], cluster_start[i]);

      }
   }
   ptr_root->N3 = node_indx ;

*/
/*
   for(i=0; i<ptr_root->N2; i++) {     /*** Re-adjust the radius of each L2 node ***/
 /*     for(j=0; j<dim; j++) current_center[j] = ptr_root->L2CentersRadii[i*(dim+1)+j] ; 
 	
      child_start = ptr_root->L2node_StartChild[i] ;
      child_end = child_start + ptr_root->L2node_NumChldrn[i] ;
      for(k=child_start; k<child_end; k++) {
         for(j=0; j<DIM; j++) child_center[j]= ptr_root->L3CentersRadii[k*(dim+1)+j] ;
         dist_ss = calc_dist_square(dim, current_center, child_center) ;
         dist_ss = sqrt(dist_ss);
         tmp = ptr_root->L3CentersRadii[k*(dim+1)+dim] + dist_ss ;
         if(ptr_root->L2CentersRadii[i*(dim+1)+dim] < tmp )   
            ptr_root->L2CentersRadii[i*(dim+1)+dim] = tmp ;   
      }
   }
*/
   for(i=0; i<ptr_root->N1; i++) {    /*** Re-adjust the radius of each L1 node ***/
      for(j=0; j<DIM; j++) current_center[j] = ptr_root->L1CentersRadii[i*(dim+1)+j] ; 

      child_start = ptr_root->L1node_StartChild[i] ;
      child_end = child_start + ptr_root->L1node_NumChldrn[i] ;
      for(k=child_start; k<child_end; k++) {
         for(j=0; j<DIM; j++) child_center[j]= ptr_root->L2CentersRadii[k*(dim+1)+j] ;
         dist_ss = calc_dist_square(dim, current_center, child_center) ;
         dist_ss=sqrt(dist_ss);
         tmp = ptr_root->L2CentersRadii[k*(dim+1)+dim] + dist_ss ;
         if(ptr_root->L1CentersRadii[i*(dim+1)+dim] < tmp)   
            ptr_root->L1CentersRadii[i*(dim+1)+dim] = tmp ; 
      }
   }

#if 0
   FILE *L1_radius;
   FILE *L1_size;
   FILE *L1_r_s;
   L1_radius = fopen("SST_L1_radius.xvg", "w");
   L1_size   = fopen("SST_L1_size.xvg", "w");
   L1_r_s    = fopen("SST_L1_r_s.xvg", "w");

   for(k=0; k<ptr_root->N1; k++) {
      fprintf(L1_radius, "%d   %f\n", k, ptr_root->L1CentersRadii[k*(dim+1)+dim] );
      fprintf(L1_size, "%d    %d\n", k, ptr_root->L1node_DataSize[k]);
      fprintf(L1_r_s, "%d    %f\n", k, ptr_root->L1CentersRadii[k*(DIM+1)+DIM]/ptr_root->L1node_DataSize[k]);
   }

   fclose(L1_radius);
   fclose(L1_size);
   fclose(L1_r_s);


   FILE *L2_radius;
   FILE *L2_size;
   FILE *L2_r_s;
   L2_radius = fopen("LP_L2_radius.xvg", "w");
   L2_size   = fopen("LP_L2_size.xvg", "w");
   L2_r_s    = fopen("LP_L2_r_s.xvg", "w");

   for(k=0; k<ptr_root->N2; k++) {
      fprintf(L2_radius, "%d   %f\n", k, ptr_root->L2CentersRadii[k*(DIM+1)+DIM] );
      fprintf(L2_size, "%d    %d\n", k, ptr_root->L2node_DataSize[k]);
      fprintf(L2_r_s, "%d    %f\n", k, ptr_root->L2CentersRadii[k*(DIM+1)+DIM]/ptr_root->L2node_DataSize[k]);
   }

   fclose(L2_radius);
   fclose(L2_size);
   fclose(L2_r_s);

/*
   FILE *L3_radius;
   FILE *L3_size;
   FILE *L3_r_s;
   L3_radius = fopen("LP_L3_radius.xvg", "w");
   L3_size   = fopen("LP_L3_size.xvg", "w");
   L3_r_s    = fopen("LP_L3_r_s.xvg", "w");

   for(k=0; k<ptr_root->N3; k++) {
      fprintf(L3_radius, "%d   %f\n", k, ptr_root->L3CentersRadii[k*(DIM+1)+DIM] );
      fprintf(L3_size, "%d    %d\n", k, ptr_root->L3node_DataSize[k]);
      fprintf(L3_r_s, "%d    %f\n", k, ptr_root->L3CentersRadii[k*(DIM+1)+DIM]/ptr_root->L3node_DataSize[k]);

   }

   fclose(L3_radius);
   fclose(L3_size);
   fclose(L3_r_s);*/
# endif 




   /*** free memory: buf, cluster_assign ***/
   free(datum);
   free(buf) ; 
   free(cluster_assign); 
   free(cluster_center); 
   free(current_center);    
   free(child_center) ;

   /*** free_memory: ptr_root->ChldrnCentersRadii, ptr_root->ChildNodeArrays  in main() ***/

} /************ End of function SStree_construc() ************/


















int allneighbors_search(int dim, int ndata, double *data, struct SphereTree4L *ptr_root, 
                         double *query, double delta, int *num_outdata, int *outdata) 
{
   int    i, j, k, kk, n1_nearclusters, n2_nearclusters, n3_nearclusters,count,
          *nearclusters1, *nearclusters2, *nearclusters3, cindx, start, end, offset ;
   double dist,radius, *center, *datum;

   center= (double *) calloc(dim, sizeof(double)) ;
   datum = (double *) calloc(dim, sizeof(double)) ;

   count = 0 ;
   nearclusters1 = (int *)calloc( ptr_root->N1, sizeof(int) ) ;

   offset = 0 ;
   for(k=0; k<ptr_root->N1; k++) {
      for(j=0; j<dim; j++) center[j] = ptr_root->L1CentersRadii[k*(dim+1)+j] ;
      radius = ptr_root->L1CentersRadii[k*(dim+1)+dim] ;
      dist = calc_dist_square(dim, query, center);
      dist = sqrt(dist);
      count++ ;
      if(dist < radius + delta) {
         nearclusters1[offset] = k ;
         offset++ ;
      }
   }
   n1_nearclusters = offset ;
   nearclusters2 = (int *)calloc( n1_nearclusters*KMAX, sizeof(int) ) ;

   printf("L1 nodes checked = %d\n", count);
   //printf("n1_clusters = %d\n",n1_nearclusters);



   offset = 0 ;
   for(i=0; i<n1_nearclusters; i++) {
      cindx = nearclusters1[i] ;
      start = ptr_root->L1node_StartChild[cindx] ;
      end = start + ptr_root->L1node_NumChldrn[cindx] ;
      for(k=start; k<end; k++) {
         for(j=0; j<dim; j++) center[j] = ptr_root->L2CentersRadii[k*(dim+1)+j] ;
         radius = ptr_root->L2CentersRadii[k*(dim+1)+dim] ;
         dist = calc_dist_square(dim, query, center);
         dist = sqrt(dist);
         count++ ;
         if(dist < radius + delta) {
            nearclusters2[offset] = k ;
            offset++ ;
         }
      }
   }
   n2_nearclusters = offset ;
   nearclusters3 = (int *)calloc( n2_nearclusters*KMAX, sizeof(int) ) ;


   printf("L2 nodes checked = %d\n", count);
  // printf("n2_clusters = %d\n",n2_nearclusters);


   offset = 0 ;
   for(i=0; i<n2_nearclusters; i++) {
      cindx = nearclusters2[i] ;
      start = ptr_root->L2node_StartChild[cindx] ;
      end = start + ptr_root->L2node_NumChldrn[cindx] ;
      for(k=start; k<end; k++) {
         for(j=0; j<dim; j++) center[j]= ptr_root->L3CentersRadii[k*(dim+1)+j] ;
         radius = ptr_root->L3CentersRadii[k*(dim+1)+dim] ;
         dist = calc_dist_square(dim, query, center);
         dist = sqrt(dist);
         count++ ;
         if(dist < radius + delta) {
            nearclusters3[offset] = k ;
            offset++ ;
         }
      }
   }
/*   n3_nearclusters = offset ;
   outdata = (int *)calloc( n3_nearclusters*LEAFSIZEMAX, sizeof(int) ) ;


   printf("L3 nodes checked = %d\n", count);
   //printf("n3_clusters = %d\n",n3_nearclusters);

   offset = 0 ;
   for(i=0; i<n3_nearclusters; i++) {
      cindx = nearclusters3[i] ;
      start = ptr_root->L3node_StartDatum[cindx] ;
      end = start + ptr_root->L3node_DataSize[cindx] ;
      for(k=start; k<end; k++) {
         for(j=0; j<dim; j++) datum[j]= data[k*dim+j] ;
         dist = calc_dist_square(dim, query, datum);
         dist = sqrt(dist);
         count++ ;
         if(dist < delta) {
            outdata[offset] = k ;
            offset++ ;
         }
      }
   }*/
   *num_outdata = offset ;
   printf("n_outdata:=%d,  count=%d\n", *num_outdata, count);
   
   free(nearclusters1);
   free(nearclusters2);
   free(nearclusters3);
   free(outdata);
   free(center);
   free(datum) ;

   return count ;
} /****** End of function allneighbors_search() ******/








