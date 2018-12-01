/******** File: main.c **********/

#include "SphereTree4L.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#define DATASET		1
#define PI 3.14


void main()
{
   int i, j, k, count, num_outdata, *outdata, kk=KMAX, dim, ndata ,nclusters;
   double delta, *data, *query;

   FILE *test_output_file;
   test_output_file = fopen("Test_LPT.out", "w");


#if(DATASET==1)
   double length, U1, U2, *buf, tmp;
   double *centers[100], stdevs[100],
          clustersizes[100], cluster_start, cluster_end;

   dim = 128;   ndata = 60123456;   nclusters=100 ;
   data= (double *) calloc(ndata*dim, sizeof(double));
   buf = (double *) calloc(ndata, sizeof(double));
   for(k=0; k<100; k++) centers[k] = calloc(dim, sizeof(double)) ;

   for(k=0; k<nclusters/10; k++) {
      clustersizes[k]    =   10*ndata/(2*nclusters);
      clustersizes[k+nclusters/10]   = 10*ndata/(4*nclusters);
      clustersizes[k+2*nclusters/10] = 10*ndata/(8*nclusters);
      clustersizes[k+3*nclusters/10] = 10*ndata/(16*nclusters);
      clustersizes[k+4*nclusters/10] = 10*ndata/(32*nclusters);
      clustersizes[k+5*nclusters/10] = 10*ndata/(64*nclusters);
      clustersizes[k+6*nclusters/10] = 10*ndata/(128*nclusters);
      clustersizes[k+7*nclusters/10] = 10*ndata/(256*nclusters);
      clustersizes[k+8*nclusters/10] = 10*ndata/(512*nclusters);
      clustersizes[k+9*nclusters/10] = 10*ndata/(512*nclusters);
   }

   srand(1) ;
   for(k=0; k<nclusters; k++) {
      for(j=0; j<dim; j++)
         centers[k][j] = 40.0*((double)rand()/RAND_MAX - 0.5);
      stdevs[k] = 2.0 + 3.0*((double)rand()/RAND_MAX); /* stdev in (2, 5) */
   }

   tmp = 1.0/RAND_MAX ; /*** Generating data[i] in [-0.5, 05] ***/
   for(i=0;i<ndata*dim;i++) data[i] = (double)rand()*tmp - 0.5;

   for(i=0;i<ndata;i++) {/*** Make each datum unit length ***/
      length = 0.0 ;
      for(j=0; j<dim; j++) length += data[i*dim+j]*data[i*dim+j];
      length = 1.0/sqrt(length) ;
      for(j=0; j<dim; j++) data[i*dim+j] = data[i*dim+j]*length ;
   }

   for(i=0;i<ndata/2;i++){
      U1=(double)rand()/RAND_MAX;
      U2=(double)rand()/RAND_MAX;
      buf[2*i]=sqrt(-2*log(U1))*cos(2*PI*U2);
      buf[2*i+1]=sqrt(-2*log(U1))*sin(2*PI*U2);
   }

   cluster_start = 0 ;
   for(k=0; k<nclusters; k++) {
      cluster_end = cluster_start + clustersizes[k] ;
      for(i=cluster_start; i<cluster_end; i++) {
         buf[i] = stdevs[k]*buf[i] ;
         for(j=0; j<dim; j++) {
            data[i*dim+j] *= buf[i] ;
            data[i*dim+j] += centers[k][j];
         }
      }
      cluster_start = cluster_end ;
   }

   free(buf) ;
   for(k=0; k<200; k++) free(centers[k]) ;
#endif


#if (DATASET == 2) /*** Read HIGGS.dat binary file of double floating-pt data ***/
     FILE *fp;

     fp=fopen("~/Datasets/HIGGS/HIGGS.dat","rb");
     dim = 29;   ndata = 11000000;
     data = (double *) calloc(dim*ndata, sizeof(double)) ;
     if(fp!= NULL) fread(data, sizeof(float), dim*ndata, fp);
     fclose(fp);
#endif


#if 1 /****** Construct LP tree ******/
   struct SphereTree4L root;   clock_t start = clock(); 
   LPtree_construc(dim, ndata, data, &root); 
   clock_t finish = clock();
   double duration = (double)(finish - start) / CLOCKS_PER_SEC; 
  // printf("LPT4L: DIM=%d, KMAX=%d, NDATA=%d, LEAFSIZEMAX=%d, DATASET=%d, delta=%f\n",   
  //        dim,KMAX,ndata,LEAFSIZEMAX,DATASET,delta);
   printf("LPT4L: Used Time: %lf s\n",duration);
#else /****** Construct SS tree ******/
   clock_t start = clock(); 
   SStree_construc(dim, data, data, &root); 
   clock_t finish = clock();
   double duration = (double)(finish - start) / CLOCKS_PER_SEC; 
  // printf("SST4L: DIM=%d, KMAX=%d, NDATA=%d, LEAFSIZEMAX=%d, DATASET=%d, delta=%f\n",   
  //        dim,KMAX,ndata,LEAFSIZEMAX,DATASET,delta);
   printf("SST4L: Used Time: %lf s\n",duration);
#endif



/****** Call search() for queries ******/
  /* delta=0.1;
   query = (double *) calloc(dim, sizeof(double)) ;

   printf("\n Query 1\n");
   query[0] = -0.1; 
   for(j=1; j<dim; j++) query[j] = -query[j-1];
   count= allneighbors_search(dim,ndata,data,&root, query,delta, &num_outdata,outdata);
   //fprintf(Test_result,"Query 1: num_outdata=%d,  Count =%d\n", num_outdata, count);


   printf("\n Query 2\n");
   query[0] = -0.05; 
   for(j=1; j< dim; j++) query[j] = -query[j-1]*j;
   count= allneighbors_search(dim,ndata,data,&root, query,delta, &num_outdata, outdata);
   //fprintf(Test_result,"Query 2: num_outdata=%d,  Count =%d\n", num_outdata, count);


   printf("\n Query 3\n");
   for (j=0; j<dim; j++) query[j] = 0.1;
   count= allneighbors_search(dim,ndata,data,&root, query,delta, &num_outdata, outdata);
   //fprintf(Test_result,"Query 3: num_outdata=%d,  Count =%d\n", num_outdata, count);


   printf("\n Query 4\n");
   for (j=0; j<dim; j++) query[j] = 0.0;
   count=allneighbors_search(dim,ndata,data,&root, query,delta, &num_outdata, outdata);
   //fprintf(Test_result,"Query 4: num_outdata=%d,  Count =%d\n",num_outdata, count);

   free(query) ;

   close(test_output_file);
*/

   free(root.L1node_StartChild)  ;
   free(root.L1node_NumChldrn ) ;
   free(root.L1node_StartDatum) ;
   free(root.L1node_DataSize) ;
   free(root.L1CentersRadii);

   free(root.L2node_StartChild) ;
   free(root.L2node_NumChldrn) ;
   free(root.L2node_StartDatum);
   free(root.L2node_DataSize);
   free(root.L2CentersRadii );
/*
   free(root.L3node_StartChild) ;
   free(root.L3node_NumChldrn);
   free(root.L3node_StartDatum);
   free(root.L3node_DataSize);
   free(root.L3CentersRadii);
*/
   free(data) ;
 
} /******************** End of main() ********************/

/******************** End of File: main.c **********************/
