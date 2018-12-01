#if !defined LPT4L_H
#define LPT4L_H

#define KMAX        128  /* Maximum number of children of each tree node */




/****************************************************************
                          o                                       L0

    o              o                o                o            L1
   
   o o o          o o o o         o o o o        o o o o o        L2

 ooo oo ooo     ooo ooo oo oo   oo oo oo oo    ooo oo oo oo oo    L3=Leaves=data

 In this example, NumLevels=4, N1=4, N2=16, N3=37
    L1node_NumChldrn[4]:  3,        4,         4,        5
    L2node_NumChldrn[16]:3,2,3,  3,3,2,2,   2,2,2,2,  3,2,2,2,2  
 ****************************************************************/

struct SphereTree4L {    /* Number of levels is 4: root is at L0        */

   int dim,       /* dimension of each datum                     */
       NumData,   /* total num of data items                     */
       NumLevels, /* num of levels of the tree                   */
       N1,        /* N1: num of L1 clusters. N0=1 and is omitted */
       N2;        /* N2: total num of L2 clusters.               */
     //N3;        /* N3: total num of L3(i.e. leaf) clusters.    */


   int *L1node_NumChldrn, /* L1…[N1]: num of children of each L1 node.
                             Total number of L1 children = N2 =
			           L1NumChldrn[0]+L1NumChldrn[1]+...+L1NumChldrn[N1-1] */
       *L2node_NumChldrn, /* L2…[N2]: num of children of each L2 node. 
                             Total number of L2 children = N3 = 
                             L2NumChldrn[0]+L2NumChldrn[1]+...+L2NumChldrn[N2-1] */
       //*L3node_NumChldrn, /* L3…[N3]=L3NodeSize: num of children 
                            // (L3 nodes are leaves) on each L3 node               */
       *L1node_StartChild,/* L1…[N1]: L2 cluster index of the 
                             start child of each L1 node                         */
       *L2node_StartChild,/* L2…[N2]: L3 cluster index of the 
                             start child of each L2 node                         */
       //*L3node_StartChild,/* L3node_StartChild[N3]=L3node_StartDatum: 
                             //index of the start datum of each L3 node            */
       *L1node_DataSize,   /* L1…[N1]: num of data of each L1 cluster */
       *L2node_DataSize,   /* L2…[N2]: num of data of each L2 cluster */
       *L3node_DataSize,   /* L3…[N3]: num of data of each L3 cluster */
       *L1node_StartDatum, /* L1…[N1]: start datum of each L1 cluster */
       *L2node_StartDatum; /* L2…[N2]: start datum of each L2 cluster */
       //*L3node_StartDatum; /* L3…[N3]: start datum of each L3 cluster */

   double *L1CentersRadii, /* L1centerRadii[N1*(DIM+1)]               */ 	
          *L2CentersRadii, /* L2centerRadii[N2*(DIM+1)]               */
          //*L3CentersRadii, /* L3centerRadii[N3*(DIM+1)]               */	
          *data;
}; /**********************  end of struct LPT4L  **********************/


int LPtree_construc(int dim, int ndata, double *data, struct SphereTree4L *ptr_root) ;
int SStree_construc(int dim, int ndata, double *data, struct SphereTree4L *ptr_root) ;


int allneighbors_search(int dim, int ndata, double *data, struct SphereTree4L *ptr_root, 
                        double *query, double delta, int *num_outdata, int *outdata); 
/* Returns the number of points to which the distances from the query ere calculated */





#endif 


