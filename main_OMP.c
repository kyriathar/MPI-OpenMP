#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>

#include "grid.h"
#include "filter.h"
#include "reduce.h"
#include "neighbors.h"
#include "read_write.h"
#include "my_functions.h"
#include "send_receive.h"

int  coord[2];

int main(int argc, char** argv) {

    /*VARS*/
    int i,j;
    int my_rank =0;
    int number_of_threads = 2;
    int processes =0;
    int inner_i =0;
    int inner_j =0;
    int array_elements=0;
    int check_value; 
    int return_value = 1;
    unsigned char * temp = NULL ;
    unsigned char * local_array_in = NULL ;
    unsigned char * local_array_out = NULL ;
    int grid_dimensions ,dims[2], period[2], reorder;
    MPI_Comm new_comm;
    int Neighbors[8];
    MPI_Datatype row,column;
    MPI_Request send[8],receive[8];
    Filter filter ;
    
    double time_in;
    double time_out;
    
    /*FILTER*/
    init_filter(&filter);
    create_s_filter(&filter,1);
    value_filter(&filter);
    
    /*FUNCTION PARAMS*/
     int x=0 ;
     int y=0 ;
     int fx=0 ,fy=0 ;
     int sum=0 ;
     int * stop = NULL ;
    
    
    MPI_Init(NULL, NULL);
    //MPI_Init_thread( 0, 0, MPI_THREAD_FUNNELED, NULL );
    MPI_Comm_size(MPI_COMM_WORLD,&processes);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank); 
    
    
    /*LOCAL ARRAY*/
    inner_i=    IMAGE_I / sqrt(processes);
    inner_j =   IMAGE_J / sqrt(processes);
    array_elements = (inner_i+2)*(inner_j+2) ;
    local_array_in = (unsigned char *)my_malloc(array_elements*sizeof(unsigned char)); 
    memset(local_array_in,255,array_elements*sizeof(unsigned char));
    
    local_array_out = (unsigned char *)my_malloc(array_elements*sizeof(unsigned char)); 
    memset(local_array_out,0,array_elements*sizeof(unsigned char));
    
    /*GRID*/
    set_grid(&grid_dimensions,dims ,processes,period,&reorder);
    MPI_Cart_create(MPI_COMM_WORLD, grid_dimensions,dims, period, reorder, &new_comm);
    
    /*FIND NEIGHBORS AND MY COORDS*/
    for(i=0;i<8;i++)
        find_my_neighbor(my_rank,i ,2,dims,new_comm,(Neighbors+i));
    
    /*PARALLEL READ*/
    read_all("../Photo/waterfall_grey_1920_2520.raw",new_comm,inner_i,inner_j,local_array_in); 
    
    /*DATATYPES*/
    MPI_Type_contiguous( inner_j, MPI_UNSIGNED_CHAR, &row );
    MPI_Type_commit(&row);
    
    MPI_Type_vector( inner_i,1,inner_j+2,MPI_UNSIGNED_CHAR, &column);
    MPI_Type_commit(&column);
    
    
    /*SEND RECEIVE + FILTER*/


stop = (int *)malloc(number_of_threads*sizeof(int));
   for(j=0;j<number_of_threads;j++)
        stop[j] = 0 ; 

time_in = MPI_Wtime();

omp_set_num_threads(number_of_threads);
#pragma omp parallel shared(stop,return_value) private(i)
{    
    
   
    
    for(i=0;(i<1000) && return_value ;i++)
    {
        //printf("i = %d\n",i);
        //sleep(2);
        //MASTER
#pragma omp master
{        
        i_send_all(local_array_in,inner_i,inner_j,row,column,Neighbors,new_comm,send);
        i_receive_all(local_array_in,inner_i,inner_j,row,column,Neighbors,new_comm,receive); 
}
#pragma omp barrier
        
        //filterise(&filter,inner_i,inner_j,local_array_in,local_array_out,INNER); // NA MPEI EDW O KWDIKAS
        //INNER
        x=0 ;
        y=0 ;
        fx=0 ,fy=0 ;
        sum=0 ;

#pragma omp for private(y,fx,fy,sum)
        for(x=2;x<inner_i;x++)
        {
            for(y=2; y<inner_j ;y++)          //y < 1920
            {
                sum = 0 ;
                for(fx=-1;fx<=1;fx++)
                    for(fy=-1;fy<=1;fy++)
                        sum += local_array_in[(x-fx)*(inner_j+2)+(y-fy)]*filter.array[fx+1][fy+1];
                local_array_out[x*(inner_j+2) + y] = (int)sum/filter.sum ;
            }
        }
        //MASTER
#pragma omp master
{
        wait_receive(receive);
}
#pragma omp barrier
        //filterise(&filter,inner_i,inner_j,local_array_in,local_array_out,OUTER); // NA MPEI EDW O KWDIKAS ? ?
        //OUTER        
        /*first row*/
        y = 1 ;
#pragma omp for private(fx,fy,sum)
        for(x=1;x<inner_i+1;x++)
        {
            sum = 0 ;
                for(fx=-1;fx<=1;fx++)
                    for(fy=-1;fy<=1;fy++)
                        sum += local_array_in[(x-fx)*(inner_j+2)+(y-fy)]*filter.array[fx+1][fy+1];
                local_array_out[x*(inner_j+2) + y] = (int)sum/filter.sum ;
        }
        /*last row*/
        y = inner_j ;
#pragma omp for private(fx,fy,sum)
        for(x=1;x<inner_i+1;x++)
        {
            sum = 0 ;
                for(fx=-1;fx<=1;fx++)
                    for(fy=-1;fy<=1;fy++)
                        sum += local_array_in[(x-fx)*(inner_j+2)+(y-fy)]*filter.array[fx+1][fy+1];
                local_array_out[x*(inner_j+2) + y] = (int)sum/filter.sum ;
        }
        /*first column - 2 corners */
        x = 1 ;
#pragma omp for private(fx,fy,sum)
        for(y=2;y<inner_j;y++)
        {
            sum = 0 ;
                for(fx=-1;fx<=1;fx++)
                    for(fy=-1;fy<=1;fy++)
                        sum += local_array_in[(x-fx)*(inner_j+2)+(y-fy)]*filter.array[fx+1][fy+1];
                local_array_out[x*(inner_j+2) + y] = (int)sum/filter.sum ;
        }
        /*last column - 2 corners */
        x = inner_i ;
#pragma omp for private(fx,fy,sum)
        for(y=2;y<inner_j;y++)
        {
            sum = 0 ;
                for(fx=-1;fx<=1;fx++)
                    for(fy=-1;fy<=1;fy++)
                        sum += local_array_in[(x-fx)*(inner_j+2)+(y-fy)]*filter.array[fx+1][fy+1];
                local_array_out[x*(inner_j+2) + y] = (int)sum/filter.sum ;
        }
#pragma omp master
{            
        temp = local_array_in;
        local_array_in = local_array_out;
        local_array_out = temp;

        
}
#pragma omp barrier
        //check_value=check_convergence(inner_i,inner_j,local_array_in,local_array_out);
        //int x,y ;

        
#pragma omp for private(y) 
        for(x=1; x<inner_i+1 ;x++)
        {
            //printf("thread no -> %d\n",omp_get_thread_num());
            if(stop[omp_get_thread_num()] == 0)
            {
                for(y=1; y<inner_j+1 ;y++)          //y < 1920
                {
                    if(local_array_in[x*(inner_j+2) + y] != local_array_out[x*(inner_j+2) + y])
                    {
    //#pragma omp critical
                        stop[omp_get_thread_num()] = 1;
                    }
                }
            }
        }
        //printf("KOLAW EDW\n");
        
#pragma omp barrier
        
#pragma omp master
{
       check_value = 0 ;     
       for(j=0;j<number_of_threads;j++)
       {
           //printf("stop[%d] = %d\n",j,stop[j]);
           
           if(stop[j]== 1)
           {
               check_value = 1 ;
               stop[j] = 0;
               
           }
        }
       
        //return 0 ;     
            
        MPI_Allreduce( &check_value,&return_value, 1, MPI_INT, MPI_SUM, new_comm );  
        
        if(return_value == 0)
            printf("Epanalipseis  = %d \n",i);
        //printf("i = %d\n",i);
        
}
#pragma omp barrier
        
        
  
        //MASTER
#pragma omp master
{
        wait_send(send);
}
#pragma omp barrier
    }
}//#pragma omp parallel
    
        
        
    time_out = MPI_Wtime();
    
    if(my_rank ==0)
    {
        printf("Time taken: %1.5f seconds\n", time_out-time_in);
    }
    
    /*PARALLEL WRITE*/
    
    compose_pic("../Photo/out.raw",new_comm,inner_i,inner_j,local_array_in);    
    
    
    MPI_Finalize();
    
    /*FREE FOR ALL*/
    free(local_array_in);
    free(local_array_out);
    free(stop);
    
    return (EXIT_SUCCESS);
}