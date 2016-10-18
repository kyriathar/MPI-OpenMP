CC = mpicc 
OBJS = filter_OMP.o grid.o main_OMP.o my_functions.o neighbors.o read_write.o reduce_OMP.o send_receive.o
CFLAGS = -c 
LFLAGS = -o 


all : $(OBJS)
		$(CC) $(OBJS) $(LFLAGS) mpiopenmp -lm -fopenmp

filter_OMP.o : filter_OMP.c
			$(CC) $(CFLAGS) filter_OMP.c -fopenmp

grid.o : grid.c
			$(CC) $(CFLAGS) grid.c -fopenmp

main_OMP.o : main_OMP.c
			$(CC) $(CFLAGS) main_OMP.c -fopenmp

my_functions.o : my_functions.c
					$(CC) $(CFLAGS) my_functions.c -fopenmp

neighbors.o : neighbors.c
					$(CC) $(CFLAGS) neighbors.c -fopenmp

read_write.o : read_write.c
					$(CC) $(CFLAGS) read_write.c -fopenmp

reduce_OMP.o : reduce_OMP.c
			$(CC) $(CFLAGS) reduce_OMP.c -fopenmp

send_receive.o : send_receive.c
					$(CC) $(CFLAGS) send_receive.c -fopenmp


clean :
		rm -rf *.o mpiopenmp ../Photo/out.raw
