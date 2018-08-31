INC=-I/home/r900-1/milos/hdf5-1.6.7/include
LIB=-L/home/r900-1/milos/hdf5-1.6.7/lib -lhdf5 -lm

extract_particles: extract_particles.c
	gcc -g $(LIB) $(INC) -o extract_particles extract_particles.c

clean:
	rm extract_particles

