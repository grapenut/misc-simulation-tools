INC=-I/home/r900-1/milos/hdf5-1.6.7/include
LIB=-L/home/r900-1/milos/hdf5-1.6.7/lib -lhdf5 -lm

flash2gadget: flash2gadget.c
	gcc -g $(LIB) $(INC) -o flash2gadget flash2gadget.c

flash2text: flash2text.c match_tags.c
	gcc -g $(LIB) $(INC) -o flash2text flash2text.c
	gcc -g $(LIB) $(INC) -o tag match_tags.c

gadget2flash: gadget2flash.c
	gcc -g $(LIB) $(INC) -o gadget2flash gadget2flash.c

clean:
	rm flash2gadget flash2text tag gadget2flash


