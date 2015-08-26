.PHONY: clean

CC=cc
CFLAGS=-g3 -std=c99 -Irandamu/include
GDB=gdb
AR=ar
ARFLAGS=rcs
OBJECTS=objects/test.o objects/tinymcmc.o objects/rnd.o

all: libtinymcmc.a 

test: libtinymcmc.a $(OBJECTS)
	$(CC) $(CFLAGS) objects/*.o -o tests/test
	cd tests; $(GDB) ./test

libtinymcmc.a: $(OBJECTS)
	$(AR) $(ARFLAGS) libtinymcmc.a objects/*.o

objects/rnd.o: randamu/src/rng.c
	$(CC) $(CFLAGS) -c randamu/src/rng.c -o objects/rnd.o

objects/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@	

clean:
	rm -f tinymcmc.o libtinymcmc.a
