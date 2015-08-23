.PHONY: clean

CC=cc
CFLAGS=-g3 -std=c99
AR=ar
ARFLAGS=rcs

all: libtinymcmc.a

libtinymcmc.a: objects/*.o
	$(AR) $(ARFLAGS) libtinymcmc.a objects/*.o

objects/%.o: src/%.c
		$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f tinymcmc.o libtinymcmc.a
