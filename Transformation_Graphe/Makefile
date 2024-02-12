SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)
CFLAGS = -W -Wall -Wextra -pedantic

run: main
	./main

main: $(OBJ) 
	gcc -o $@ $^  

# main.o: main.c
#	gcc -c  main.c $(CFLAGS)

%.o: %.c 
	gcc -o $@ -c $< $(CFLAGS)

clean:
	rm -f *.o main

valgrind: main
	valgrind ./main