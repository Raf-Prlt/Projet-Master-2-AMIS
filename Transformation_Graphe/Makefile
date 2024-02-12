SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)
CFLAGS = -W -Wall -Wextra -pedantic

run: main
	./main

main: $(OBJ) 
	gcc -o $@ -Lnauty2_8_8/ $^  nauty2_8_8/nauty.a 

# main.o: main.c
# 	gcc -c  main.c $(CFLAGS)

%.o: %.c 
	gcc -g -o $@ -c $< $(CFLAGS)

clean:
	rm -f *.o main
	rm -f *.o test

runtest: test
	./test

test: test.o
	gcc -o test test.o 

test.o: test.c
	gcc -c test.c 

valgrind: main
	valgrind --track-origins=yes ./main

debug: $(OBJ) 
	gcc -g -o $@ -Lnauty2_8_8/ $^  nauty2_8_8/nauty.a 
	gdb ./debug

