SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)
CFLAGS = -W -Wall -Wextra -pedantic

run: main
	./main parametres.txt

main: $(OBJ) 
	gcc -o $@ -Lnauty2_8_8/ $^  nauty2_8_8/nauty.a 

# main.o: main.c
# 	gcc -c  main.c $(CFLAGS)

%.o: %.c 
	gcc -g -o $@ -c $< $(CFLAGS)

clean:
	rm -f *.o main

valgrind: main
	valgrind --track-origins=yes  ./main
#--exit-on-first-error=yes

debug: $(OBJ) 
	gcc -g -o $@ -Lnauty2_8_8/ $^  nauty2_8_8/nauty.a 
	gdb ./debug

