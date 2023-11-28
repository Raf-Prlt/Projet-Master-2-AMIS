run: main
	./main

main: main.o
	gcc -o main main.o 

main.o: main.c
	gcc -c main.c 

clean:
	rm -f *.o main

valgrind: main
	valgrind ./main

