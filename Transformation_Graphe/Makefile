run: main
	./main

main: main.o
	gcc -o main main.o 

main.o: main.c
	gcc -c main.c 

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
	valgrind ./main --track-origins=yes

