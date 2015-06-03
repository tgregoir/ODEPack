all: demo

demo: test.c
	@gcc -c test.c
	@gcc test.o ../cfiles/*.o -o demo -lm

clean:
	@rm -f test.o demo
