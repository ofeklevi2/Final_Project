CC = gcc 
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm

spkmeans: spkmeans.o 
	$(CC) -o spkmeans spkmeans.o $(CFLAGS)

spkmeans.o: spkmeans.c 
	$(CC) -c spkmeans.c $(CFLAGS)

