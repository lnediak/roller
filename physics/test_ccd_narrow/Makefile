all: test

test: main.o
	g++ -o $@ $+

main.o: main.cpp ../*.hpp ../../*.hpp
	g++ --std=c++11 -Wall -Wextra -pedantic -g3 -O2 -I.. -I../.. -c -o $@ $<

clean:
	/bin/rm -f test main.o

