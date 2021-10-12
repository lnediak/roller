all: test

test: main.o
	g++ -o $@ $+ -ldl -lGL -lglfw -lX11 -lXxf86vm -lXrandr -lpthread -lXi -lXinerama -lXcursor

main.o: main.cpp *.hpp ../*.hpp
	g++ --std=c++11 -Wall -Wextra -pedantic -g3 -O0 -I.. -c -o $@ $<

clean:
	/bin/rm -f test main.o

