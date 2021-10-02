all: roller

roller: main.o
	g++ -o $@ $+ -ldl -lGL -lglfw -lX11 -lXxf86vm -lXrandr -lpthread -lXi -lXinerama -lXcursor

main.o: main.cpp *.hpp
	g++ --std=c++11 -Wall -Wextra -pedantic -g3 -O2 -c -o $@ $<

clean:
	/bin/rm -f roller main.o

