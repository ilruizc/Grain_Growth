All: 1.x 2.x
1: Paper1.cpp
	g++ -g -fsanitize=address  -o 1.x Paper1.cpp
	./1.x
2: Paper2.cpp 
	g++  -g -fsanitize=address $<  -o 2.x 
	time ./2.x
clean: 
	rm ./a.out *~ *.dat *.x
