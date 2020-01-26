gpu : main.cpp
	nvcc -x cu main.cpp -o gm

cpu: main.cpp
	g++ -o m main.cpp

clean: gm m 
	rm -f gm m