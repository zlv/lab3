SOURCEOBJS = PolStr.o eigenvector.cpp
main :
	g++ -o lab3 $(SOURCEOBJS)
main-debug :
	g++ -g -O1 -o lab3 $(SOURCEOBJS)
