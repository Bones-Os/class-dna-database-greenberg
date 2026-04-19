CXX = g++
CXXFLAGS = -Wall

proj4: dnadb.o driver.o
	$(CXX) $(CXXFLAGS) dnadb.o driver.o -o proj4

dnadb.o: dnadb.cpp dnadb.h
	$(CXX) $(CXXFLAGS) -c dnadb.cpp

driver.o: driver.cpp
	$(CXX) $(CXXFLAGS) -c driver.cpp

mytest.o: dnadb.o driver.o mytest.cpp
	$(CXX) $(CXXFLAGS) -c mytest.cpp

clean:
	rm *.o*
	rm *~

run:	./proj4

test1: mytest.o
	$(CXX) $(CXXFLAGS) -g mytest.o dnadb.o -o mytest
	./mytest

test2: mytest.o
	$(CXX) $(CXXFLAGS) -g mytest.o dnadb.o -o mytest #dnadb.cpp mytest.o -o mytest
	valgrind -s --track-origins=yes ./mytest

debug: mytest.o
	$(CXX) $(CXXFLAGS) dnadb.cpp mytest.o -o mytest
	gdb dnadb.h dnadb.cpp mytest.cpp --args mytest
