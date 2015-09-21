CXX=g++
#HOME=/zhenyuhan/
#FASTJETLOCATION=$(HOME)/fastjet
#FASTJETLIB=`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --libs --plugins`
FASTJETLIB=`fastjet-config --cxxflags --libs --plugins`
all: example
example: double_cone.hh example.cc
	$(CXX) $(FASTJETLIB) -o example example.cc
