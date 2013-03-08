CXX=icpc
CC=icc

# use -lipo instead of -ipo for gcc

CXXFLAGS = -openmp -Wno-deprecated -Wcheck -Wall -fno-common -std=c++0x -fno-gnu-keywords -parallel -Ofast -xavx -no-prec-div -vec-report=1 -no-prec-sqrt -i-static -fno-alias -fno-exceptions -ipo -I/usr/local/include/vtk-5.10

OBJS = ACA2-2013.o Mesh.o Smooth.o SVD2x2.o

LIBS = -L/usr/lib64/vtk -L/usr/local/lib/vtk-5.10 -lvtkIO -lvtkFiltering -lvtkCommon -lvtkzlib -lvtkexpat -lvtksys -ldl -lpthread

# File names
EXEC = ACA2-2013
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Main target
$(EXEC): $(OBJECTS)
	$(CXX) $(OBJS) $(CXXFLAGS) $(LIBS) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(LIBS) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJS)
