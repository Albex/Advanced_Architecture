CXX=icpc
CC=icc

# use -lipo instead of -ipo for gcc

CXXFLAGS = -g -O3 -complex-limited-range -Wall -Wno-deprecated -I/usr/local/include/vtk-5.10 -ipo -static -fno-alias -fno-exceptions -prof-use -opt-report-file=opt.txt

OBJS = ACA2-2013.o Mesh.o Smooth.o SVD2x2.o

LIBS = -L/usr/lib64/vtk -L/usr/local/lib/vtk-5.10 -lvtkIO -lvtkFiltering -lvtkCommon -lvtkzlib -lvtkexpat -lvtksys -ldl -lpthread

TARGET = ACA2-2013

$(TARGET):	$(OBJS)
	@$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	@$(TARGET)

clean:
	@rm -f $(OBJS) $(TARGET)
