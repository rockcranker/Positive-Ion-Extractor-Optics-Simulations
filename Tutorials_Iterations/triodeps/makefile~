CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

triode3d: triode3d.o
	$(CC) -o triode3d triode3d.o $(LDFLAGS)

triode3d.o: triode3d.cpp
	$(CC) -c -o triode3d.o triode3d.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o triode3d
