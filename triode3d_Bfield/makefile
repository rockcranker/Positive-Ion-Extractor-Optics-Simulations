CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

triode3dB: triode3dB.o
	$(CC) -o triode3dB triode3dB.o $(LDFLAGS)

triode3d.o: triode3dB.cpp
	$(CC) -c -o triode3dB.o triode3dB.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o triode3dB
