CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

triode2: triode2.o
	$(CC) -o triode2 triode2.o $(LDFLAGS)

triode2.o: triode2.cpp
	$(CC) -c -o triode2.o triode2.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o triode2
