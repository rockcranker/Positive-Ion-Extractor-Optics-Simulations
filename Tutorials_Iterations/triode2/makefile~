CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

plasmacyl: triode.o
	$(CC) -o triode triode.o $(LDFLAGS)

plasmacyl.o: plasmacyl.cpp
	$(CC) -c -o triode.o triode.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o triode
