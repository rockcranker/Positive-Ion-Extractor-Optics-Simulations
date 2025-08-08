CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

triodeps: triodeps.o
	$(CC) -o triodeps triodeps.o $(LDFLAGS)

triodeps.o: triodeps.cpp
	$(CC) -c -o triodeps.o triodeps.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o triodeps
