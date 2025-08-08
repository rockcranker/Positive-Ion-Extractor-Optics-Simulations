CC = g++
LIBS = `pkg-config --libs ibsimu-1.0.6`
LDFLAGS = -Wall -g
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

analysis: analysis.o
	$(CC) $(LDFLAGS) $^ $(LIBS) -o $@

.cpp.o:
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) *.o  analysis *~
