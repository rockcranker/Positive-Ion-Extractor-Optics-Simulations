CC = g++
LIBS = `pkg-config --libs ibsimu-1.0.6`
LIBS += -L/usr/lib64/root -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
LDFLAGS = -Wall -g
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`
CXXFLAGS += `root-config --cflags`
LDFLAGS += `root-config --ldflags`

objects = BES_sim analysis
all: $(objects)

.cpp.o:
	$(CC) $(CXXFLAGS) -c $< -o $@ 

$(objects): %: %.o
	$(CC) $(LDFLAGS) $^ $(LIBS)  -o $@ 

clean:
	$(RM) *.o  $(objects) *~
