CC = g++
LIBS = `pkg-config --libs ibsimu-1.0.6`
LIBS += -L/home/imorozov/src/root/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lROOTNTuple -lMultiProc -lROOTDataFrame -lROOTNTupleUtil -Wl,-rpath,/home/imorozov/src/root/lib -pthread -lm -ldl -rdynamic
LDFLAGS = -Wall -g
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`
CXXFLAGS += `root-config --cflags`
LDFLAGS += `root-config --ldflags`



BES_sim_v2: BES_sim_v2.o
	$(CC) $(LDFLAGS) $^ $(LIBS)  -o $@ 

.cpp.o:
	$(CC) $(CXXFLAGS) -c $< -o $@ 

clean:
	$(RM) *.o  BES_sim_v2 *~
