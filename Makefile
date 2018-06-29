TARGET = g2
INCL = -I/opt/boost-1.57.0/include/
CXX = g++
CXXFLAGS = -O3 -std=c++11

default: $(TARGET)

all: default

test:
	./$(TARGET) -m 0.9 example/SIM.NE_20000.MATCH_FREQ.SHAPEIT.haps example/SIM.NE_20000.MATCH_FREQ.SHAPEIT.sample example/genMap.1KG.b37.chr1.map example/SIM.NE_20000.MATCH_FREQ.INFERRED.match
	diff -s -q example/SIM.NE_20000.MATCH_FREQ.INFERRED.match example/SIM.NE_20000.MATCH_FREQ.INFERRED.match.TESTED
	bash example/accuracy.sh example/SIM.NE_20000.TRUE.match example/SIM.NE_20000.MATCH_FREQ.INFERRED.match example/SIM.NE_20000.MATCH_FREQ.INFERRED.accuracy
clean:
	rm -f $(TARGET) example/SIM.NE_20000.MATCH_FREQ.INFERRED.match
