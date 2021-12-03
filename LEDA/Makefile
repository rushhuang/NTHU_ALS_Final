LEDAROOT=./LEDA_lib/LEDA

.PHONY: all clean

all:	test map

test:	test.cpp
	# $(CXX) -o $@ $< -I$(LEDAROOT)/incl -L$(LEDAROOT) -lG -lL -lm
	$(CXX) -o $@ $< -I$(LEDAROOT)/incl -L$(LEDAROOT) -lG -lL -lm -no-pie

map:	map.cpp
	# $(CXX) -o $@ $< -I$(LEDAROOT)/incl -L$(LEDAROOT) -lG -lL -lm
	$(CXX) -o $@ $< -I$(LEDAROOT)/incl -L$(LEDAROOT) -lG -lL -lm -no-pie

clean:
	$(RM) test

