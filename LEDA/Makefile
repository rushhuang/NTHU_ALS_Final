LEDAROOT=./LEDA_lib/LEDA

.PHONY: all clean

all:	test

test:	test.cpp
	$(CXX) -o $@ $< -I$(LEDAROOT)/incl -L$(LEDAROOT) -lG -lL -lm

clean:
	$(RM) test

