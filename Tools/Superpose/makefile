THEFILE = superpose
SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h*)
OBJECTS = $(SOURCES:.cpp=.o)

all: $(THEFILE)

$(THEFILE): $(OBJECTS)
	gcc -L/usr/local/Cellar/boost/1.57.0/lib/ -lboost_program_options-mt -lboost_filesystem-mt -lboost_thread-mt -O $(shell root-config --libs) -lTreePlayer     $(OBJECTS) -o $(THEFILE)
%.o: %.cpp $(HEADERS)
	gcc -I/usr/local/Cellar/boost/1.57.0/include $(shell root-config --cflags) -c $< -o $@

clean:
	rm -f $(THEFILE) $(OBJECTS)