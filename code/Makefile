MPICXX = mpic++
CXXFLAGS = -Wall

SRC = main.cpp custom_collectives.h custom_collectives.cpp
EXEC = primitives

all: $(EXEC)

$(EXEC): $(SRC)
	$(MPICXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(EXEC)