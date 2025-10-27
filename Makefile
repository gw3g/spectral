CC	= g++
CFLAGS	= -O3 -std=c++11 -lm -lgsl -lgslcblas
CFLAGS += -fopenmp
TARGET	= bin/rho
SRCEXT	= cpp
SDIR	= src
ODIR	= build
OUT	= out

SRC	= $(wildcard $(SDIR)/*.$(SRCEXT))
OBJ	= $(patsubst $(SDIR)/%,$(ODIR)/%,$(SRC:.$(SRCEXT)=.o))
INC	= -I inc

$(TARGET): $(OBJ)
	@mkdir -p bin
	@mkdir -p $(OUT)/data
	$(CC) $^ $(INC) main.cpp -o $(TARGET) $(CFLAGS)

$(ODIR)/%.o: $(SDIR)/%.$(SRCEXT)
	@mkdir -p $(ODIR)
	$(CC) $(INC) -c -o $@ $< $(CFLAGS)

dileptons: $(OBJ)
	$(CC) $^ $(INC) examples/dileptons.cpp -o bin/dil $(CFLAGS)

neutrinos: $(OBJ)
	$(CC) $^ $(INC) examples/neutrinos.cpp -o bin/nu $(CFLAGS)

average: $(OBJ)
	$(CC) $^ $(INC) examples/average.cpp -o bin/ave $(CFLAGS)

integrand: $(OBJ)
	$(CC) $^ $(INC) examples/integrand.cpp -o bin/int $(CFLAGS)

relations: $(OBJ)
	$(CC) $^ $(INC) examples/relations.cpp -o bin/rel $(CFLAGS)

clean:
	$(RM) -r $(ODIR) $(TARGET); $(RM) bin/*;

.PHONY: clean

