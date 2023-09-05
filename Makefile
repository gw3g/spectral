CC	= g++
CFLAGS	= -std=c++11 -lm -lgsl -lgslcblas
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
	$(CC) $(CFLAGS) $^ $(INC) main.cpp -o $(TARGET)

$(ODIR)/%.o: $(SDIR)/%.$(SRCEXT)
	@mkdir -p $(ODIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

dileptons: $(OBJ)
	$(CC) $(CFLAGS) $^ $(INC) examples/dileptons.cpp -o bin/dil

neutrinos: $(OBJ)
	$(CC) $(CFLAGS) $^ $(INC) examples/neutrinos.cpp -o bin/nu

average: $(OBJ)
	$(CC) $(CFLAGS) $^ $(INC) examples/average.cpp -o bin/ave

integrand: $(OBJ)
	$(CC) $(CFLAGS) $^ $(INC) examples/integrand.cpp -o bin/int

relations: $(OBJ)
	$(CC) $(CFLAGS) $^ $(INC) examples/relations.cpp -o bin/rel

clean:
	$(RM) -r $(ODIR) $(TARGET); $(RM) bin/*;

.PHONY: clean

