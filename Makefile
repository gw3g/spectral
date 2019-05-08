CC	= g++
CFLAGS	= -lm -lgsl -lgslcblas
SDIR	= src
ODIR	= build
OUT	= out
TARGET	= bin/rho

# need to make a choice of ONE suffix
SRCEXT	= cpp

# find ALL *.SRCEXT files in ~/src directory
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

clean:
	$(RM) -r $(ODIR) $(TARGET); $(RM) bin/*;

dileptons: $(OBJ)
	$(CC) $(CFLAGS) $^ $(INC) spike/dileptons.cpp -o bin/dileptons

.PHONY: clean

