CC	= g++
CFLAGS	= -lm -lgsl -lgslcblas
SDIR	= ../src
ODIR	= ../build

SRCEXT	= cpp

SRC	= $(wildcard $(SDIR)/*.$(SRCEXT))
OBJ	= $(patsubst $(SDIR)/%,$(ODIR)/%,$(SRC:.$(SRCEXT)=.o))
INC	= -I ../inc

$(ODIR)/%.o: $(SDIR)/%.$(SRCEXT)
	@mkdir -p $(ODIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

dileptons: $(OBJ)
	$(CC) $(CFLAGS) $^ $(INC) dileptons.cpp

average: $(OBJ)
	$(CC) $(CFLAGS) $^ $(INC) average.cpp


