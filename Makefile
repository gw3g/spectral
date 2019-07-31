CC	= g++
CFLAGS	= -lm -lgsl -lgslcblas
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

clean:
	$(RM) -r $(ODIR) $(TARGET); $(RM) bin/*;

.PHONY: clean

