GCC = gcc
TARGET = src/main.c
EXE = solver.out

INCLUDE_PATH = -I./include
LINKS = -lgsl-lc

MODULES = src/fe_section.c
ANC_MODULES = src/shape_functions.c

DEBUG = 1
DEBUG_FLAGS = -fsanitize=leak \
			  -fsanitize=address \
			  -g -O1

all: $(EXE)

$(EXE): $(MODULES:.c=.o) $(TARGET:.c=.o) 
	$(GCC) $(DEBUG_FLAGS) -o $@ $^

$(MODULES:.c=.o): $(MODULES) $(ANC_MODULES)
	$(GCC) $(INCLUDE_PATH) -c $< -o $@ $(LINKS)

$(TARGET:.c=.o): $(TARGET)
	$(GCC) $(INCLUDE_PATH) -c $< -o $@

clean:
	rm $(EXE) ./src/*.o







