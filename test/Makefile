INCLUDE_PATH = -I../include
CC = gcc

MODULES = ../src/fe_section.c 
SUBMODULES = ../src/shape_functions.c \
			 ../src/composition_functions.c
MODULES_BASE = $(notdir $(MODULES))

LIBS = -lgsl -lcheck -lsubunit -lrt -lm -pthread
DEBUG_FLAGS = -fsanitize=leak \
			  -fsanitize=address \
			  -g -O1

U_SOLVER = unit/test_solver_comp.c
U_ELEMENT = unit/test_element_creation.c

EXE_SOLVER = test_comp_and_assembly.out
EXE_ELEMENT = test_element.out

all: $(EXE_SOLVER) $(EXE_ELEMENT)

# Compile modules first into this directory
./unit/$(MODULES_BASE:.c=.o): $(MODULES) $(SUBMODULES)
	$(CC) $(INCLUDE_PATH) $(DEBUG_FLAGS) -c $< -o $@

# Object files for the test executables
$(U_SOLVER:.c=.o): $(U_SOLVER)
	$(CC) $(INCLUDE_PATH) $(DEBUG_FLAGS) -c $< -o $@

$(U_ELEMENT:.c=.o): $(U_ELEMENT)
	$(CC) $(INCLUDE_PATH) $(DEBUG_FLAGS) -c $< -o $@

# Solver Exectuable Instructions
$(EXE_SOLVER): $(U_SOLVER:.c=.o) ./unit/$(MODULES_BASE:.c=.o)
	$(CC) $(INCLUDE_PATH) $^ $(LIBS) $(DEBUG_FLAGS) -o $@

$(EXE_ELEMENT): $(U_ELEMENT:.c=.o) ./unit/$(MODULES_BASE:.c=.o)
	$(CC) $(INCLUDE_PATH) $^ $(LIBS) $(DEBUG_FLAGS) -o $@

clean:
	rm ./unit/*.o *.out


