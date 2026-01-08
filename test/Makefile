INCLUDE_PATH = -I../include
CC = gcc

MODULES = ../src/fe_section.c 
SUBMODULES = ../src/shape_functions.c \
			 ../src/composition_functions.c
MODULES_BASE = $(notdir $(MODULES))

LIBS = -lgsl -lcheck -lsubunit -lrt -lm -pthread
DEBUG_FLAGS = -fsanitize=address,undefined \
			  -fsanitize=leak \
			  -g -O0

U_SOLVER = unit/test_solver_comp.c
U_ELEMENT = unit/test_element_creation.c
I_PARSER = integration/test_parser.c
I_SOLVER = integration/test_solver.c

EXE_ASSEMBLY = test_comp_and_assembly.out
EXE_ELEMENT = test_element.out
EXE_PARSER = test_parser.out
EXE_SOLVER = test_solver.out

all: $(EXE_ASSEMBLY) $(EXE_ELEMENT) $(EXE_PARSER) $(EXE_SOLVER)

# Compile modules first into this directory
./unit/$(MODULES_BASE:.c=.o): $(MODULES) $(SUBMODULES)
	$(CC) $(INCLUDE_PATH) $(DEBUG_FLAGS) -c $< -o $@

./integration/$(MODULES_BASE:.c=.o): $(MODULES) $(SUBMODULES)
	$(CC) $(INCLUDE_PATH) $(DEBUG_FLAGS) -c $< -o $@

# Object files for the test executables
$(U_SOLVER:.c=.o): $(U_SOLVER)
	$(CC) $(INCLUDE_PATH) $(DEBUG_FLAGS) -c $< -o $@

$(U_ELEMENT:.c=.o): $(U_ELEMENT)
	$(CC) $(INCLUDE_PATH) $(DEBUG_FLAGS) -c $< -o $@

$(I_PARSER:.c=.o): $(I_PARSER)
	$(CC) $(INCLUDE_PATH) $(DEBUG_FLAGS) -c $< -o $@

$(I_SOLVER:.c=.o): $(I_SOLVER)
	$(CC) $(INCLUDE_PATH) $(DEBUG_FLAGS) -c $< -o $@

# Solver Exectuable Instructions
$(EXE_ASSEMBLY): $(U_SOLVER:.c=.o) ./unit/$(MODULES_BASE:.c=.o)
	$(CC) $(INCLUDE_PATH) $^ $(LIBS) $(DEBUG_FLAGS) -o $@

$(EXE_ELEMENT): $(U_ELEMENT:.c=.o) ./unit/$(MODULES_BASE:.c=.o)
	$(CC) $(INCLUDE_PATH) $^ $(LIBS) $(DEBUG_FLAGS) -o $@

$(EXE_PARSER): $(I_PARSER:.c=.o) ./integration/$(MODULES_BASE:.c=.o)
	$(CC) $(INCLUDE_PATH) $^ $(LIBS) $(DEBUG_FLAGS) -o $@

$(EXE_SOLVER): $(I_SOLVER:.c=.o) ./integration/$(MODULES_BASE:.c=.o)
	$(CC) $(INCLUDE_PATH) $^ $(LIBS) $(DEBUG_FLAGS) -o $@

clean:
	rm ./unit/*.o ./integration/*.o *.out


