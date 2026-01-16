CC = gcc
AR = ar
AR_FLAGS = rcs

BUILD_DIR = build

INCLUDE_PATH = -I./include
LIB_PATH = -L./$(BUILD_DIR)

LIB_OUTPUT = $(BUILD_DIR)/libode.a

SOURCE = src/fe_section.c \
		 src/function_field.c
BUILD_OBJ = $(SOURCE:src/%.c=./$(BUILD_DIR)/%.o)
AUX_SOURCE = src/composition_functions.c \
			 src/shape_functions.c

EXE = solver.out
EXE_SOURCE = src/main.c
EXE_OBJECT = $(BUILD_DIR)/$(notdir $(EXE_SOURCE:.c=.o))
LINK_FLAG = -lode -lgsl -lm

CC_FLAGS = -g -O0

.PHONY = lib clean_all main

# Rules for main exectuable build
main: $(EXE)

$(EXE): $(EXE_OBJECT) $(LIB_OUTPUT) | $(BUILD_DIR)
	$(CC) $(INCLUDE_PATH) $(LIB_PATH) -O2 -Wall -Wextra $< -o $@ $(LINK_FLAG)

$(EXE_OBJECT): $(EXE_SOURCE) | $(BUILD_DIR)
	$(CC) $(INCLUDE_PATH) $(CC_FLAGS) -c $< -o $@

# Rules for the static library build
lib: $(LIB_OUTPUT)

$(LIB_OUTPUT): $(BUILD_OBJ) | $(BUILD_DIR)
	$(AR) $(AR_FLAGS) $@ $^

$(BUILD_DIR)/%.o: src/%.c $(AUX_SOURCE) | $(BUILD_DIR)
	$(CC) $(INCLUDE_PATH) $(CC_FLAGS) -c $< -o $@

# Ensure that the build directory exists
$(BUILD_DIR):
	mkdir -p $@

clean:
	rm *.out ./$(BUILD_DIR)/*





