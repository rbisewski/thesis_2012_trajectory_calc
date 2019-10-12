#
# Makefile for the project
#
PROJECT_NAME=thesis_2012_trajectory_calc

# Location of the source code
SRC_DIR=.

# Includes
INCLUDES= -I/usr/include

# Standard Libraries
LIBS= -L/usr/lib

# Flags
CFLAGS = -std=c++17 \
         -O2 \
         -fpic \
         -Wall \
         -Wextra \
         -Wpedantic \
         -Wno-missing-braces

# Compiler
CC = g++

# C sources
SRC = $(SRC_DIR)/trajectory.cpp

# C objects
OBJ = ${SRC:.c=.o}


#
# Makefile options
#


# State the "phony" targets
.PHONY: all clean


all: clean ${PROJECT_NAME}

.c.o:
	@echo CC $<
	@${CC} -c ${CFLAGS} $< -o ${<:.c=.o}

${PROJECT_NAME}: ${OBJ}
	@mkdir -p bin
	@echo ${CC} ${CFLAGS} ${INCLUDES} -o ./bin/$@ ${LIBS} $^
	@${CC} ${CFLAGS} ${INCLUDES} -o ./bin/$@ ${LIBS} $^

clean:
	@echo "Cleaning away old build..."
	@rm -fr bin $(SRC_DIR)*.o