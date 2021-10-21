# 
# ------ MAKEFILE FILE Unified ----------------
# 

CC     = gcc 
CXX = g++


# CXXFLAGS = -w -Wall -Werror -g -O3 -static -Dusing_linux 
CXXFLAGS = -w  -g -O3 -static -Dusing_linux 

PGM    = ../bin/unified

LIB    = $(PGM).a

.PRECIOUS :   $(LIB) $(PGM)

ELEMENTS =  \
 $(LIB)(main.o) 

$(PGM) : $(LIB) 
	 $(CXX) $(LIB) \
	        -g  \
	        -lm \
	        -o $(PGM)

$(LIB) :  $(ELEMENTS)
