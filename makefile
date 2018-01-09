# PATH
EXECPATH = ./bin/

## CLEANBLASTP
# source files
SOURCESCB = ./src/cleanblastp/*.cpp

# executable name and path def
EXECNAMECB = cleanblastp

## COMPOSITESEARCH
# igraph library
PATH_IGRAPH_LIB = -L/usr/local/lib
PATH_IGRAPH_INC = -I/usr/local/include/igraph

# source files
SOURCESCS = ./src/compositesearch/*.cpp

# executable name and path def
EXECNAMECS = compositeSearch

# parameters for g++ compiler
GCC_NAME    = g++
GCC_FLAGS   = -std=c++11
GCC_LIBS    = -ligraph -pthread

.SUFFIXES: .cpp

exec: $(OBJECTS)
	$(GCC_NAME) $(SOURCESCS) $(PATH_IGRAPH_INC) $(PATH_IGRAPH_LIB) $(GCC_LIBS) $(GCC_FLAGS) -o $(EXECPATH)$(EXECNAMECS)
	$(GCC_NAME) $(SOURCESCB) $(GCC_FLAGS) -o $(EXECPATH)$(EXECNAMECB)
clean:
	rm $(EXECNAMECS)
	rm $(EXECNAMECB)
# END
