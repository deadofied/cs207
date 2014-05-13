#
# 'make'        build executable file
# 'make clean'  removes all .o and executable files
#

# Executables to build
EXEC    += primes
SDLEXEC += viewer
EXEC += test_nodes
EXEC += test_edges
SDLEXEC += subgraph
SDLEXEC += shortest_path
EXEC    += mtl_test
SDLEXEC += mass_spring
SDLEXEC += poisson

# Get the shell name to determine the OS
UNAME := $(shell uname)

# Define the C++ compiler to use
CXX := $(shell which g++) -std=gnu++0x
#CXX = clang++ -std=c++0x
ifeq ($(UNAME), Darwin)
CC := $(shell which gcc) -O3
SDLOBJS := CS207/SDLMain.o
endif

# Dependency directory and flags
DEPSDIR := $(shell mkdir -p .deps; echo .deps)
# MD: Dependency as side-effect of compilation
# MF: File for output
# MP: Include phony targets
DEPSFILE = $(DEPSDIR)/$(notdir $*.d)
DEPSFLAGS = -MD -MF $(DEPSFILE) #-MP

# Define any directories containing header files
#   To include directories use -Ipath/to/files
INCLUDES += -I.
INCLUDES += -IMTL-4.0.9507-Linux/usr/include

# Define CXX compile flags
CXXFLAGS += -fopenmp -funroll-loops -O3 -W -Wall -Wextra -Wfatal-errors

# Define any directories containing libraries
#   To include directories use -Lpath/to/files
LDFLAGS +=

# Define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX
ifeq ($(UNAME), Linux)
LDLIBS += -lSDL -lGL -lGLU
endif
ifeq ($(UNAME), Darwin)
LDLIBS += -framework SDL -framework OpenGL -framework Cocoa
endif

##################
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
##################

# 'make' - default rule
all: $(EXEC) $(SDLEXEC)

# Default rule for creating an exec of $(EXEC) from a .o file
$(EXEC): % : %.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Default rule for creating an exec of $(EXEC) from a .o file
$(SDLEXEC): % : %.o $(SDLOBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Default rule for creating a .o file from a .cpp file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEPSFLAGS) -c -o $@ $<

# Extra dependencies for executables
#   Nothing here

# 'make clean' - deletes all .o files, exec, and dependency files
clean:
	-$(RM) *.o $(EXEC) $(SDLEXEC) $(SDLOBJS)
	$(RM) -r $(DEPSDIR)

# Define rules that do not actually generate the corresponding file
.PHONY: clean all

# Include the dependency files
-include $(wildcard $(DEPSDIR)/*.d)
