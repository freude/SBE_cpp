CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main_SBE.cpp gain.cpp polarization1.cpp int_matrix.cpp exchange.cpp Fermi_level.cpp E_field.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=gain.out

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
