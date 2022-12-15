CC := g++
CFLAGS := -Wall -std=c++11
OPTFLAGS = -O3

SRCDIR := src
BUILDDIR := build
OUTPUTSDIR := outputs
GRAPHSDIR := graphs
TARGET := bin/csswm


SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ $(OPTFLAGS) -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c $(OPTFLAGS) -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET) $(OUTPUTSDIR)/* $(GRAPHSDIR)/*"; $(RM) -r $(BUILDDIR) $(TARGET) $(OUTPUTSDIR)/* $(GRAPHSDIR)/*