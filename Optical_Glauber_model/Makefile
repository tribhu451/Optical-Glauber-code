
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

#EXTRA_FLAGS = -D SIMPLE # EoS p=e/3
EXTRA_FLAGS  = -D TABLE # Laine EoS, tabulated

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -march=native
LD            = g++
LDFLAGS       = -O3 -march=native

CXXFLAGS     += $(ROOTCFLAGS) $(EXTRA_FLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)

vpath %.cpp src
objdir     = obj

SRC        = Function.cpp main.cpp Eccentricity.cpp
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
              
TARGET	   = a
#-------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)
	
obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
