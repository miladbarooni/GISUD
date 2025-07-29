# Detect platform
UNAME_S := $(shell uname -s)

# Set default paths based on the platform
ifeq ($(UNAME_S),Linux)
    CPLEXROOT := /home/ibm/cplex-studio/22.1.1
    EIGENROOT := /home/gencol/prebuilt/x86_64-pc-linux_el9-gnu/eigen-3.4.0
    CPLEXLIBSUBDIR := x86-64_linux
else ifeq ($(UNAME_S),Darwin)
    CPLEXROOT := /Applications/CPLEX_Studio2211
    EIGENROOT := /opt/homebrew/include/eigen3
    CPLEXLIBSUBDIR := arm64_osx
else
    $(error Unsupported platform: $(UNAME_S))
endif

SRC := src
BUILD := build
PRG := main

# Include and Library Paths
EIGENINCPATH := $(EIGENROOT)
CPLEXINCPATH := $(CPLEXROOT)/concert/include \
                $(CPLEXROOT)/cplex/include
CPLEXLIBPATH := $(CPLEXROOT)/concert/lib/$(CPLEXLIBSUBDIR)/static_pic \
                $(CPLEXROOT)/cplex/lib/$(CPLEXLIBSUBDIR)/static_pic

# Compiler and Flags
CXX := g++
CPPFLAGS := -DIL_STD -DNDEBUG \
            $(patsubst %,-I%,$(CPLEXINCPATH) $(EIGENINCPATH))
CXXFLAGS := -std=gnu++17 -O3 \
            -fno-math-errno -m64 -finline-functions \
            -fexceptions -fno-strict-aliasing \
            -Wpedantic -Wreturn-type -Wunused -Wswitch -Wcomment \
            -Wtrigraphs -Wformat -Wchar-subscripts -Wparentheses \
            -Wwrite-strings -Woverloaded-virtual -Wuninitialized
CXXDEPFLAGS := -MMD -MP

LDFLAGS := $(patsubst %,-L%,$(CPLEXLIBPATH))
LDLIBS := -lconcert -lilocplex -lcplex -ldl -lm -lpthread

# Source and Object files
SRCLIBFILES := $(shell find $(SRC) -type f -name '*.cpp' -not -name "$(PRG).cpp")
SRCFILES := $(SRC)/$(PRG).cpp $(SRCLIBFILES)
OBJFILES := $(patsubst $(SRC)/%.cpp,$(BUILD)/%.o, $(SRCFILES))
DEPFILES := $(patsubst %.o,%.d, $(OBJFILES))

.PHONY: all clean obj

all: $(PRG)

clean:
	rm -rf $(BUILD) $(PRG)

$(PRG): $(BUILD)/$(PRG)
	ln -s -f $< $@

$(BUILD)/$(PRG): $(OBJFILES)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

obj: $(OBJFILES)

$(BUILD)/%.o: $(SRC)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXXDEPFLAGS) -c $< -o $@

-include $(DEPFILES)
