# Site dependant: external libraries. Must be configured by user.
CPLEXROOT := /home/ibm/cplex-studio/22.1.1
EIGENROOT := /home/gencol/prebuilt/x86_64-pc-linux_el9-gnu/eigen-3.4.0

SRC := src
BUILD := build

EIGENINCPATH := $(EIGENROOT)
CPLEXINCPATH := $(CPLEXROOT)/concert/include \
                $(CPLEXROOT)/cplex/include
CPLEXLIBPATH := $(CPLEXROOT)/concert/lib/x86-64_linux/static_pic \
                $(CPLEXROOT)/cplex/lib/x86-64_linux/static_pic \

CXX := g++
CPPFLAGS := -DIL_STD -DNDEBUG \
            $(patsubst %,-I%,$(CPLEXINCPATH) $(EIGENINCPATH))
CXXFLAGS := -std=gnu++17 -O3 \
            -fno-math-errno -msse2 -mfpmath=sse -m64 -finline-functions \
            -fexceptions -fno-strict-aliasing \
            -Wpedantic -Wreturn-type -Wunused -Wswitch -Wcomment \
            -Wtrigraphs -Wformat -Wchar-subscripts -Wparentheses \
            -Wwrite-strings -Woverloaded-virtual -Wuninitialized
CXXDEPFLAGS := -MMD -MP

LDFLAGS := $(patsubst %,-L%,$(CPLEXLIBPATH))
LDLIBS := -lconcert -lilocplex -lcplex -ldl -lm -lrt -lpthread

PRG := main
SRCLIBFILES := $(shell find $(SRC) -type f -name '*.cpp' \
                                   -not -name "$(PRG).cpp" -print)
SRCFILES := $(SRC)/$(PRG).cpp $(SRCLIBFILES)
OBJFILES := $(patsubst $(SRC)/%.cpp,$(BUILD)/%.o, $(SRCFILES))
DEPFILES := $(patsubst %.o,%.d, $(OBJFILES))

.PHONY: all
all: $(PRG)

.PHONY: clean
clean:
	rm -rf $(BUILD) $(PRG)

$(PRG): $(BUILD)/$(PRG)
	ln -s -f $< $@

$(BUILD)/$(PRG): $(OBJFILES)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

.PHONY: obj
obj: $(OBJFILES)

$(BUILD)/%.o: $(SRC)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXXDEPFLAGS) -c $< -o $@

-include $(DEPFILES)
