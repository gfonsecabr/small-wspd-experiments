# Experiments in R^3 (separate from the planar code; nothing is shared). exp3d uses CPLEX 22.2 (Concert) for OPT;
# check_lb3 and lb_swap3 need no library.
#   make            # build all
CPLEX   ?= /opt/ibm/ILOG/CPLEX_Studio222
CXX     ?= g++
CXXFLAGS ?= -std=c++23 -O2 -DNDEBUG -Wall -Wextra -DIL_STD
INCLUDES := -I$(CPLEX)/cplex/include -I$(CPLEX)/concert/include
LIBS     := -L$(CPLEX)/cplex/lib/x86-64_linux/static_pic -L$(CPLEX)/concert/lib/x86-64_linux/static_pic \
            -lilocplex -lconcert -lcplex -lm -lpthread -ldl
.PHONY: all clean
all: exp3d check_lb3 lb_swap3
exp3d: exp3d.cpp wspd3.hpp ls3.hpp partition3.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LIBS)
check_lb3: check_lb3.cpp
	$(CXX) -std=c++17 -O2 -Wall -Wextra -o $@ $<
lb_swap3: lb_swap3.cpp
	$(CXX) -std=c++17 -O2 -Wall -Wextra -o $@ $<
clean:
	rm -f exp3d check_lb3 lb_swap3
