# Builds the two binaries of the experiments, the lower bound, and its checker:
#   exp        experiment driver (sweep, Callahan-Kosaraju, exact OPT, local search, partition WSPDs)
#   optx       exact OPT for larger n (set-cover reductions + CPLEX per component)
#   lb_swap    the lower bound of the paper: local search on the anchors of the sweep (no dependency)
#   check_lb   independent check of a lower bound (no dependency: make check_lb needs only a compiler)
# Requirements of exp and optx: a C++23 compiler, IBM ILOG CPLEX 22.2 (Concert), and the header-only pgl library (see
# ../get_deps.sh).
#   make                         # build all four
#   make CPLEX=/path/to/CPLEX_Studio222 PGL=/path/to/pgl/include
CPLEX   ?= /opt/ibm/ILOG/CPLEX_Studio222
PGL     ?= ../deps/pgl/include
CXX     ?= g++
CXXFLAGS ?= -std=c++23 -O2 -DNDEBUG -Wall -Wextra -DIL_STD
INCLUDES := -I$(PGL) -I$(CPLEX)/cplex/include -I$(CPLEX)/concert/include
LIBS     := -L$(CPLEX)/cplex/lib/x86-64_linux/static_pic -L$(CPLEX)/concert/lib/x86-64_linux/static_pic \
            -lilocplex -lconcert -lcplex -lm -lpthread -ldl
HEADERS := wspd.hpp ls.hpp partition.hpp inputs.hpp
.PHONY: all clean
all: exp optx lb_swap check_lb
check_lb lb_swap: %: %.cpp
	$(CXX) -std=c++17 -O2 -Wall -Wextra -o $@ $<
%: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LIBS)
clean:
	rm -f exp optx lb_swap check_lb
