CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic
LDFLAGS  ?=
LDLIBS   ?= -pthread

SRC := $(wildcard src/*.cpp)
OBJ := $(patsubst src/%.cpp,build/%.o,$(SRC))
DEP := $(OBJ:.o=.d)
PY  := .venv/bin/python

crnsim: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -MMD -MP -c -o $@ $<

build:
	@mkdir -p build

test: crnsim
	@$(PY) tests/test_crnsim.py

bench: crnsim
	@./scripts/bench.sh

venv:
	python3 -m venv .venv && .venv/bin/pip install numpy matplotlib pandas

clean:
	rm -rf build crnsim

.PHONY: clean test bench venv

-include $(DEP)
