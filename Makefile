CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic
LDFLAGS  ?=

SRC := $(wildcard src/*.cpp)
OBJ := $(patsubst src/%.cpp,build/%.o,$(SRC))
DEP := $(OBJ:.o=.d)
PY  := .venv/bin/python

crnsim: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -MMD -MP -c -o $@ $<

build:
	@mkdir -p build

test: crnsim
	@$(PY) tests/test_crnsim.py

venv:
	python3 -m venv .venv && .venv/bin/pip install numpy matplotlib pandas

clean:
	rm -rf build crnsim

.PHONY: clean test venv

-include $(DEP)
