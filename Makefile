
CXX = g++
CXXFLAGS = -Wall -std=c++17

SOURCES = src/structures/structures.cpp src/gates/gates.cpp src/Experiments/myLab.cpp src/main.cpp
OBJECTS = $(SOURCES:src/%.cpp=build/%.o)
TARGET = build/app

all: $(TARGET)

build:
	mkdir -p build/gates build/structures build/Experiments

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

clean:
	rm -rf build/app.o build/structures/*.o build/gates/*.o build/Experiments/*.o

.PHONY: all clean

