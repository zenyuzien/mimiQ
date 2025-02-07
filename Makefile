# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -std=c++17

# Source and object files
SOURCES = src/structures/structures.cpp src/gates/gates.cpp src/Experiments/myLab.cpp src/main.cpp
OBJECTS = $(SOURCES:src/%.cpp=build/%.o)
TARGET = build/app

# Default target
all: $(TARGET)

# Create build directories if they don't exist
build:
	mkdir -p build/gates build/structures build/Experiments

# Compile object files from source files
build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link the target executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

# Clean build
clean:
	rm -rf build/app.o build/structures/*.o build/gates/*.o build/Experiments/*.o

# Phony targets
.PHONY: all clean

