CXX = mpicxx
CXXFLAGS = -O3
LDFLAGS = -L/opt/mkl/mkl/lib/intel64
LDLIBS = -lmkl_rt
INCLUDES = -I/opt/mkl/mkl/include

SRC_DIR = src
OBJ_DIR = $(SRC_DIR)/obj
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))
EXECUTABLE = H2DTOMO

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDES)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
	rm -rf $(OBJ_DIR)
