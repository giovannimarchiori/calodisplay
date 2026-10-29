# Select compiler command depending on environment
OS:=$(shell uname)
CXX=$(shell root-config --cxx)
ifeq ($(OS),Darwin)
	CXX=clang++
endif

# Variables for compiling and linking
CXXFLAGS=$(shell root-config --cflags) -fPIC
LDFLAGS=$(shell root-config --ldflags)
LIBS=$(shell root-config --glibs --evelibs)
# OpenGL is called directly (palette legend): on macOS the framework must be linked explicitly
ifeq ($(OS),Darwin)
	LIBS+=-framework OpenGL
endif
INCDIR=$(shell root-config --incdir)

# Define variables for directories
SRC_DIR=src
OBJ_DIR=build
BIN_DIR=dist

# Set the name of your executable
APP_NAME=calodisplay

# Dictionary .cxx (generated in build/) and .pcm (must sit next to the shared library)
DICT_NAME=$(APP_NAME)-dictionary
DICT_CXX=$(OBJ_DIR)/$(DICT_NAME).cxx
DICT_OBJ=$(OBJ_DIR)/$(DICT_NAME).o
DICT_PCM_FILENAME=$(DICT_NAME)_rdict.pcm

# Shared library (all code except main + dictionary) and its debug symbols
SHARED_LIBRARY=$(APP_NAME)-library.so
SHARED_LIBRARY_PATH=$(BIN_DIR)/$(SHARED_LIBRARY)
SHARED_LIBRARY_DS=$(SHARED_LIBRARY).dSYM

# the executable finds the library by its bare name (via LD_LIBRARY_PATH,
# or via the rpath for the -local executable)
ifeq ($(OS),Darwin)
	SONAME_FLAG=-Wl,-install_name,$(SHARED_LIBRARY)
else
	SONAME_FLAG=-Wl,-soname,$(SHARED_LIBRARY)
endif

# File recording the ROOT version used for the build (checked by setup.sh)
BUILD_INFO=$(BIN_DIR)/build-info.txt

# Binary executable file names with path
EXECUTABLE=$(BIN_DIR)/$(APP_NAME)
EXECUTABLE_LOCAL=$(BIN_DIR)/$(APP_NAME)-local

# Define list of header files (.h), source files (.cpp) and object files (.o)
H_EXT = h
HEADERS = $(shell find $(SRC_DIR) -type f -name "*.$(H_EXT)")
HEADERS := $(filter-out $(SRC_DIR)/LinkDef.h,$(HEADERS))

SRC_EXT = cpp
SOURCES = $(shell find $(SRC_DIR) -type f -name "*.$(SRC_EXT)")

OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
MAIN_OBJECT = $(OBJ_DIR)/main.o
LIB_OBJECTS = $(filter-out $(MAIN_OBJECT),$(OBJECTS))

# Convenience variable for creating new directories
dir_guard=@mkdir -p $(@D)

# Path for installing the shared library and dictionary .pcm
# Tip: macOS system integrity protection blocks access to $(DYLD_LIBRARY_PATH)
#      so instead we use $(ROOTSYS)/lib
# https://github.com/nteract/nteract/issues/1523#issuecomment-284027093
DYNAMIC_LIBRARY_PATH = $(ROOTSYS)/lib

# Path for installing the executable binary
# https://stackoverflow.com/questions/39892692/how-to-implement-make-install-in-a-makefile
ifeq ($(PREFIX),)
	PREFIX := /usr/local
endif

# Default target builds production executable
all: production

# Add -O3 optimization level for the production release
# Note: objects are not rebuilt when switching between production/release/debug,
# run 'make clean' first
production: CXXFLAGS+=-O3
production: $(EXECUTABLE)

# Add -O3 optimization level for the test release
release: CXXFLAGS+=-O3
release: $(EXECUTABLE_LOCAL)

# Also might add flags for debug optimizations: -Og -ggdb -DDEBUG
debug: CXXFLAGS+=-g
debug: $(EXECUTABLE_LOCAL) move_debug_symbols

# Target for dictionary generation (creates .cxx and .pcm in build/)
$(DICT_CXX): $(HEADERS) $(SRC_DIR)/LinkDef.h
	$(dir_guard)
	rootcling -f $@ $^

# Compile the dictionary (it includes the headers as src/...)
$(DICT_OBJ): $(DICT_CXX)
	@echo "Compiling "$@
	$(CXX) -c $(CXXFLAGS) -I. -I$(SRC_DIR) $< -o $@

# Target for compiling the object files, with header dependency files (.d)
# (the generated headers are in build/)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(dir_guard)
	@echo "Compiling "$@
	$(CXX) -c $(CXXFLAGS) -I$(OBJ_DIR) -MMD -MP $< -o $@

# Header with the hash of the code that decides which volumes of the geometry
# are drawn, which identifies the geometry caches it wrote: rewritten only if
# the hash changes
GEOMETRY_CODE_VERSION=$(OBJ_DIR)/GeometryCodeVersion.h
$(GEOMETRY_CODE_VERSION): $(SRC_DIR)/EventDisplay.cpp $(SRC_DIR)/GeometryCache.cpp scripts/geometry_code_version.sh
	$(dir_guard)
	@scripts/geometry_code_version.sh $(SRC_DIR) > $@.tmp
	@if cmp -s $@.tmp $@; then rm -f $@.tmp; else mv $@.tmp $@; fi
$(OBJ_DIR)/EventDisplay.o: $(GEOMETRY_CODE_VERSION)

# Target for linking the shared library from the object files and the dictionary
# Official ROOT docs are missing $(LIBS) parameter https://root.cern.ch/interacting-shared-libraries-rootcling
# The dictionary .pcm is copied next to the library, where ROOT looks for it
$(SHARED_LIBRARY_PATH): $(LIB_OBJECTS) $(DICT_OBJ)
	$(dir_guard)
	@echo "Linking "$@
	$(CXX) -shared -o $@ $(LDFLAGS) $(SONAME_FLAG) $^ $(LIBS)
	cp $(OBJ_DIR)/$(DICT_PCM_FILENAME) $(BIN_DIR)/$(DICT_PCM_FILENAME)
	@echo "ROOT_VERSION=$$(root-config --version)" > $(BUILD_INFO)
	@echo "ROOT_PREFIX=$$(root-config --prefix)" >> $(BUILD_INFO)
	@echo "CXX=$$($(CXX) --version | head -1)" >> $(BUILD_INFO)

# Target for linking the production executable (main only, rest from the library)
$(EXECUTABLE): $(MAIN_OBJECT) $(SHARED_LIBRARY_PATH)
	@echo "Linking "$@
	$(CXX) -o $@ $(MAIN_OBJECT) $(SHARED_LIBRARY_PATH) $(LIBS)

# Target for linking the test release and development executable
# We assign relative .so search path (./) to the executable
$(EXECUTABLE_LOCAL): $(MAIN_OBJECT) $(SHARED_LIBRARY_PATH)
	@echo "Linking "$@
ifeq ($(OS),Darwin)
	# for macOS just link against the shared library
	$(CXX) -o $@ $(MAIN_OBJECT) $(SHARED_LIBRARY_PATH) $(LIBS)
	# then change search location of the .so library in the executable - set as same directory
	install_name_tool -change $(SHARED_LIBRARY) @executable_path/$(SHARED_LIBRARY) $@
else
	# for Linux add runtime shared library search path ./ relative to the executable
	# https://stackoverflow.com/questions/38058041/correct-usage-of-rpath-relative-vs-absolute
	$(CXX) -o $@ $(MAIN_OBJECT) $(SHARED_LIBRARY_PATH) $(LIBS) -Wl,-rpath,'$$ORIGIN'
endif

# Target creates the .so.dSYM debug symbols in /dist folder
move_debug_symbols: $(SHARED_LIBRARY_PATH)
ifeq ($(OS),Darwin)
	dsymutil $(SHARED_LIBRARY_PATH) -o $(BIN_DIR)/$(SHARED_LIBRARY_DS)
endif

# Target installs the executable and moves shared library to $ROOTSYS/lib
install:
	# run 'make install' or 'sudo -E make install' to preserve user environment
	sudo install -m 755 $(EXECUTABLE) $(DESTDIR)$(PREFIX)/bin/
	sudo install -m 755 $(SHARED_LIBRARY_PATH) $(DYNAMIC_LIBRARY_PATH)/
	sudo install -m 755 $(BIN_DIR)/$(DICT_PCM_FILENAME) $(DYNAMIC_LIBRARY_PATH)/

# Target to copy the application icon and create the launcher on Linux
#install-linux-launcher:
#	xdg-icon-resource install --context apps --size 128 ./resources/$(APP_NAME).png $(APP_NAME)
#	xdg-desktop-menu install ./resources/$(APP_NAME).desktop

# Target that cleans the build
clean:
	rm -f -r $(OBJ_DIR)
	rm -f -r $(BIN_DIR)
	rm -f $(DICT_NAME).cxx $(DICT_PCM_FILENAME) $(SHARED_LIBRARY)
	rm -f -r $(SHARED_LIBRARY_DS)
	rm -f *~ src/*~
	rm -f src/*_cpp.* src/*.pcm
	rm -f *_dict*.h

# Target for test output
echo:
	$(info SOURCES: $(SOURCES))
	$(info HEADERS: $(HEADERS))
	$(info LIB_OBJECTS: $(LIB_OBJECTS))

# Header dependencies generated by -MMD
-include $(OBJECTS:.o=.d)

# List of special targets that do not generate files
.PHONY: all production release debug clean move_debug_symbols install echo
