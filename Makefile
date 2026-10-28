CC = g++
#CFLAGS = -std=c++17 -O2 -DNDEBUG -pedantic-errors -Wall -Wextra -Wconversion -Wno-sign-compare -Iinclude -fsanitize=address -g
CFLAGS = -std=c++17 -O2 -DNDEBUG -pedantic-errors -Wall -Wextra -Wconversion -Wno-sign-compare -Iinclude 

SRC1_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

SRC_FILES = $(wildcard $(SRC1_DIR)/*.cpp) $(wildcard $(SRC1_DIR)/*.c)
OBJ_FILES = $(patsubst $(SRC1_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(filter %.cpp,$(SRC_FILES))) \
            $(patsubst $(SRC1_DIR)/%.c,$(BUILD_DIR)/%.o,$(filter %.c,$(SRC_FILES)))

TARGET = deknot
LIB_TARGET = libdeknot.a

$(TARGET): $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $(OBJ_FILES)

$(LIB_TARGET): $(filter-out $(BUILD_DIR)/main.o,$(OBJ_FILES))
	ar rcs $@ $^

lib: $(LIB_TARGET)

$(BUILD_DIR)/%.o: $(SRC1_DIR)/%.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC1_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

-include $(OBJ_FILES:.o=.d)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# ---------------------------------------------------------------------------
# WebAssembly build (requires Emscripten: https://emscripten.org/docs/getting_started/)
# Produces web/deknot.js + web/deknot.wasm, loaded by web/index.html via web/worker.js
# ---------------------------------------------------------------------------
EMCC = em++
WASM_DIR = web
WASM_SRC = $(filter-out $(SRC1_DIR)/main.cpp,$(SRC_FILES)) $(WASM_DIR)/wasm_main.cpp
WASM_FLAGS = -std=c++17 -O2 -DNDEBUG -Wall -Wextra -Wno-sign-compare -Iinclude \
	-sMODULARIZE=1 -sEXPORT_NAME=createDeKnot -sENVIRONMENT=web,worker,node \
	-sALLOW_MEMORY_GROWTH=1 -sMAXIMUM_MEMORY=4GB -sSTACK_SIZE=8MB \
	-sINVOKE_RUN=0 -sEXIT_RUNTIME=0 -sFILESYSTEM=1 \
	-sEXPORTED_FUNCTIONS=_dk_wasm_assemble,_malloc,_free \
	-sEXPORTED_RUNTIME_METHODS=FS,ccall,cwrap

wasm: $(WASM_DIR)/deknot.js

$(WASM_DIR)/deknot.js: $(WASM_SRC) $(wildcard $(INCLUDE_DIR)/*.h)
	$(EMCC) $(WASM_FLAGS) -o $@ $(WASM_SRC)

wasm-clean:
	rm -f $(WASM_DIR)/deknot.js $(WASM_DIR)/deknot.wasm
