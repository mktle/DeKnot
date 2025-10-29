CC = g++
CFLAGS = -std=c++17 -O2 -DNEDEBUG -pedantic-errors -Wall -Wextra -Wconversion -Wno-sign-compare -Iinclude -fsanitize=address -g
#CFLAGS = -std=c++17 -O2 -DNEDEBUG -pedantic-errors -Wall -Wextra -Wconversion -Wno-sign-compare -Iinclude 

# Directories
SRC1_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

SRC_FILES = $(wildcard $(SRC1_DIR)/*.cpp) $(wildcard $(SRC1_DIR)/*.c)
OBJ_FILES = $(patsubst $(SRC1_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(filter %.cpp,$(SRC_FILES))) \
            $(patsubst $(SRC1_DIR)/%.c,$(BUILD_DIR)/%.o,$(filter %.c,$(SRC_FILES)))

TARGET = deknot

# Rule to build the target
$(TARGET): $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $(OBJ_FILES)

$(BUILD_DIR)/%.o: $(SRC1_DIR)/%.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC1_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET)
