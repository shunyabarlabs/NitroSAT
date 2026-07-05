# NitroSAT Makefile
# ==================
# Supports building v1, v2, and disk-streaming v3 on various hardware configurations
#
# Usage:
#   make                    # Build all (v1 + v2 + v3)
#   make v1                # Build only v1
#   make v2                # Build only v2
#   make v3                # Build disk-streaming v3
#   make install           # Install to ~/.local/bin
#   make clean             # Remove build artifacts
#
# Options:
#   CC=gcc|clang          # Compiler (default: gcc)
#   march=x86-64|armv8|znver4|skylake|native  # CPU architecture
#   OPT=release|debug|size # Optimization level
#   CFLAGS=...            # Additional C flags
#
# Examples:
#   make v1 march=znver4          # Build v1 for AMD Zen 4
#   make CC=clang v2               # Build v2 with clang
#   make v1 OPT=debug              # Build v1 with debug symbols

# Default values
CC ?= gcc
OPT ?= release
PREFIX ?= $(HOME)/.local/bin
BUILD_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Source files
V1_SRC = $(BUILD_DIR)v1/nitrosat.c
V2_SRC = $(BUILD_DIR)v2/nitrosatv2.c
V3_SRC = $(BUILD_DIR)v3/nitrosatv3.c

# Output binaries
V1_BIN = $(BUILD_DIR)v1/nitrosat
V2_BIN = $(BUILD_DIR)v2/nitrosatv2
V3_BIN = $(BUILD_DIR)v3/nitrosatv3

# Architecture flags
MARCH = native

# Optimization flags
ifeq ($(OPT),debug)
    CFLAGS += -O0 -g -DNITROSAT_DEBUG
else ifeq ($(OPT),size)
    CFLAGS += -Os
else
    # release (default)
    CFLAGS += -O3
endif

# Add march flag
CFLAGS += -march=$(MARCH) -std=c99

# Common libraries
LDFLAGS = -lm

# Targets
.PHONY: all v1 v2 v3 test-v3 clean install uninstall help

all: v1 v2 v3

v1: $(V1_BIN)

v2: $(V2_BIN)

v3: $(V3_BIN)

test-v3: $(V3_BIN)
	bash $(BUILD_DIR)../../tests/v3/test_v3.sh $(V3_BIN)

$(V1_BIN): $(V1_SRC)
	@echo "Building NitroSAT v1..."
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)
	@echo "Built: $@"

$(V2_BIN): $(V2_SRC)
	@echo "Building NitroSAT v2..."
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)
	@echo "Built: $@"

$(V3_BIN): $(V3_SRC)
	@echo "Building NitroSAT v3 (disk streaming)..."
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)
	@echo "Built: $@"

clean:
	rm -f $(V1_BIN) $(V2_BIN) $(V3_BIN)
	@echo "Cleaned build artifacts"

install: v1 v2 v3
	@mkdir -p $(PREFIX)
	cp $(V1_BIN) $(PREFIX)/nitrosat
	cp $(V2_BIN) $(PREFIX)/nitrosatv2
	cp $(V3_BIN) $(PREFIX)/nitrosatv3
	@echo "Installed to $(PREFIX): nitrosat, nitrosatv2, nitrosatv3"

uninstall:
	rm -f $(PREFIX)/nitrosat $(PREFIX)/nitrosatv2 $(PREFIX)/nitrosatv3
	@echo "Uninstalled from $(PREFIX)"

help:
	@echo "NitroSAT Makefile"
	@echo "================="
	@echo ""
	@echo "Targets:"
	@echo "  all         Build v1, v2, and v3 (default)"
	@echo "  v1          Build NitroSAT v1"
	@echo "  v2          Build NitroSAT v2"
	@echo "  v3          Build disk-streaming NitroSAT v3"
	@echo "  test-v3     Run V3 end-to-end tests"
	@echo "  clean       Remove build artifacts"
	@echo "  install     Install to ~/.local/bin"
	@echo "  uninstall   Remove installed binaries"
	@echo "  help        Show this help"
	@echo ""
	@echo "Options:"
	@echo "  CC=gcc|clang      Compiler (default: gcc)"
	@echo "  march=...         CPU architecture (default: native)"
	@echo "  OPT=release|debug|size  Optimization (default: release)"
	@echo "  CFLAGS=...        Additional C flags"
	@echo ""
	@echo "Architecture options: native, x86-64, armv8-a, znver4, skylake"
	@echo ""
	@echo "Examples:"
	@echo "  make v1 march=znver4"
	@echo "  make CC=clang v2"
	@echo "  make v1 OPT=debug"
