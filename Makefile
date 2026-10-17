SOURCES := $(shell find src -type f -name '*.cpp')
FORMAT_FILES := $(shell find src include tests -type f \( -name '*.cpp' -o -name '*.hpp' \))

.PHONY: release debug lint format

release debug:
	cmake --workflow $@

lint:
	cmake --preset debug
	clang-tidy -p build/debug $(SOURCES)

format:
	clang-format -i $(FORMAT_FILES)
