SOURCES := $(shell find src -path src/bindings -prune -o -type f -name '*.cpp' -print)
FORMAT_FILES := $(shell find src include tests benchmarks -type f \( -name '*.cpp' -o -name '*.hpp' \))
BENCH_ARGS ?= 2048 500000 50

.PHONY: release debug python flamegraph lint format

release debug python:
	cmake --workflow --preset $@

flamegraph:
	cmake --workflow --preset profile
	cd build/profile && flamegraph -c "record -F 997 -e cpu-clock:u --call-graph dwarf" --post-process "tee flamegraph.folded" -o flamegraph.svg -- ./stabaal_bench $(BENCH_ARGS)
	cp build/profile/flamegraph.svg flamegraph.svg
	cp build/profile/flamegraph.folded flamegraph.folded

lint:
	cmake --preset debug
	clang-tidy -p build/debug $(SOURCES)

format:
	clang-format -i $(FORMAT_FILES)

dry-format:
	clang-format --dry-run $(FORMAT_FILES)

test:
	cmake --workflow --preset $@
	build/test/stabaal_tests