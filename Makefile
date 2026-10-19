SOURCES := $(shell find src -type f -name '*.cpp')
FORMAT_FILES := $(shell find src include tests benchmarks -type f \( -name '*.cpp' -o -name '*.hpp' \))
BENCH_ARGS ?= 256 1000000 10

.PHONY: release debug flamegraph lint format

release debug:
	cmake --workflow --preset $@

flamegraph:
	cmake --workflow --preset profile
	cd build/profile && flamegraph -c "record -F 997 -e cpu-clock:u --call-graph dwarf" -o flamegraph.svg -- ./stabaal_bench $(BENCH_ARGS)
	cp build/profile/flamegraph.svg flamegraph.svg

lint:
	cmake --preset debug
	clang-tidy -p build/debug $(SOURCES)

format:
	clang-format -i $(FORMAT_FILES)
