#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "stabaal/simulator.hpp"
#include "stabaal/types.hpp"

namespace {

using stabaal::QubitId;
using stabaal::QubitSize;

// Exercises id routing without implementing quantum discard or compaction in the real backend.
struct IndexBackend {
    explicit IndexBackend(const QubitSize count)
        : alive(count, true), ids(count), active_count(count), next_id(count) {
        for (QubitSize i{0}; i < count; ++i) ids[i] = i;
        compact = false;
        fail_add = false;
        fail_discard = false;
    }

    [[nodiscard]] QubitSize NumQubits() const noexcept { return active_count; }

    QubitId Add() {
        if (fail_add) throw std::runtime_error{"Backend add failed"};
        const auto free{std::find(alive.begin(), alive.end(), false)};
        if (free != alive.end()) {
            const auto index{static_cast<QubitSize>(free - alive.begin())};
            alive[index] = true;
            ++active_count;
            ids[index] = next_id++;
            return ids[index];
        }
        const auto index{static_cast<QubitSize>(alive.size())};
        alive.push_back(true);
        ids.push_back(next_id++);
        ++active_count;
        return ids[index];
    }

    void Discard(const QubitId id) {
        if (fail_discard) throw std::runtime_error{"Backend discard failed"};
        const auto q{Resolve(id)};
        alive[q] = false;
        --active_count;
        if (compact) {
            QubitSize destination{0};
            for (QubitSize source{0}; source < alive.size(); ++source) {
                if (!alive[source]) continue;
                if (source != destination) {
                    ids[destination] = ids[source];
                    alive[destination] = true;
                    alive[source] = false;
                }
                ++destination;
            }
            alive.resize(active_count);
            ids.resize(active_count);
        }
    }

    QubitSize Resolve(const QubitId id) const {
        const auto found{std::find(ids.begin(), ids.end(), id)};
        if (found == ids.end()) throw std::out_of_range{"Qubit id is not active"};
        const auto index{static_cast<QubitSize>(found - ids.begin())};
        if (!alive[index]) throw std::out_of_range{"Qubit id is not active"};
        return index;
    }

    void Record(const QubitId q) const { last_q = Resolve(q); }
    void Record(const QubitId c, const QubitId t) const {
        const auto control{Resolve(c)};
        const auto target{Resolve(t)};
        last_q = control;
        last_target = target;
    }

    bool MeasureX(const QubitId q) {
        Record(q);
        return true;
    }
    bool MeasureZ(const QubitId q) {
        Record(q);
        return false;
    }
    bool IsDeterministicMeasureX(const QubitId q) const {
        Record(q);
        return true;
    }
    bool IsDeterministicMeasureZ(const QubitId q) const {
        Record(q);
        return false;
    }
    bool ForceMeasureX(const QubitId q, const bool outcome) {
        Record(q);
        return outcome;
    }
    bool ForceMeasureZ(const QubitId q, const bool outcome) {
        Record(q);
        return outcome;
    }
    void ApplyH(const QubitId q) { Record(q); }
    void ApplyS(const QubitId q) { Record(q); }
    void ApplySdg(const QubitId q) { Record(q); }
    void ApplySx(const QubitId q) { Record(q); }
    void ApplyI(const QubitId q) { Record(q); }
    void ApplyX(const QubitId q) { Record(q); }
    void ApplyZ(const QubitId q) { Record(q); }
    void ApplyY(const QubitId q) { Record(q); }
    void ApplyCx(const QubitId c, const QubitId t) { Record(c, t); }
    void ApplyCz(const QubitId c, const QubitId t) { Record(c, t); }
    void ApplyCy(const QubitId c, const QubitId t) { Record(c, t); }
    void Print() const {}
    [[nodiscard]] std::size_t MemoryUsage() const noexcept { return sizeof(*this); }

    std::vector<bool> alive;
    std::vector<QubitId> ids;
    QubitSize active_count;
    QubitId next_id;
    static inline QubitSize last_q{};
    static inline QubitSize last_target{};
    static inline bool compact{};
    static inline bool fail_add{};
    static inline bool fail_discard{};
};

using IndexSim = stabaal::Simulator<IndexBackend>;

}  // namespace

TEST_CASE("New qubits return IDs that survive growth", "[simulator][ids]") {
    stabaal::Simulator<> sim{2};
    const auto added{sim.Add()};
    CHECK(added == QubitId{2});
    sim.X(added);
    for (unsigned i{0}; i < 150; ++i) sim.Add();
    CHECK(sim.MeasureZ(added));
    CHECK_FALSE(sim.MeasureZ(0));
    CHECK_FALSE(sim.MeasureZ(1));
    CHECK(sim.NumQubits() == 153);
    stabaal::detail::InverseTableau backend{153};
    CHECK(sizeof(sim) == sizeof(backend));
}

TEST_CASE("Append-only simulators allocate no id arrays", "[simulator][ids][memory]") {
    for (const auto n : {QubitSize{0}, QubitSize{1}, QubitSize{64}, QubitSize{129}}) {
        stabaal::Simulator<> sim{n};
        stabaal::detail::InverseTableau backend{n};
        const auto object_overhead{sizeof(sim) - sizeof(backend)};
        CHECK(sim.MemoryUsage() == backend.MemoryUsage() + object_overhead);
        for (unsigned i{0}; i < 150; ++i) {
            CHECK(sim.Add() == backend.Add());
            CHECK(sim.MemoryUsage() == backend.MemoryUsage() + object_overhead);
        }
        auto copy{sim};
        const auto copied_backend{backend};
        CHECK(copy.MemoryUsage() == copied_backend.MemoryUsage() + object_overhead);
        copy = sim;
        CHECK(copy.MemoryUsage() == copied_backend.MemoryUsage() + object_overhead);
    }
}

TEST_CASE("Reusing a middle storage slot gives a fresh id for every operation",
          "[simulator][ids]") {
    IndexSim sim{3};
    sim.Discard(1);
    const auto replacement{sim.Add()};
    CHECK(replacement != 1);
    CHECK(replacement != 2);
    CHECK(sim.NumQubits() == 3);

    const auto routes = [&](auto operation) {
        IndexBackend::last_q = std::numeric_limits<QubitSize>::max();
        operation();
        CHECK(IndexBackend::last_q == 1);
    };
    routes([&] { sim.H(replacement); });
    routes([&] { sim.S(replacement); });
    routes([&] { sim.Sdg(replacement); });
    routes([&] { sim.Sx(replacement); });
    routes([&] { sim.I(replacement); });
    routes([&] { sim.X(replacement); });
    routes([&] { sim.Y(replacement); });
    routes([&] { sim.Z(replacement); });
    routes([&] { CHECK(sim.MeasureX(replacement)); });
    routes([&] { CHECK_FALSE(sim.MeasureZ(replacement)); });
    routes([&] { CHECK(std::as_const(sim).IsDeterministicMeasureX(replacement)); });
    routes([&] { CHECK_FALSE(std::as_const(sim).IsDeterministicMeasureZ(replacement)); });
    routes([&] { CHECK(sim.ForceMeasureX(replacement, true)); });
    routes([&] { CHECK_FALSE(sim.ForceMeasureZ(replacement, false)); });
    routes([&] { sim.Cx(replacement, 2); });
    CHECK(IndexBackend::last_target == 2);
    routes([&] { sim.Cz(replacement, 2); });
    CHECK(IndexBackend::last_target == 2);
    routes([&] { sim.Cy(replacement, 2); });
    CHECK(IndexBackend::last_target == 2);
    sim.Cx(2, replacement);
    CHECK(IndexBackend::last_q == 2);
    CHECK(IndexBackend::last_target == 1);
    sim.H(2);
    CHECK(IndexBackend::last_q == 2);

    const auto rejects = [&](auto operation) {
        IndexBackend::last_q = std::numeric_limits<QubitSize>::max();
        CHECK_THROWS_AS(operation(), std::out_of_range);
        CHECK(IndexBackend::last_q == std::numeric_limits<QubitSize>::max());
    };
    rejects([&] { sim.H(1); });
    rejects([&] { sim.S(1); });
    rejects([&] { sim.Sdg(1); });
    rejects([&] { sim.Sx(1); });
    rejects([&] { sim.I(1); });
    rejects([&] { sim.X(1); });
    rejects([&] { sim.Y(1); });
    rejects([&] { sim.Z(1); });
    rejects([&] { sim.MeasureX(1); });
    rejects([&] { sim.MeasureZ(1); });
    rejects([&] { static_cast<void>(sim.IsDeterministicMeasureX(1)); });
    rejects([&] { static_cast<void>(sim.IsDeterministicMeasureZ(1)); });
    rejects([&] { sim.ForceMeasureX(1, false); });
    rejects([&] { sim.ForceMeasureZ(1, false); });
    rejects([&] { sim.Cx(1, replacement); });
    rejects([&] { sim.Cx(replacement, 1); });
    rejects([&] { sim.Cz(1, replacement); });
    rejects([&] { sim.Cz(replacement, 1); });
    rejects([&] { sim.Cy(1, replacement); });
    rejects([&] { sim.Cy(replacement, 1); });
    rejects([&] { sim.Discard(1); });
    CHECK(sim.NumQubits() == 3);
}

TEST_CASE("Compaction relocations preserve IDs across subsequent slot reuse", "[simulator][ids]") {
    IndexSim sim{4};
    IndexBackend::compact = true;
    sim.Discard(1);
    CHECK(sim.NumQubits() == 3);
    sim.H(2);
    CHECK(IndexBackend::last_q == 1);
    sim.H(3);
    CHECK(IndexBackend::last_q == 2);
    const auto added{sim.Add()};
    CHECK(added != 1);
    sim.H(added);
    CHECK(IndexBackend::last_q == 3);
    sim.Discard(2);
    CHECK_THROWS_AS(sim.H(2), std::out_of_range);
    sim.Cx(3, added);
    CHECK(IndexBackend::last_q == 1);
    CHECK(IndexBackend::last_target == 2);
    sim.H(0);
    CHECK(IndexBackend::last_q == 0);
}

TEST_CASE("Backend failures leave qubit identities unchanged", "[simulator][ids]") {
    IndexSim sim{1};
    IndexBackend::fail_add = true;
    CHECK_THROWS_AS(sim.Add(), std::runtime_error);
    CHECK(sim.NumQubits() == 1);
    IndexBackend::fail_add = false;
    CHECK(sim.Add() == QubitId{1});
    IndexBackend::fail_discard = true;
    CHECK_THROWS_AS(sim.Discard(0), std::runtime_error);
    sim.H(0);
    CHECK(IndexBackend::last_q == 0);
    CHECK(sim.NumQubits() == 2);
}

TEST_CASE("Copies preserve IDs and keep their lifetimes independent", "[simulator][ids][copy]") {
    IndexSim original{3};
    original.Discard(1);
    auto copy{original};
    const auto replacement{original.Add()};
    CHECK_THROWS_AS(copy.H(replacement), std::out_of_range);
    CHECK_THROWS_AS(copy.H(1), std::out_of_range);
    CHECK(copy.Add() == replacement);
    copy.Discard(replacement);
    original.H(replacement);
    CHECK(IndexBackend::last_q == 1);
    CHECK_THROWS_AS(copy.H(replacement), std::out_of_range);
    copy = original;
    auto moved{std::move(copy)};
    moved.H(replacement);
    CHECK(IndexBackend::last_q == 1);
}

TEST_CASE("Repeated slot reuse rejects stale IDs without growing storage", "[simulator][ids]") {
    stabaal::Simulator<> sim{3};
    QubitId current{1};
    sim.Discard(current);
    current = sim.Add();
    const auto memory{sim.MemoryUsage()};
    for (unsigned iteration{0}; iteration < 1000; ++iteration) {
        const auto previous{current};
        sim.Discard(current);
        CHECK_THROWS_AS(sim.MeasureZ(previous), std::out_of_range);
        current = sim.Add();
        CHECK(current != previous);
        CHECK_THROWS_AS(sim.MeasureZ(previous), std::out_of_range);
        CHECK_THROWS_AS(sim.MeasureZ(1), std::out_of_range);
        CHECK_FALSE(sim.MeasureZ(current));
        CHECK_FALSE(sim.MeasureZ(2));
        CHECK(sim.MemoryUsage() == memory);
    }
}

TEST_CASE("Compound backend gates validate both IDs before changing the state",
          "[simulator][ids][bounds]") {
    stabaal::detail::InverseTableau backend{1};
    CHECK_THROWS_AS(backend.ApplyCz(1, 0), std::out_of_range);
    CHECK(backend.IsDeterministicMeasureZ(0));
    CHECK_FALSE(backend.MeasureZ(0));
    backend.ApplyH(0);
    CHECK_THROWS_AS(backend.ApplyCy(1, 0), std::out_of_range);
    CHECK(backend.IsDeterministicMeasureX(0));
    CHECK_FALSE(backend.MeasureX(0));
}

TEST_CASE("Add prefers reusable pure holes over retained purification slots",
          "[simulator][ids][mixed]") {
    stabaal::Simulator<> sim{4};
    sim.H(0);
    sim.Cx(0, 1);
    sim.Discard(3);
    sim.Discard(1);
    const auto memory{sim.MemoryUsage()};
    const auto first{sim.Add()};
    CHECK(sim.MemoryUsage() == memory + 4 * sizeof(std::uint32_t));
    CHECK(sim.IsDeterministicMeasureZ(first));
    CHECK_FALSE(sim.MeasureZ(first));
    CHECK_FALSE(sim.IsDeterministicMeasureZ(0));
    const auto second{sim.Add()};
    CHECK(second != first);
    CHECK(sim.IsDeterministicMeasureZ(second));
    CHECK_FALSE(sim.MeasureZ(second));
    CHECK_FALSE(sim.IsDeterministicMeasureZ(0));
    CHECK_THROWS_AS(sim.H(1), std::out_of_range);
    CHECK_THROWS_AS(sim.H(3), std::out_of_range);
}

TEST_CASE("Discarding all qubits contracts storage and never revives old IDs",
          "[simulator][ids][discard][memory]") {
    for (const auto n : {QubitSize{1}, QubitSize{64}, QubitSize{65}, QubitSize{129}}) {
        stabaal::Simulator<> sim{n};
        const auto original_memory{sim.MemoryUsage()};
        // Descending order exercises the bitmap search hint when earlier holes appear.
        for (QubitId q{n}; q > 0;) {
            --q;
            sim.Discard(q);
            CHECK(sim.MemoryUsage() <= original_memory);
            CHECK_THROWS_AS(sim.H(q), std::out_of_range);
        }
        CHECK(sim.NumQubits() == 0);
        auto copy{sim};
        CHECK(copy.MemoryUsage() == sizeof(sim));
        std::vector<QubitId> replacements;
        for (QubitId q{0}; q < n; ++q) {
            const auto id{sim.Add()};
            replacements.push_back(id);
            CHECK(id != q);
            CHECK(copy.Add() == id);
            CHECK(sim.MemoryUsage() == copy.MemoryUsage());
            CHECK(sim.IsDeterministicMeasureZ(id));
            CHECK_FALSE(sim.MeasureZ(id));
            CHECK_THROWS_AS(sim.H(q), std::out_of_range);
        }
        // Generation counters must survive both storage growth and copies after reuse.
        const auto extra{sim.Add()};
        sim.X(extra);
        auto grown_copy{sim};
        grown_copy = sim;
        CHECK(grown_copy.MeasureZ(extra));
        for (const auto id : replacements) CHECK_FALSE(grown_copy.MeasureZ(id));
    }
}

TEST_CASE("Entangled discard and fresh append need no generation array",
          "[simulator][ids][mixed][memory]") {
    stabaal::Simulator<> sim{4};
    sim.H(0);
    sim.Cx(0, 1);
    const auto memory{sim.MemoryUsage()};
    sim.Discard(1);
    CHECK(sim.MemoryUsage() == memory);
    CHECK_FALSE(sim.IsDeterministicMeasureZ(0));
    const auto added{sim.Add()};
    stabaal::Simulator<> same_extent{4};
    same_extent.Add();
    CHECK(sim.MemoryUsage() == same_extent.MemoryUsage());
    CHECK_FALSE(sim.MeasureZ(added));
    CHECK_FALSE(sim.IsDeterministicMeasureZ(0));
    CHECK_THROWS_AS(sim.H(1), std::out_of_range);
}
