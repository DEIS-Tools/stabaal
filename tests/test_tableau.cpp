#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <array>
#include <complex>
#include <random>
#include <set>
#include <stdexcept>
#include <string>

#include "tableau.h"

// ---------------------------------------------------------------------------
// Test helpers
// ---------------------------------------------------------------------------

// Friend of Tableau: gives the tests access to private state and helpers.
struct TableauTestAccess {
    static basis_table& comp(Tableau& t) { return t.computational; }
    static basis_table& had(Tableau& t) { return t.hadamard; }
    static uint64_t qubits(Tableau& t) { return t.qubits; }

    static void clean_up(Tableau& t, int pivot, Generator& gen) { t.clean_up(pivot, gen); }
    static void collapse(Tableau& t, int pivot, bool p) { t.collapse(pivot, p); }
    static bool random_measurement(Tableau& t, int q, Tableau::Basis b) { return t.random_measurement(q, b); }
    static int e(Tableau& t, Generator a, Generator b) { return t.e(a, b); }
    static bool phase(Tableau& t, Generator a, Generator b, int rule) { return t.phase(a, b, rule); }
    static void swap_bit(Tableau& t, bit_ref a, bit_ref b) { t.swap_bit(a, b); }
};
using A = TableauTestAccess;

// Builds a generator from a Pauli string, e.g. pauli("XIZY", true) = -X⊗I⊗Z⊗Y.
// Character i is the Pauli acting on qubit i.
static Generator pauli(const std::string& s, bool neg = false) {
    Generator g(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        switch (s[i]) {
            case 'I': break;
            case 'X': g.x[i] = true; break;
            case 'Z': g.z[i] = true; break;
            case 'Y': g.x[i] = true; g.z[i] = true; break;
            default: throw std::invalid_argument("bad pauli char");
        }
    }
    g.r = neg;
    return g;
}

// Renders a generator as a signed Pauli string, e.g. "-XIZY".
static std::string str(const Generator& g) {
    std::string s = g.r ? "-" : "+";
    auto n = std::max(g.x.size(), g.z.size());
    for (size_t i = 0; i < n; ++i) {
        bool x = i < g.x.size() && g.x[i];
        bool z = i < g.z.size() && g.z[i];
        s += x ? (z ? 'Y' : 'X') : (z ? 'Z' : 'I');
    }
    return s;
}

static bool anticommute(const Generator& a, const Generator& b) {
    return ((a.x & b.z).count() + (a.z & b.x).count()) % 2 == 1;
}

// All rows have the right size, and the rows form a valid symplectic basis:
// comp[i] and had[j] anticommute iff i == j, everything else commutes.
static void check_symplectic(Tableau& t) {
    auto n = A::qubits(t);
    auto& comp = A::comp(t);
    auto& had  = A::had(t);
    REQUIRE(comp.size() == n);
    REQUIRE(had.size() == n);
    for (size_t i = 0; i < n; ++i) {
        REQUIRE(comp[i].x.size() == n);
        REQUIRE(comp[i].z.size() == n);
        REQUIRE(had[i].x.size() == n);
        REQUIRE(had[i].z.size() == n);
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            CHECK(anticommute(comp[i], had[j]) == (i == j));
            CHECK_FALSE(anticommute(comp[i], comp[j]));
            CHECK_FALSE(anticommute(had[i], had[j]));
        }
    }
}

// Runs `f(seed)` for a range of seeds and collects the boolean outcomes.
// The constructor seeds rand() with time(), so `f` must call srand(seed)
// after constructing its tableau.
template <class F>
static std::set<bool> outcomes(F f, int seeds = 64) {
    std::set<bool> seen;
    for (int s = 0; s < seeds; ++s) seen.insert(f(s));
    return seen;
}

// Redirects std::cout for the lifetime of the object.
struct CoutCapture {
    std::ostringstream ss;
    std::streambuf* old;
    CoutCapture() : old(std::cout.rdbuf(ss.rdbuf())) {}
    ~CoutCapture() { std::cout.rdbuf(old); }
};

// 2x2 Pauli matrices used as an independent reference for e().
using Mat = std::array<std::complex<double>, 4>;
static Mat mat(char p) {
    using c = std::complex<double>;
    switch (p) {
        case 'X': return {c(0), c(1), c(1), c(0)};
        case 'Y': return {c(0), c(0, -1), c(0, 1), c(0)};
        case 'Z': return {c(1), c(0), c(0), c(-1)};
        default:  return {c(1), c(0), c(0), c(1)};
    }
}
static Mat mul(const Mat& a, const Mat& b) {
    return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]};
}

// ===========================================================================
// Public API
// ===========================================================================

TEST_CASE("Constructor initialises |0...0> tableau", "[tableau][ctor]") {
    auto n = GENERATE(1, 2, 5);
    Tableau t(n);

    REQUIRE(A::qubits(t) == static_cast<uint64_t>(n));
    check_symplectic(t);

    for (int i = 0; i < n; ++i) {
        std::string zi(n, 'I'), xi(n, 'I');
        zi[i] = 'Z';
        xi[i] = 'X';
        CHECK(str(A::comp(t)[i]) == "+" + zi);
        CHECK(str(A::had(t)[i]) == "+" + xi);
    }
}

TEST_CASE("Constructor with zero qubits produces empty tableau", "[tableau][ctor]") {
    Tableau t(0);
    CHECK(A::qubits(t) == 0);
    CHECK(A::comp(t).empty());
    CHECK(A::had(t).empty());
}

TEST_CASE("measure_z on fresh state is deterministically 0", "[tableau][measure]") {
    auto n = GENERATE(1, 3);
    Tableau t(n);
    for (int q = 0; q < n; ++q) {
        CHECK(t.measure_z(q) == false);
        CHECK(t.measure_z(q) == false);
    }
    // Deterministic measurement must not change the tableau.
    Tableau fresh(n);
    for (int q = 0; q < n; ++q) {
        CHECK(str(A::comp(t)[q]) == str(A::comp(fresh)[q]));
        CHECK(str(A::had(t)[q]) == str(A::had(fresh)[q]));
    }
}

TEST_CASE("measure_x on fresh state is random", "[tableau][measure]") {
    auto seen = outcomes([](int s) {
        Tableau t(1);
        srand(s);
        return t.measure_x(0);
    });
    CHECK(seen.size() == 2);
}

TEST_CASE("Repeated measurement returns the same outcome", "[tableau][measure]") {
    for (int s = 0; s < 32; ++s) {
        SECTION("X basis, seed " + std::to_string(s)) {
            Tableau t(2);
            srand(s);
            bool first = t.measure_x(0);
            CHECK(t.measure_x(0) == first);
            CHECK(t.measure_x(0) == first);
            check_symplectic(t);
        }
        SECTION("Z basis after H, seed " + std::to_string(s)) {
            Tableau t(2);
            t.apply_h(1);
            srand(s);
            bool first = t.measure_z(1);
            CHECK(t.measure_z(1) == first);
            check_symplectic(t);
        }
    }
}

TEST_CASE("Repeated measurement is consistent when the measured row has a negative sign", "[tableau][measure]") {
    for (int s = 0; s < 32; ++s) {
        SECTION("Z basis, row -X, seed " + std::to_string(s)) {
            Tableau t(1);
            t.apply_z(0);
            t.apply_h(0); // comp row = -X
            REQUIRE(str(A::comp(t)[0]) == "-X");
            srand(s);
            bool first = t.measure_z(0);
            CHECK(t.measure_z(0) == first);
        }
        SECTION("X basis, row -X, seed " + std::to_string(s)) {
            Tableau t(1);
            t.apply_z(0); // had row = -X
            REQUIRE(str(A::had(t)[0]) == "-X");
            srand(s);
            bool first = t.measure_x(0);
            CHECK(t.measure_x(0) == first);
        }
    }
}

TEST_CASE("Measuring X then Z re-randomises Z", "[tableau][measure]") {
    auto seen = outcomes([](int s) {
        Tableau t(1);
        srand(s);
        t.measure_x(0);
        bool z = t.measure_z(0);
        CHECK(t.measure_z(0) == z);
        return z;
    });
    CHECK(seen.size() == 2);
}

TEST_CASE("Measuring one qubit does not disturb an independent qubit", "[tableau][measure]") {
    for (int s = 0; s < 16; ++s) {
        Tableau t(3);
        t.apply_x(1);
        t.apply_h(0);
        t.apply_h(2);
        t.apply_z(2); // qubit 2 in |->
        srand(s);
        t.measure_z(0);
        CHECK(t.measure_z(1) == true);
        CHECK(t.measure_x(2) == true);
        check_symplectic(t);
    }
}

TEST_CASE("apply_x", "[tableau][gate][x]") {
    SECTION("flips |0> to |1>") {
        Tableau t(1);
        t.apply_x(0);
        CHECK(t.measure_z(0) == true);
        CHECK(str(A::comp(t)[0]) == "-Z");
        CHECK(str(A::had(t)[0]) == "+X");
    }
    SECTION("X twice is identity") {
        Tableau t(1);
        t.apply_x(0);
        t.apply_x(0);
        CHECK(t.measure_z(0) == false);
        CHECK(str(A::comp(t)[0]) == "+Z");
    }
    SECTION("X leaves |+> unchanged") {
        Tableau t(1);
        t.apply_h(0);
        t.apply_x(0);
        CHECK(t.measure_x(0) == false);
    }
    SECTION("X turns |-> into -|->, X measurement still 1") {
        Tableau t(1);
        t.apply_x(0);
        t.apply_h(0);
        t.apply_x(0);
        CHECK(t.measure_x(0) == true);
    }
    SECTION("only the targeted qubit flips") {
        Tableau t(4);
        t.apply_x(2);
        CHECK(t.measure_z(0) == false);
        CHECK(t.measure_z(1) == false);
        CHECK(t.measure_z(2) == true);
        CHECK(t.measure_z(3) == false);
        check_symplectic(t);
    }
}

TEST_CASE("apply_z", "[tableau][gate][z]") {
    SECTION("Z on |0> has no observable effect") {
        Tableau t(1);
        t.apply_z(0);
        CHECK(t.measure_z(0) == false);
        CHECK(str(A::comp(t)[0]) == "+Z");
        CHECK(str(A::had(t)[0]) == "-X");
    }
    SECTION("Z maps |+> to |->") {
        Tableau t(1);
        t.apply_h(0);
        t.apply_z(0);
        CHECK(t.measure_x(0) == true);
    }
    SECTION("Z twice is identity") {
        Tableau t(1);
        t.apply_h(0);
        t.apply_z(0);
        t.apply_z(0);
        CHECK(t.measure_x(0) == false);
        CHECK(str(A::had(t)[0]) == "+Z");
    }
    SECTION("only the targeted qubit is affected") {
        Tableau t(3);
        for (int q = 0; q < 3; ++q) t.apply_h(q);
        t.apply_z(1);
        CHECK(t.measure_x(0) == false);
        CHECK(t.measure_x(1) == true);
        CHECK(t.measure_x(2) == false);
    }
}

TEST_CASE("apply_h", "[tableau][gate][h]") {
    SECTION("H|0> = |+>: X deterministic 0") {
        Tableau t(1);
        t.apply_h(0);
        CHECK(str(A::comp(t)[0]) == "+X");
        CHECK(str(A::had(t)[0]) == "+Z");
        CHECK(t.measure_x(0) == false);
    }
    SECTION("H|0> measured in Z is random") {
        auto seen = outcomes([](int s) {
            Tableau t(1);
            t.apply_h(0);
            srand(s);
            return t.measure_z(0);
        });
        CHECK(seen.size() == 2);
    }
    SECTION("H twice restores rows exactly") {
        Tableau t(2);
        t.apply_x(0);
        t.apply_z(0);
        auto c = str(A::comp(t)[0]), h = str(A::had(t)[0]);
        t.apply_h(0);
        t.apply_h(0);
        CHECK(str(A::comp(t)[0]) == c);
        CHECK(str(A::had(t)[0]) == h);
    }
    SECTION("H|1> = |->") {
        Tableau t(1);
        t.apply_x(0);
        t.apply_h(0);
        CHECK(t.measure_x(0) == true);
    }
    SECTION("HZH = X") {
        Tableau a(1), b(1);
        a.apply_h(0);
        a.apply_z(0);
        a.apply_h(0);
        b.apply_x(0);
        CHECK(str(A::comp(a)[0]) == str(A::comp(b)[0]));
        CHECK(str(A::had(a)[0]) == str(A::had(b)[0]));
        CHECK(a.measure_z(0) == true);
    }
    SECTION("H on one qubit leaves others untouched") {
        Tableau t(3);
        t.apply_h(1);
        CHECK(str(A::comp(t)[0]) == "+ZII");
        CHECK(str(A::comp(t)[1]) == "+IXI");
        CHECK(str(A::had(t)[1]) == "+IZI");
        CHECK(str(A::comp(t)[2]) == "+IIZ");
        check_symplectic(t);
    }
}

TEST_CASE("apply_s", "[tableau][gate][s]") {
    SECTION("S on |0> leaves Z measurement 0, X row becomes -Y") {
        Tableau t(1);
        t.apply_s(0);
        CHECK(t.measure_z(0) == false);
        CHECK(str(A::comp(t)[0]) == "+Z");
        CHECK(str(A::had(t)[0]) == "-Y");
        check_symplectic(t);
    }
    SECTION("SH|0> = |+i>: X row is +Y, X measurement random") {
        Tableau t(1);
        t.apply_h(0);
        t.apply_s(0);
        CHECK(str(A::had(t)[0]) == "+Y");
        auto seen = outcomes([](int s) {
            Tableau u(1);
            u.apply_h(0);
            u.apply_s(0);
            srand(s);
            return u.measure_x(0);
        });
        CHECK(seen.size() == 2);
    }
    SECTION("S twice equals Z: |+> -> |->") {
        Tableau t(1);
        t.apply_h(0);
        t.apply_s(0);
        t.apply_s(0);
        CHECK(str(A::had(t)[0]) == "-Z");
        CHECK(t.measure_x(0) == true);
    }
    SECTION("S twice matches apply_z on rows") {
        Tableau a(1), b(1);
        a.apply_s(0);
        a.apply_s(0);
        b.apply_z(0);
        CHECK(str(A::had(a)[0]) == str(A::had(b)[0]));
        CHECK(str(A::comp(a)[0]) == str(A::comp(b)[0]));
    }
    SECTION("S cycles X row: X -> -Y -> -X -> +Y -> X") {
        Tableau t(1);
        const char* expected[] = {"-Y", "-X", "+Y", "+X"};
        for (auto exp : expected) {
            t.apply_s(0);
            CHECK(str(A::had(t)[0]) == exp);
            CHECK(str(A::comp(t)[0]) == "+Z");
        }
    }
    SECTION("S on different qubits is independent") {
        Tableau t(2);
        t.apply_h(0);
        t.apply_h(1);
        t.apply_s(0);
        t.apply_s(0);
        CHECK(t.measure_x(0) == true);
        CHECK(t.measure_x(1) == false);
        check_symplectic(t);
    }
    SECTION("measuring |+i> in X twice is consistent") {
        for (int s = 0; s < 16; ++s) {
            Tableau t(1);
            t.apply_h(0);
            t.apply_s(0);
            srand(s);
            bool first = t.measure_x(0);
            CHECK(t.measure_x(0) == first);
        }
    }
}

TEST_CASE("apply_cx truth table on basis states", "[tableau][gate][cx]") {
    auto a = GENERATE(false, true);
    auto b = GENERATE(false, true);
    CAPTURE(a, b);
    Tableau t(2);
    if (a) t.apply_x(0);
    if (b) t.apply_x(1);
    t.apply_cx(0, 1);
    CHECK(t.measure_z(0) == a);
    CHECK(t.measure_z(1) == (a != b));
    check_symplectic(t);
}

TEST_CASE("apply_cx", "[tableau][gate][cx]") {
    SECTION("exact rows after H0 CX01") {
        Tableau t(2);
        t.apply_h(0);
        t.apply_cx(0, 1);
        CHECK(str(A::comp(t)[0]) == "+XI");
        CHECK(str(A::comp(t)[1]) == "+XZ");
        CHECK(str(A::had(t)[0]) == "+ZX");
        CHECK(str(A::had(t)[1]) == "+IX");
        check_symplectic(t);
    }
    SECTION("CX twice is identity") {
        Tableau t(2);
        t.apply_x(0);
        t.apply_h(1);
        t.apply_z(1);
        std::vector<std::string> before;
        for (int i = 0; i < 2; ++i) {
            before.push_back(str(A::comp(t)[i]));
            before.push_back(str(A::had(t)[i]));
        }
        t.apply_cx(0, 1);
        t.apply_cx(0, 1);
        for (int i = 0; i < 2; ++i) {
            CHECK(str(A::comp(t)[i]) == before[2 * i]);
            CHECK(str(A::had(t)[i]) == before[2 * i + 1]);
        }
    }
    SECTION("phase kickback: |+>|-> -> |->|->") {
        Tableau t(2);
        t.apply_h(0);
        t.apply_x(1);
        t.apply_h(1);
        t.apply_cx(0, 1);
        CHECK(t.measure_x(0) == true);
        CHECK(t.measure_x(1) == true);
    }
    SECTION("|+>|0> target |+> is unaffected: CX on |+>|+> is identity") {
        Tableau t(2);
        t.apply_h(0);
        t.apply_h(1);
        t.apply_cx(0, 1);
        CHECK(t.measure_x(0) == false);
        CHECK(t.measure_x(1) == false);
    }
    SECTION("reversed control/target (c > t)") {
        Tableau t(2);
        t.apply_x(1);
        t.apply_cx(1, 0);
        CHECK(t.measure_z(0) == true);
        CHECK(t.measure_z(1) == true);
    }
    SECTION("non-adjacent qubits leave the middle untouched") {
        Tableau t(3);
        t.apply_x(2);
        t.apply_cx(2, 0);
        CHECK(t.measure_z(0) == true);
        CHECK(t.measure_z(1) == false);
        CHECK(t.measure_z(2) == true);
        check_symplectic(t);
    }
}

TEST_CASE("Bell state correlations", "[tableau][gate][cx][entangle]") {
    SECTION("Z outcomes always equal and both occur") {
        auto seen = outcomes([](int s) {
            Tableau t(2);
            t.apply_h(0);
            t.apply_cx(0, 1);
            srand(s);
            bool a = t.measure_z(0);
            bool b = t.measure_z(1);
            CHECK(a == b);
            return a;
        });
        CHECK(seen.size() == 2);
    }
    SECTION("measuring the target first also correlates") {
        for (int s = 0; s < 32; ++s) {
            Tableau t(2);
            t.apply_h(0);
            t.apply_cx(0, 1);
            srand(s);
            bool b = t.measure_z(1);
            CHECK(t.measure_z(0) == b);
        }
    }
    SECTION("X outcomes always equal and both occur") {
        auto seen = outcomes([](int s) {
            Tableau t(2);
            t.apply_h(0);
            t.apply_cx(0, 1);
            srand(s);
            bool a = t.measure_x(0);
            bool b = t.measure_x(1);
            CHECK(a == b);
            return a;
        });
        CHECK(seen.size() == 2);
    }
    SECTION("|Phi-> (Z on control) gives anti-correlated X outcomes") {
        for (int s = 0; s < 32; ++s) {
            Tableau t(2);
            t.apply_h(0);
            t.apply_z(0);
            t.apply_cx(0, 1);
            srand(s);
            bool a = t.measure_x(0);
            CHECK(t.measure_x(1) != a);
        }
    }
    SECTION("|Psi+> (X on target) gives anti-correlated Z outcomes") {
        for (int s = 0; s < 32; ++s) {
            Tableau t(2);
            t.apply_h(0);
            t.apply_cx(0, 1);
            t.apply_x(1);
            srand(s);
            bool a = t.measure_z(0);
            CHECK(t.measure_z(1) != a);
        }
    }
}

TEST_CASE("GHZ state correlations", "[tableau][gate][cx][entangle]") {
    auto seen = outcomes([](int s) {
        Tableau t(3);
        t.apply_h(0);
        t.apply_cx(0, 1);
        t.apply_cx(0, 2);
        check_symplectic(t);
        srand(s);
        bool a = t.measure_z(0);
        CHECK(t.measure_z(1) == a);
        CHECK(t.measure_z(2) == a);
        return a;
    });
    CHECK(seen.size() == 2);
}

// ---------------------------------------------------------------------------
// Generated gates: find_gate_effect, apply<M> and the GATE_LIST methods
// ---------------------------------------------------------------------------

// Renders a GenEffect as a signed Pauli string, same format as str(Generator).
template <size_t N>
static std::string str(const GenEffect<N>& g) {
    std::string s = g.r ? "-" : "+";
    for (size_t i = 0; i < N; ++i)
        s += g.x[i] ? (g.z[i] ? 'Y' : 'X') : (g.z[i] ? 'Z' : 'I');
    return s;
}

// Every row of a and b, signs included, is identical.
static void check_same_rows(Tableau& a, Tableau& b) {
    REQUIRE(A::qubits(a) == A::qubits(b));
    for (size_t i = 0; i < A::qubits(a); ++i) {
        CAPTURE(i);
        CHECK(str(A::comp(a)[i]) == str(A::comp(b)[i]));
        CHECK(str(A::had(a)[i])  == str(A::had(b)[i]));
    }
}

// Puts a and b in the same pseudo-random state using the hand-written gates.
static void scramble(Tableau& a, Tableau& b, unsigned seed, int steps = 20) {
    std::mt19937 rng(seed);
    auto n = int(A::qubits(a));
    for (int i = 0; i < steps; ++i) {
        int q = rng() % n;
        int t = (q + 1 + rng() % (n - 1)) % n;
        switch (rng() % 5) {
            case 0: a.apply_h(q);     b.apply_h(q);     break;
            case 1: a.apply_s(q);     b.apply_s(q);     break;
            case 2: a.apply_x(q);     b.apply_x(q);     break;
            case 3: a.apply_z(q);     b.apply_z(q);     break;
            case 4: a.apply_cx(q, t); b.apply_cx(q, t); break;
        }
    }
}

TEST_CASE("find_gate_effect runs at compile time", "[tableau][codegen][gate_effect]") {
    constexpr auto h = find_gate_effect<1>(H_mat);
    STATIC_REQUIRE(h.computational[0].x[0]);
    STATIC_REQUIRE(!h.computational[0].z[0]);
    STATIC_REQUIRE(h.hadamard[0].z[0]);
    STATIC_REQUIRE(!h.hadamard[0].x[0]);
}

// computational[j] is the image of Z_j and hadamard[j] the image of X_j, both as U^dagger P U.
TEST_CASE("find_gate_effect single-qubit images", "[tableau][codegen][gate_effect]") {
    SECTION("I") {
        constexpr auto g = find_gate_effect<1>(I_mat);
        CHECK(str(g.computational[0]) == "+Z");
        CHECK(str(g.hadamard[0])      == "+X");
    }
    SECTION("H") {
        constexpr auto g = find_gate_effect<1>(H_mat);
        CHECK(str(g.computational[0]) == "+X");
        CHECK(str(g.hadamard[0])      == "+Z");
    }
    SECTION("S: X -> -Y") {
        constexpr auto g = find_gate_effect<1>(S_mat);
        CHECK(str(g.computational[0]) == "+Z");
        CHECK(str(g.hadamard[0])      == "-Y");
    }
    SECTION("S^dagger: X -> +Y") {
        constexpr auto g = find_gate_effect<1>(S_mat.adjoint());
        CHECK(str(g.computational[0]) == "+Z");
        CHECK(str(g.hadamard[0])      == "+Y");
    }
    SECTION("X") {
        constexpr auto g = find_gate_effect<1>(X_mat);
        CHECK(str(g.computational[0]) == "-Z");
        CHECK(str(g.hadamard[0])      == "+X");
    }
    SECTION("Z") {
        constexpr auto g = find_gate_effect<1>(Z_mat);
        CHECK(str(g.computational[0]) == "+Z");
        CHECK(str(g.hadamard[0])      == "-X");
    }
    SECTION("Y") {
        constexpr auto g = find_gate_effect<1>(Y_mat);
        CHECK(str(g.computational[0]) == "-Z");
        CHECK(str(g.hadamard[0])      == "-X");
    }
}

TEST_CASE("find_gate_effect CX images", "[tableau][codegen][gate_effect]") {
    constexpr auto g = find_gate_effect<2>(CX_mat);
    CHECK(str(g.computational[0]) == "+ZI");
    CHECK(str(g.computational[1]) == "+ZZ");
    CHECK(str(g.hadamard[0])      == "+XX");
    CHECK(str(g.hadamard[1])      == "+IX");
}

TEST_CASE("generated methods on a fresh tableau", "[tableau][codegen][gate]") {
    SECTION("h(0) gives |+>") {
        Tableau t(1);
        t.h(0);
        CHECK(str(A::comp(t)[0]) == "+X");
        CHECK(str(A::had(t)[0])  == "+Z");
        CHECK(t.measure_x(0) == false);
    }
    SECTION("s(0) X row becomes -Y") {
        Tableau t(1);
        t.s(0);
        CHECK(str(A::comp(t)[0]) == "+Z");
        CHECK(str(A::had(t)[0])  == "-Y");
    }
    SECTION("h(0) z(0) gives |->") {
        Tableau t(1);
        t.h(0);
        t.z(0);
        CHECK(t.measure_x(0) == true);
    }
    SECTION("exact rows after h(0) cx(0, 1)") {
        Tableau t(2);
        t.h(0);
        t.cx(0, 1);
        CHECK(str(A::comp(t)[0]) == "+XI");
        CHECK(str(A::comp(t)[1]) == "+XZ");
        CHECK(str(A::had(t)[0])  == "+ZX");
        CHECK(str(A::had(t)[1])  == "+IX");
        check_symplectic(t);
    }
}

TEST_CASE("generated methods match the hand-written gates", "[tableau][codegen][gate]") {
    auto seed = GENERATE(range(0u, 32u));
    CAPTURE(seed);
    Tableau a(4), b(4);
    scramble(a, b, seed);
    check_same_rows(a, b);

    SECTION("h") { for (int q = 0; q < 4; ++q) { a.apply_h(q); b.h(q); } }
    SECTION("s") { for (int q = 0; q < 4; ++q) { a.apply_s(q); b.s(q); } }
    SECTION("z") { for (int q = 0; q < 4; ++q) { a.apply_z(q); b.z(q); } }
    SECTION("cx adjacent")     { a.apply_cx(0, 1); b.cx(0, 1); }
    SECTION("cx reversed")     { a.apply_cx(3, 1); b.cx(3, 1); }
    SECTION("cx non-adjacent") { a.apply_cx(0, 3); b.cx(0, 3); }

    check_same_rows(a, b);
    check_symplectic(b);
}

TEST_CASE("apply<M> for gates without a generated method", "[tableau][codegen][gate]") {
    auto seed = GENERATE(range(0u, 32u));
    CAPTURE(seed);
    Tableau a(3), b(3);
    scramble(a, b, seed);

    SECTION("I leaves the tableau unchanged") {
        b.apply<I_mat>(1);
    }
    SECTION("X matches apply_x") {
        a.apply_x(2);
        b.apply<X_mat>(2);
    }
    SECTION("Y matches Z then X (equal up to global phase)") {
        a.apply_z(0);
        a.apply_x(0);
        b.apply<Y_mat>(0);
    }
    SECTION("S^dagger matches three S") {
        for (int i = 0; i < 3; ++i) a.apply_s(1);
        b.apply<S_mat.adjoint()>(1);
    }

    check_same_rows(a, b);
    check_symplectic(b);
}

TEST_CASE("generated methods match the hand-written gates on random circuits", "[tableau][codegen][gate]") {
    auto seed = GENERATE(range(0u, 64u));
    CAPTURE(seed);
    std::mt19937 rng(seed);
    Tableau a(5), b(5);
    for (int i = 0; i < 40; ++i) {
        int q = rng() % 5;
        int t = (q + 1 + rng() % 4) % 5;
        switch (rng() % 4) {
            case 0: a.apply_h(q);     b.h(q);     break;
            case 1: a.apply_s(q);     b.s(q);     break;
            case 2: a.apply_z(q);     b.z(q);     break;
            case 3: a.apply_cx(q, t); b.cx(q, t); break;
        }
    }
    check_same_rows(a, b);
    check_symplectic(b);
}

TEST_CASE("add appends a fresh |0> qubit", "[tableau][add]") {
    SECTION("all rows resized and new rows are Z_n / X_n") {
        Tableau t(1);
        t.add();
        REQUIRE(A::qubits(t) == 2);
        check_symplectic(t);
        CHECK(str(A::comp(t)[0]) == "+ZI");
        CHECK(str(A::had(t)[0]) == "+XI");
        CHECK(str(A::comp(t)[1]) == "+IZ");
        CHECK(str(A::had(t)[1]) == "+IX");
    }
    SECTION("existing state is preserved") {
        Tableau t(1);
        t.apply_x(0);
        t.add();
        check_symplectic(t);
        CHECK(t.measure_z(0) == true);
        CHECK(t.measure_z(1) == false);
    }
    SECTION("gates act on the new qubit") {
        Tableau t(1);
        t.add();
        check_symplectic(t);
        t.apply_x(1);
        t.apply_cx(1, 0);
        CHECK(t.measure_z(0) == true);
        CHECK(t.measure_z(1) == true);
    }
    SECTION("repeated add from zero qubits") {
        Tableau t(0);
        for (int n = 1; n <= 4; ++n) {
            t.add();
            REQUIRE(A::qubits(t) == static_cast<uint64_t>(n));
            check_symplectic(t);
        }
    }
    SECTION("add after entangling keeps correlations") {
        Tableau t(2);
        t.apply_h(0);
        t.apply_cx(0, 1);
        t.add();
        check_symplectic(t);
        srand(3);
        bool a = t.measure_z(0);
        CHECK(t.measure_z(1) == a);
        CHECK(t.measure_z(2) == false);
    }
}

TEST_CASE("print", "[tableau][print]") {
    SECTION("fresh 1-qubit tableau") {
        Tableau t(1);
        CoutCapture cap;
        t.print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 1 | 0 | 0 |"
              "\nZ0: | 0 | 1 | 0 |"
              "\n");
    }
    SECTION("fresh 2-qubit tableau") {
        Tableau t(2);
        CoutCapture cap;
        t.print();
        CHECK(cap.ss.str() ==
              "\n    | x x | z z | r |"
              "\nX0: | 1 0 | 0 0 | 0 |"
              "\nX1: | 0 1 | 0 0 | 0 |"
              "\nZ0: | 0 0 | 1 0 | 0 |"
              "\nZ1: | 0 0 | 0 1 | 0 |"
              "\n");
    }
    SECTION("sign column after X") {
        Tableau t(1);
        t.apply_x(0);
        CoutCapture cap;
        t.print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 1 | 0 | 0 |"
              "\nZ0: | 0 | 1 | 1 |"
              "\n");
    }
    SECTION("rows swapped after H") {
        Tableau t(1);
        t.apply_h(0);
        CoutCapture cap;
        t.print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 0 | 1 | 0 |"
              "\nZ0: | 1 | 0 | 0 |"
              "\n");
    }
    SECTION("print does not modify the tableau") {
        Tableau t(2);
        t.apply_h(0);
        auto before = str(A::comp(t)[0]) + str(A::had(t)[0]);
        {
            CoutCapture cap;
            t.print();
        }
        CHECK(str(A::comp(t)[0]) + str(A::had(t)[0]) == before);
    }
}

// ===========================================================================
// Private helpers
// ===========================================================================

TEST_CASE("swap_bit", "[tableau][private][swap_bit]") {
    Tableau t(1);
    SECTION("differing bits in two bitsets") {
        Generator g = pauli("XI");
        A::swap_bit(t, g.x[0], g.z[0]);
        CHECK(str(g) == "+ZI");
        A::swap_bit(t, g.x[0], g.z[0]);
        CHECK(str(g) == "+XI");
    }
    SECTION("both bits 0") {
        Generator g = pauli("IZ");
        A::swap_bit(t, g.x[0], g.z[0]);
        CHECK(str(g) == "+IZ");
    }
    SECTION("both bits 1") {
        Generator g = pauli("YZ");
        A::swap_bit(t, g.x[0], g.z[0]);
        CHECK(str(g) == "+YZ");
    }
    SECTION("positions within the same bitset") {
        boost::dynamic_bitset<> b(4);
        b[1] = true;
        A::swap_bit(t, b[1], b[3]);
        CHECK(b[0] == false);
        CHECK(b[1] == false);
        CHECK(b[2] == false);
        CHECK(b[3] == true);
    }
    SECTION("same position is a no-op") {
        boost::dynamic_bitset<> b(2);
        b[0] = true;
        A::swap_bit(t, b[0], b[0]);
        CHECK(b[0] == true);
        CHECK(b[1] == false);
    }
    SECTION("does not touch the sign or other bits") {
        Generator g = pauli("XZY", true);
        A::swap_bit(t, g.x[1], g.z[1]);
        CHECK(str(g) == "-XXY");
    }
}

TEST_CASE("e: single-qubit table", "[tableau][private][e]") {
    Tableau t(1);
    struct Case { const char* a; const char* b; int e; };
    auto c = GENERATE(values<Case>({
        {"I", "I", 0}, {"I", "X", 0}, {"I", "Y", 0}, {"I", "Z", 0},
        {"X", "I", 0}, {"X", "X", 0}, {"X", "Y", 1}, {"X", "Z", -1},
        {"Y", "I", 0}, {"Y", "X", 3}, {"Y", "Y", 4}, {"Y", "Z", 1},
        {"Z", "I", 0}, {"Z", "X", 1}, {"Z", "Y", 3}, {"Z", "Z", 0},
    }));
    CAPTURE(c.a, c.b);
    CHECK(A::e(t, pauli(c.a), pauli(c.b)) == c.e);
}

TEST_CASE("e: matches Pauli matrix multiplication", "[tableau][private][e]") {
    Tableau t(1);
    const std::complex<double> I(0, 1);
    for (char a : std::string("IXYZ")) {
        for (char b : std::string("IXYZ")) {
            CAPTURE(a, b);
            Generator ga = pauli(std::string(1, a));
            Generator gb = pauli(std::string(1, b));
            Generator prod(1);
            prod.x = ga.x ^ gb.x;
            prod.z = ga.z ^ gb.z;
            char p = str(prod)[1];

            int e = ((A::e(t, ga, gb) % 4) + 4) % 4;
            std::complex<double> factor = std::pow(I, e);
            Mat lhs = mul(mat(a), mat(b));
            Mat rhs = mat(p);
            for (int k = 0; k < 4; ++k) {
                CHECK(std::abs(lhs[k] - factor * rhs[k]) < 1e-12);
            }
        }
    }
}

TEST_CASE("e: is additive over tensor products", "[tableau][private][e]") {
    Tableau t(1);
    auto pair = GENERATE(table<std::string, std::string>({
        {"XY", "ZZ"}, {"YYYY", "XZXZ"}, {"IXYZ", "ZYXI"}, {"ZZZ", "XXX"}, {"YIY", "IYI"},
    }));
    auto [a, b] = pair;
    int sum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += A::e(t, pauli(std::string(1, a[i])), pauli(std::string(1, b[i])));
    }
    CAPTURE(a, b);
    CHECK(A::e(t, pauli(a), pauli(b)) == sum);
}

TEST_CASE("e: ignores sign bits", "[tableau][private][e]") {
    Tableau t(1);
    CHECK(A::e(t, pauli("XZ", true), pauli("ZY", false)) == A::e(t, pauli("XZ"), pauli("ZY")));
    CHECK(A::e(t, pauli("Y", true), pauli("Y", true)) == 4);
}

TEST_CASE("phase: sign of product for all signs and rules", "[tableau][private][phase]") {
    Tableau t(1);
    auto ab = GENERATE(table<std::string, std::string>({
        {"X", "Z"}, {"Z", "X"}, {"Y", "Y"}, {"X", "X"}, {"Z", "Y"}, {"Y", "X"},
        {"XZ", "ZX"}, {"YI", "IY"}, {"XX", "ZZ"},
    }));
    auto ra   = GENERATE(false, true);
    auto rb   = GENERATE(false, true);
    auto rule = GENERATE(0, 1, 2, 3);
    auto [a, b] = ab;
    CAPTURE(a, b, ra, rb, rule);

    Generator ga = pauli(a, ra), gb = pauli(b, rb);
    int total = A::e(t, ga, gb) + 2 * ra + 2 * rb + rule;
    int mod   = ((total % 4) + 4) % 4;

    if (mod % 2 == 1) {
        CHECK_THROWS_AS(A::phase(t, ga, gb, rule), std::runtime_error);
    } else {
        CHECK(A::phase(t, ga, gb, rule) == (mod == 2));
    }
}

TEST_CASE("phase: hand-picked values", "[tableau][private][phase]") {
    Tableau t(1);
    // Z*Z = I, no sign
    CHECK(A::phase(t, pauli("Z"), pauli("Z"), 0) == false);
    // -Z * Z = -I
    CHECK(A::phase(t, pauli("Z", true), pauli("Z"), 0) == true);
    // -Z * -Z = I
    CHECK(A::phase(t, pauli("Z", true), pauli("Z", true), 0) == false);
    // Y*Y = I (e = 4)
    CHECK(A::phase(t, pauli("Y"), pauli("Y"), 0) == false);
    // S rule: X*Z with rule 3 => e = -1 + 3 = 2 => negative
    CHECK(A::phase(t, pauli("X"), pauli("Z"), 3) == true);
    // Z*X with rule 3 => 1 + 3 = 4 => positive
    CHECK(A::phase(t, pauli("Z"), pauli("X"), 3) == false);
    // (XZ)(ZX) = (-iY)(iY) = +YY; bits are YY so e should be 0
    CHECK(A::phase(t, pauli("XZ"), pauli("ZX"), 0) == false);
    // (XX)(ZZ) = (-iY)(-iY) = -YY
    CHECK(A::phase(t, pauli("XX"), pauli("ZZ"), 0) == true);
}

TEST_CASE("phase: throws on odd total", "[tableau][private][phase]") {
    Tableau t(1);
    CHECK_THROWS_AS(A::phase(t, pauli("X"), pauli("Z"), 0), std::runtime_error);
    CHECK_THROWS_AS(A::phase(t, pauli("Z"), pauli("X"), 0), std::runtime_error);
    CHECK_THROWS_AS(A::phase(t, pauli("Z"), pauli("Z"), 1), std::runtime_error);
    CHECK_THROWS_AS(A::phase(t, pauli("I"), pauli("I"), 3), std::runtime_error);
    CHECK_NOTHROW(A::phase(t, pauli("X"), pauli("Z"), 1));
}

TEST_CASE("collapse", "[tableau][private][collapse]") {
    SECTION("after H, collapse(0, p) gives (-1)^p Z and X") {
        auto p = GENERATE(false, true);
        Tableau t(1);
        t.apply_h(0); // comp = X, had = Z
        A::collapse(t, 0, p);
        CHECK(str(A::comp(t)[0]) == (p ? "-Z" : "+Z"));
        CHECK(str(A::had(t)[0]) == "+X");
        CHECK(t.measure_z(0) == p);
        check_symplectic(t);
    }
    SECTION("Y column picks up H sign and X^p sign") {
        auto p = GENERATE(false, true);
        Tableau t(1);
        A::comp(t)[0] = pauli("Y");
        A::had(t)[0]  = pauli("Z");
        A::collapse(t, 0, p);
        // H Y H = -Y, then X(-Y)X = +Y
        CHECK(str(A::comp(t)[0]) == (p ? "+Y" : "-Y"));
        CHECK(str(A::had(t)[0]) == "+X");
    }
    SECTION("only the pivot column is changed") {
        Tableau t(3);
        A::comp(t)[0] = pauli("XZY");
        A::comp(t)[1] = pauli("ZXI", true);
        A::comp(t)[2] = pauli("IYZ");
        A::had(t)[0]  = pauli("ZIX");
        A::had(t)[1]  = pauli("IZZ");
        A::had(t)[2]  = pauli("YXI", true);
        A::collapse(t, 1, false);
        CHECK(str(A::comp(t)[0]) == "+XXY");
        CHECK(str(A::comp(t)[1]) == "-ZZI");
        CHECK(str(A::comp(t)[2]) == "-IYZ");
        CHECK(str(A::had(t)[0]) == "+ZIX");
        CHECK(str(A::had(t)[1]) == "+IXZ");
        CHECK(str(A::had(t)[2]) == "-YZI");
    }
    SECTION("p = 1 flips sign of rows with Z on pivot after swap") {
        Tableau t(2);
        A::comp(t)[0] = pauli("IX");
        A::comp(t)[1] = pauli("IZ");
        A::had(t)[0]  = pauli("XI");
        A::had(t)[1]  = pauli("ZI");
        A::collapse(t, 1, true);
        CHECK(str(A::comp(t)[0]) == "-IZ");
        CHECK(str(A::comp(t)[1]) == "+IX");
        CHECK(str(A::had(t)[0]) == "+XI");
        CHECK(str(A::had(t)[1]) == "+ZI");
    }
}

TEST_CASE("clean_up", "[tableau][private][clean_up]") {
    SECTION("physical state H0 H1 CX01: Z1 row X0X1 reduced to X0") {
        Tableau t(2);
        A::comp(t)[0] = pauli("XI");
        A::comp(t)[1] = pauli("XX");
        A::had(t)[0]  = pauli("ZZ");
        A::had(t)[1]  = pauli("IZ");
        check_symplectic(t);

        A::clean_up(t, 0, A::comp(t)[1]);
        CHECK(str(A::comp(t)[0]) == "+XX");
        CHECK(str(A::comp(t)[1]) == "+XI");
        CHECK(str(A::had(t)[0]) == "+IZ");
        CHECK(str(A::had(t)[1]) == "+ZZ");
        check_symplectic(t);
    }
    SECTION("column CNOT rule including sign updates") {
        Tableau t(2);
        A::comp(t)[0] = pauli("XZ");
        A::comp(t)[1] = pauli("XX");
        A::had(t)[0]  = pauli("YZ", true);
        A::had(t)[1]  = pauli("IY");
        A::clean_up(t, 0, A::comp(t)[1]);
        CHECK(str(A::comp(t)[0]) == "-YY"); // X0Z1 -> -Y0Y1
        CHECK(str(A::comp(t)[1]) == "+XI");
        CHECK(str(A::had(t)[0]) == "-XY");  // Y0Z1 -> X0Y1, sign kept
        CHECK(str(A::had(t)[1]) == "+ZY");  // Y1 -> Z0Y1
    }
    SECTION("only columns where gen has X are touched") {
        Tableau t(3);
        A::comp(t)[0] = pauli("XZZ");
        A::comp(t)[1] = pauli("XIX");
        A::comp(t)[2] = pauli("IXI");
        A::had(t)[0]  = pauli("ZII");
        A::had(t)[1]  = pauli("IZI");
        A::had(t)[2]  = pauli("IIZ");
        A::clean_up(t, 0, A::comp(t)[1]);
        CHECK(str(A::comp(t)[0]) == "-YZY");
        CHECK(str(A::comp(t)[1]) == "+XII");
        CHECK(str(A::comp(t)[2]) == "+IXI");
        CHECK(str(A::had(t)[0]) == "+ZII");
        CHECK(str(A::had(t)[1]) == "+IZI");
        CHECK(str(A::had(t)[2]) == "+ZIZ");
    }
    SECTION("non-zero pivot leaves earlier columns alone") {
        Tableau t(3);
        A::comp(t)[0] = pauli("IXX");
        A::comp(t)[1] = pauli("XIZ");
        A::comp(t)[2] = pauli("ZII");
        A::had(t)[0]  = pauli("III");
        A::had(t)[1]  = pauli("III");
        A::had(t)[2]  = pauli("III");
        A::clean_up(t, 1, A::comp(t)[0]);
        CHECK(str(A::comp(t)[0]) == "+IXI");
        CHECK(str(A::comp(t)[1]) == "+XZZ");
        CHECK(str(A::comp(t)[2]) == "+ZII");
    }
    SECTION("single X on pivot with no Z is a no-op") {
        Tableau t(2);
        t.apply_h(0);
        std::vector<std::string> before;
        for (int i = 0; i < 2; ++i) before.push_back(str(A::comp(t)[i]) + str(A::had(t)[i]));
        A::clean_up(t, 0, A::comp(t)[0]);
        for (int i = 0; i < 2; ++i) CHECK(str(A::comp(t)[i]) + str(A::had(t)[i]) == before[i]);
    }
    SECTION("Y on pivot is rotated to X in the table (S-dagger on column)") {
        Tableau t(1);
        A::comp(t)[0] = pauli("Y");
        A::had(t)[0]  = pauli("Z");
        A::clean_up(t, 0, A::comp(t)[0]);
        CHECK(str(A::comp(t)[0]) == "+X");
        CHECK(str(A::had(t)[0]) == "+Z");
    }
    SECTION("Z created on pivot by column ops is also removed") {
        Tableau t(2);
        A::comp(t)[1] = pauli("XY"); // after CNOT column op becomes +YZ
        A::comp(t)[0] = pauli("ZI");
        A::had(t)[0]  = pauli("XX");
        A::had(t)[1]  = pauli("IZ");
        A::clean_up(t, 0, A::comp(t)[1]);
        CHECK(A::comp(t)[1].x.count() == 1);
        CHECK(A::comp(t)[1].x[0] == true);
        CHECK(A::comp(t)[1].z[0] == false);
    }
}

TEST_CASE("random_measurement", "[tableau][private][random_measurement]") {
    SECTION("computational basis on |+>") {
        auto seen = outcomes([](int s) {
            Tableau t(1);
            t.apply_h(0);
            srand(s);
            bool p = A::random_measurement(t, 0, Tableau::COMPUTATIONAL);
            CHECK(A::comp(t)[0].x.none());
            CHECK(A::comp(t)[0].r == p);
            CHECK(t.measure_z(0) == p);
            check_symplectic(t);
            return p;
        });
        CHECK(seen.size() == 2);
    }
    SECTION("hadamard basis on |0>") {
        auto seen = outcomes([](int s) {
            Tableau t(1);
            srand(s);
            bool p = A::random_measurement(t, 0, Tableau::HADAMARD);
            CHECK(A::had(t)[0].x.none());
            CHECK(A::had(t)[0].r == p);
            CHECK(t.measure_x(0) == p);
            check_symplectic(t);
            return p;
        });
        CHECK(seen.size() == 2);
    }
    SECTION("multi-X row: pivot is first X, other qubit untouched") {
        for (int s = 0; s < 16; ++s) {
            // State |++> written as H0 H1 CX01: Z1 row = X0X1.
            Tableau t(2);
            A::comp(t)[0] = pauli("XI");
            A::comp(t)[1] = pauli("XX");
            A::had(t)[0]  = pauli("ZZ");
            A::had(t)[1]  = pauli("IZ");
            srand(s);
            bool p = A::random_measurement(t, 1, Tableau::COMPUTATIONAL);
            CHECK(A::comp(t)[1].x.none());
            CHECK(A::comp(t)[1].r == p);
            CHECK(t.measure_z(1) == p);
            // Qubit 0 is still |+>.
            CHECK(t.measure_x(0) == false);
            check_symplectic(t);
        }
    }
    SECTION("Y row on measured qubit") {
        for (int s = 0; s < 16; ++s) {
            Tableau t(1);
            t.apply_h(0);
            t.apply_s(0); // had row = +Y
            srand(s);
            bool p = A::random_measurement(t, 0, Tableau::HADAMARD);
            CHECK(A::had(t)[0].x.none());
            CHECK(A::had(t)[0].r == p);
            check_symplectic(t);
        }
    }
    SECTION("negative-sign row: returned outcome matches the collapsed row") {
        auto seen = outcomes([](int s) {
            Tableau t(1);
            t.apply_z(0);
            t.apply_h(0); // comp row = -X
            srand(s);
            bool p = A::random_measurement(t, 0, Tableau::COMPUTATIONAL);
            CHECK(A::comp(t)[0].x.none());
            CHECK(A::comp(t)[0].r == p);
            return p;
        });
        CHECK(seen.size() == 2);
    }
    SECTION("same seed gives same outcome") {
        Tableau a(1), b(1);
        srand(7);
        bool pa = A::random_measurement(a, 0, Tableau::HADAMARD);
        srand(7);
        bool pb = A::random_measurement(b, 0, Tableau::HADAMARD);
        CHECK(pa == pb);
    }
}
