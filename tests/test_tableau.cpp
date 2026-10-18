#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>

#include "tableau.hpp"
#include "types.hpp"

namespace stabaal {

// ---------------------------------------------------------------------------
// Test helpers
// ---------------------------------------------------------------------------

// Friend of Tableau: gives the tests access to private state and helpers.
struct TableauTestAccess {
    static void Seed(std::uint32_t value) { Tableau::RandomEngine().seed(value); }
    static BasisTable& Comp(Tableau& t) { return t.computational_; }
    static BasisTable& Had(Tableau& t) { return t.hadamard_; }
    static QubitIdx Qubits(Tableau& t) { return t.num_qubits_; }

    static void CleanUp(Tableau& t, QubitIdx pivot, const PauliString& gen) {
        t.CleanUp(pivot, gen);
    }
    static void Collapse(Tableau& t, QubitIdx pivot, bool p) { t.Collapse(pivot, p); }
    static bool RandomMeasurement(Tableau& t, QubitIdx q, Tableau::Basis b) {
        return t.RandomMeasurement(q, b);
    }
    static std::int64_t E(Tableau& t, PauliString a, PauliString b) { return t.E(a, b); }
    static bool Phase(Tableau& t, PauliString a, PauliString b, std::uint8_t rule) {
        return t.Phase(a, b, rule);
    }
    static void SwapBit(Tableau& t, BitRef a, BitRef b) { t.SwapBit(a, b); }
};
using A = TableauTestAccess;

// Builds a generator from a Pauli string, e.g. Pauli("XIZY", true) = -X⊗I⊗Z⊗Y.
// Character i is the Pauli acting on qubit i.
static PauliString Pauli(const std::string& s, bool neg = false) {
    PauliString g(static_cast<QubitIdx>(s.size()));
    for (std::size_t i = 0; i < s.size(); ++i) {
        switch (s[i]) {
            case 'I':
                break;
            case 'X':
                g.x_[i] = true;
                break;
            case 'Z':
                g.z_[i] = true;
                break;
            case 'Y':
                g.x_[i] = true;
                g.z_[i] = true;
                break;
            default:
                throw std::invalid_argument("bad pauli char");
        }
    }
    g.sign_ = neg;
    return g;
}

// Renders a generator as a signed Pauli string, e.g. "-XIZY".
static std::string Str(const PauliString& g) {
    std::string s = g.sign_ ? "-" : "+";
    auto n = std::max(g.x_.size(), g.z_.size());
    for (std::size_t i = 0; i < n; ++i) {
        bool x = i < g.x_.size() && g.x_[i];
        bool z = i < g.z_.size() && g.z_[i];
        s += x ? (z ? 'Y' : 'X') : (z ? 'Z' : 'I');
    }
    return s;
}

static bool Anticommute(const PauliString& a, const PauliString& b) {
    return ((a.x_ & b.z_).count() + (a.z_ & b.x_).count()) % 2 == 1;
}

// All rows have the right size, and the rows form a valid symplectic basis:
// comp[i] and had[j] anticommute iff i == j, everything else commutes.
static void CheckSymplectic(Tableau& t) {
    auto n = A::Qubits(t);
    auto& comp = A::Comp(t);
    auto& had = A::Had(t);
    REQUIRE(comp.size() == n);
    REQUIRE(had.size() == n);
    for (std::size_t i = 0; i < n; ++i) {
        REQUIRE(comp[i].x_.size() == n);
        REQUIRE(comp[i].z_.size() == n);
        REQUIRE(had[i].x_.size() == n);
        REQUIRE(had[i].z_.size() == n);
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK(Anticommute(comp[i], had[j]) == (i == j));
            CHECK_FALSE(Anticommute(comp[i], comp[j]));
            CHECK_FALSE(Anticommute(had[i], had[j]));
        }
    }
}

// Runs `f(seed)` for a range of seeds and collects the boolean outcomes.
// Each callback seeds the C++ engine through the test helper.
template <class F>
static std::set<bool> Outcomes(F f, std::uint32_t seeds = 64) {
    std::set<bool> seen;
    for (std::uint32_t s = 0; s < seeds; ++s) seen.insert(f(s));
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
static Mat PauliMatrix(char p) {
    using Complex = std::complex<double>;
    switch (p) {
        case 'X':
            return {Complex(0), Complex(1), Complex(1), Complex(0)};
        case 'Y':
            return {Complex(0), Complex(0, -1), Complex(0, 1), Complex(0)};
        case 'Z':
            return {Complex(1), Complex(0), Complex(0), Complex(-1)};
        default:
            return {Complex(1), Complex(0), Complex(0), Complex(1)};
    }
}
static Mat Multiply(const Mat& a, const Mat& b) {
    return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3], a[2] * b[0] + a[3] * b[2],
            a[2] * b[1] + a[3] * b[3]};
}

// ===========================================================================
// Public API
// ===========================================================================

TEST_CASE("Constructor initialises |0...0> tableau", "[tableau][ctor]") {
    QubitIdx n = GENERATE(1, 2, 5);
    Tableau t(n);

    REQUIRE(A::Qubits(t) == n);
    CheckSymplectic(t);

    for (std::size_t i = 0; i < n; ++i) {
        std::string zi(n, 'I'), xi(n, 'I');
        zi[i] = 'Z';
        xi[i] = 'X';
        CHECK(Str(A::Comp(t)[i]) == "+" + zi);
        CHECK(Str(A::Had(t)[i]) == "+" + xi);
    }
}

TEST_CASE("Constructor with zero qubits produces empty tableau", "[tableau][ctor]") {
    Tableau t(0);
    CHECK(A::Qubits(t) == 0);
    CHECK(A::Comp(t).empty());
    CHECK(A::Had(t).empty());
}

TEST_CASE("MeasureZ on fresh state is deterministically 0", "[tableau][measure]") {
    QubitIdx n = GENERATE(1, 3);
    Tableau t(n);
    for (QubitIdx q = 0; q < n; ++q) {
        CHECK(t.MeasureZ(q) == false);
        CHECK(t.MeasureZ(q) == false);
    }
    // Deterministic measurement must not change the tableau.
    Tableau fresh(n);
    for (QubitIdx q = 0; q < n; ++q) {
        CHECK(Str(A::Comp(t)[q]) == Str(A::Comp(fresh)[q]));
        CHECK(Str(A::Had(t)[q]) == Str(A::Had(fresh)[q]));
    }
}

TEST_CASE("MeasureX on fresh state is random", "[tableau][measure]") {
    auto seen = Outcomes([](std::uint32_t s) {
        Tableau t(1);
        A::Seed(s);
        return t.MeasureX(0);
    });
    CHECK(seen.size() == 2);
}

TEST_CASE("Repeated measurement returns the same outcome", "[tableau][measure]") {
    for (std::uint32_t s = 0; s < 32; ++s) {
        SECTION("X basis, seed " + std::to_string(s)) {
            Tableau t(2);
            A::Seed(s);
            bool first = t.MeasureX(0);
            CHECK(t.MeasureX(0) == first);
            CHECK(t.MeasureX(0) == first);
            CheckSymplectic(t);
        }
        SECTION("Z basis after H, seed " + std::to_string(s)) {
            Tableau t(2);
            t.ApplyH(1);
            A::Seed(s);
            bool first = t.MeasureZ(1);
            CHECK(t.MeasureZ(1) == first);
            CheckSymplectic(t);
        }
    }
}

TEST_CASE("Repeated measurement is consistent when the measured row has a negative sign",
          "[tableau][measure]") {
    for (std::uint32_t s = 0; s < 32; ++s) {
        SECTION("Z basis, row -X, seed " + std::to_string(s)) {
            Tableau t(1);
            t.ApplyZ(0);
            t.ApplyH(0);  // comp row = -X
            REQUIRE(Str(A::Comp(t)[0]) == "-X");
            A::Seed(s);
            bool first = t.MeasureZ(0);
            CHECK(t.MeasureZ(0) == first);
        }
        SECTION("X basis, row -X, seed " + std::to_string(s)) {
            Tableau t(1);
            t.ApplyZ(0);  // had row = -X
            REQUIRE(Str(A::Had(t)[0]) == "-X");
            A::Seed(s);
            bool first = t.MeasureX(0);
            CHECK(t.MeasureX(0) == first);
        }
    }
}

TEST_CASE("Measuring X then Z re-randomises Z", "[tableau][measure]") {
    auto seen = Outcomes([](std::uint32_t s) {
        Tableau t(1);
        A::Seed(s);
        t.MeasureX(0);
        bool z = t.MeasureZ(0);
        CHECK(t.MeasureZ(0) == z);
        return z;
    });
    CHECK(seen.size() == 2);
}

TEST_CASE("Measuring one qubit does not disturb an independent qubit", "[tableau][measure]") {
    for (std::uint32_t s = 0; s < 16; ++s) {
        Tableau t(3);
        t.ApplyX(1);
        t.ApplyH(0);
        t.ApplyH(2);
        t.ApplyZ(2);  // qubit 2 in |->
        A::Seed(s);
        t.MeasureZ(0);
        CHECK(t.MeasureZ(1) == true);
        CHECK(t.MeasureX(2) == true);
        CheckSymplectic(t);
    }
}

TEST_CASE("ApplyX", "[tableau][gate][x]") {
    SECTION("flips |0> to |1>") {
        Tableau t(1);
        t.ApplyX(0);
        CHECK(t.MeasureZ(0) == true);
        CHECK(Str(A::Comp(t)[0]) == "-Z");
        CHECK(Str(A::Had(t)[0]) == "+X");
    }
    SECTION("X twice is identity") {
        Tableau t(1);
        t.ApplyX(0);
        t.ApplyX(0);
        CHECK(t.MeasureZ(0) == false);
        CHECK(Str(A::Comp(t)[0]) == "+Z");
    }
    SECTION("X leaves |+> unchanged") {
        Tableau t(1);
        t.ApplyH(0);
        t.ApplyX(0);
        CHECK(t.MeasureX(0) == false);
    }
    SECTION("X turns |-> into -|->, X measurement still 1") {
        Tableau t(1);
        t.ApplyX(0);
        t.ApplyH(0);
        t.ApplyX(0);
        CHECK(t.MeasureX(0) == true);
    }
    SECTION("only the targeted qubit flips") {
        Tableau t(4);
        t.ApplyX(2);
        CHECK(t.MeasureZ(0) == false);
        CHECK(t.MeasureZ(1) == false);
        CHECK(t.MeasureZ(2) == true);
        CHECK(t.MeasureZ(3) == false);
        CheckSymplectic(t);
    }
}

TEST_CASE("ApplyZ", "[tableau][gate][z]") {
    SECTION("Z on |0> has no observable effect") {
        Tableau t(1);
        t.ApplyZ(0);
        CHECK(t.MeasureZ(0) == false);
        CHECK(Str(A::Comp(t)[0]) == "+Z");
        CHECK(Str(A::Had(t)[0]) == "-X");
    }
    SECTION("Z maps |+> to |->") {
        Tableau t(1);
        t.ApplyH(0);
        t.ApplyZ(0);
        CHECK(t.MeasureX(0) == true);
    }
    SECTION("Z twice is identity") {
        Tableau t(1);
        t.ApplyH(0);
        t.ApplyZ(0);
        t.ApplyZ(0);
        CHECK(t.MeasureX(0) == false);
        CHECK(Str(A::Had(t)[0]) == "+Z");
    }
    SECTION("only the targeted qubit is affected") {
        Tableau t(3);
        for (QubitIdx q = 0; q < 3; ++q) t.ApplyH(q);
        t.ApplyZ(1);
        CHECK(t.MeasureX(0) == false);
        CHECK(t.MeasureX(1) == true);
        CHECK(t.MeasureX(2) == false);
    }
}

TEST_CASE("ApplyH", "[tableau][gate][h]") {
    SECTION("H|0> = |+>: X deterministic 0") {
        Tableau t(1);
        t.ApplyH(0);
        CHECK(Str(A::Comp(t)[0]) == "+X");
        CHECK(Str(A::Had(t)[0]) == "+Z");
        CHECK(t.MeasureX(0) == false);
    }
    SECTION("H|0> measured in Z is random") {
        auto seen = Outcomes([](std::uint32_t s) {
            Tableau t(1);
            t.ApplyH(0);
            A::Seed(s);
            return t.MeasureZ(0);
        });
        CHECK(seen.size() == 2);
    }
    SECTION("H twice restores rows exactly") {
        Tableau t(2);
        t.ApplyX(0);
        t.ApplyZ(0);
        auto c = Str(A::Comp(t)[0]), h = Str(A::Had(t)[0]);
        t.ApplyH(0);
        t.ApplyH(0);
        CHECK(Str(A::Comp(t)[0]) == c);
        CHECK(Str(A::Had(t)[0]) == h);
    }
    SECTION("H|1> = |->") {
        Tableau t(1);
        t.ApplyX(0);
        t.ApplyH(0);
        CHECK(t.MeasureX(0) == true);
    }
    SECTION("HZH = X") {
        Tableau a(1), b(1);
        a.ApplyH(0);
        a.ApplyZ(0);
        a.ApplyH(0);
        b.ApplyX(0);
        CHECK(Str(A::Comp(a)[0]) == Str(A::Comp(b)[0]));
        CHECK(Str(A::Had(a)[0]) == Str(A::Had(b)[0]));
        CHECK(a.MeasureZ(0) == true);
    }
    SECTION("H on one qubit leaves others untouched") {
        Tableau t(3);
        t.ApplyH(1);
        CHECK(Str(A::Comp(t)[0]) == "+ZII");
        CHECK(Str(A::Comp(t)[1]) == "+IXI");
        CHECK(Str(A::Had(t)[1]) == "+IZI");
        CHECK(Str(A::Comp(t)[2]) == "+IIZ");
        CheckSymplectic(t);
    }
}

TEST_CASE("ApplyS", "[tableau][gate][s]") {
    SECTION("S on |0> leaves Z measurement 0, X row becomes -Y") {
        Tableau t(1);
        t.ApplyS(0);
        CHECK(t.MeasureZ(0) == false);
        CHECK(Str(A::Comp(t)[0]) == "+Z");
        CHECK(Str(A::Had(t)[0]) == "-Y");
        CheckSymplectic(t);
    }
    SECTION("SH|0> = |+i>: X row is +Y, X measurement random") {
        Tableau t(1);
        t.ApplyH(0);
        t.ApplyS(0);
        CHECK(Str(A::Had(t)[0]) == "+Y");
        auto seen = Outcomes([](std::uint32_t s) {
            Tableau u(1);
            u.ApplyH(0);
            u.ApplyS(0);
            A::Seed(s);
            return u.MeasureX(0);
        });
        CHECK(seen.size() == 2);
    }
    SECTION("S twice equals Z: |+> -> |->") {
        Tableau t(1);
        t.ApplyH(0);
        t.ApplyS(0);
        t.ApplyS(0);
        CHECK(Str(A::Had(t)[0]) == "-Z");
        CHECK(t.MeasureX(0) == true);
    }
    SECTION("S twice matches ApplyZ on rows") {
        Tableau a(1), b(1);
        a.ApplyS(0);
        a.ApplyS(0);
        b.ApplyZ(0);
        CHECK(Str(A::Had(a)[0]) == Str(A::Had(b)[0]));
        CHECK(Str(A::Comp(a)[0]) == Str(A::Comp(b)[0]));
    }
    SECTION("S cycles X row: X -> -Y -> -X -> +Y -> X") {
        Tableau t(1);
        const char* expected[] = {"-Y", "-X", "+Y", "+X"};
        for (auto exp : expected) {
            t.ApplyS(0);
            CHECK(Str(A::Had(t)[0]) == exp);
            CHECK(Str(A::Comp(t)[0]) == "+Z");
        }
    }
    SECTION("S on different qubits is independent") {
        Tableau t(2);
        t.ApplyH(0);
        t.ApplyH(1);
        t.ApplyS(0);
        t.ApplyS(0);
        CHECK(t.MeasureX(0) == true);
        CHECK(t.MeasureX(1) == false);
        CheckSymplectic(t);
    }
    SECTION("measuring |+i> in X twice is consistent") {
        for (std::uint32_t s = 0; s < 16; ++s) {
            Tableau t(1);
            t.ApplyH(0);
            t.ApplyS(0);
            A::Seed(s);
            bool first = t.MeasureX(0);
            CHECK(t.MeasureX(0) == first);
        }
    }
}

TEST_CASE("ApplyCx truth table on basis states", "[tableau][gate][cx]") {
    auto a = GENERATE(false, true);
    auto b = GENERATE(false, true);
    CAPTURE(a, b);
    Tableau t(2);
    if (a) t.ApplyX(0);
    if (b) t.ApplyX(1);
    t.ApplyCx(0, 1);
    CHECK(t.MeasureZ(0) == a);
    CHECK(t.MeasureZ(1) == (a != b));
    CheckSymplectic(t);
}

TEST_CASE("ApplyCx", "[tableau][gate][cx]") {
    SECTION("exact rows after H0 CX01") {
        Tableau t(2);
        t.ApplyH(0);
        t.ApplyCx(0, 1);
        CHECK(Str(A::Comp(t)[0]) == "+XI");
        CHECK(Str(A::Comp(t)[1]) == "+XZ");
        CHECK(Str(A::Had(t)[0]) == "+ZX");
        CHECK(Str(A::Had(t)[1]) == "+IX");
        CheckSymplectic(t);
    }
    SECTION("CX twice is identity") {
        Tableau t(2);
        t.ApplyX(0);
        t.ApplyH(1);
        t.ApplyZ(1);
        std::vector<std::string> before;
        for (std::size_t i = 0; i < 2; ++i) {
            before.push_back(Str(A::Comp(t)[i]));
            before.push_back(Str(A::Had(t)[i]));
        }
        t.ApplyCx(0, 1);
        t.ApplyCx(0, 1);
        for (std::size_t i = 0; i < 2; ++i) {
            CHECK(Str(A::Comp(t)[i]) == before[2 * i]);
            CHECK(Str(A::Had(t)[i]) == before[2 * i + 1]);
        }
    }
    SECTION("phase kickback: |+>|-> -> |->|->") {
        Tableau t(2);
        t.ApplyH(0);
        t.ApplyX(1);
        t.ApplyH(1);
        t.ApplyCx(0, 1);
        CHECK(t.MeasureX(0) == true);
        CHECK(t.MeasureX(1) == true);
    }
    SECTION("|+>|0> target |+> is unaffected: CX on |+>|+> is identity") {
        Tableau t(2);
        t.ApplyH(0);
        t.ApplyH(1);
        t.ApplyCx(0, 1);
        CHECK(t.MeasureX(0) == false);
        CHECK(t.MeasureX(1) == false);
    }
    SECTION("reversed control/target (c > t)") {
        Tableau t(2);
        t.ApplyX(1);
        t.ApplyCx(1, 0);
        CHECK(t.MeasureZ(0) == true);
        CHECK(t.MeasureZ(1) == true);
    }
    SECTION("non-adjacent qubits leave the middle untouched") {
        Tableau t(3);
        t.ApplyX(2);
        t.ApplyCx(2, 0);
        CHECK(t.MeasureZ(0) == true);
        CHECK(t.MeasureZ(1) == false);
        CHECK(t.MeasureZ(2) == true);
        CheckSymplectic(t);
    }
}

TEST_CASE("Bell state correlations", "[tableau][gate][cx][entangle]") {
    SECTION("Z outcomes always equal and both occur") {
        auto seen = Outcomes([](std::uint32_t s) {
            Tableau t(2);
            t.ApplyH(0);
            t.ApplyCx(0, 1);
            A::Seed(s);
            bool a = t.MeasureZ(0);
            bool b = t.MeasureZ(1);
            CHECK(a == b);
            return a;
        });
        CHECK(seen.size() == 2);
    }
    SECTION("measuring the target first also correlates") {
        for (std::uint32_t s = 0; s < 32; ++s) {
            Tableau t(2);
            t.ApplyH(0);
            t.ApplyCx(0, 1);
            A::Seed(s);
            bool b = t.MeasureZ(1);
            CHECK(t.MeasureZ(0) == b);
        }
    }
    SECTION("X outcomes always equal and both occur") {
        auto seen = Outcomes([](std::uint32_t s) {
            Tableau t(2);
            t.ApplyH(0);
            t.ApplyCx(0, 1);
            A::Seed(s);
            bool a = t.MeasureX(0);
            bool b = t.MeasureX(1);
            CHECK(a == b);
            return a;
        });
        CHECK(seen.size() == 2);
    }
    SECTION("|Phi-> (Z on control) gives anti-correlated X outcomes") {
        for (std::uint32_t s = 0; s < 32; ++s) {
            Tableau t(2);
            t.ApplyH(0);
            t.ApplyZ(0);
            t.ApplyCx(0, 1);
            A::Seed(s);
            bool a = t.MeasureX(0);
            CHECK(t.MeasureX(1) != a);
        }
    }
    SECTION("|Psi+> (X on target) gives anti-correlated Z outcomes") {
        for (std::uint32_t s = 0; s < 32; ++s) {
            Tableau t(2);
            t.ApplyH(0);
            t.ApplyCx(0, 1);
            t.ApplyX(1);
            A::Seed(s);
            bool a = t.MeasureZ(0);
            CHECK(t.MeasureZ(1) != a);
        }
    }
}

TEST_CASE("GHZ state correlations", "[tableau][gate][cx][entangle]") {
    auto seen = Outcomes([](std::uint32_t s) {
        Tableau t(3);
        t.ApplyH(0);
        t.ApplyCx(0, 1);
        t.ApplyCx(0, 2);
        CheckSymplectic(t);
        A::Seed(s);
        bool a = t.MeasureZ(0);
        CHECK(t.MeasureZ(1) == a);
        CHECK(t.MeasureZ(2) == a);
        return a;
    });
    CHECK(seen.size() == 2);
}

TEST_CASE("Add appends a fresh |0> qubit", "[tableau][Add]") {
    SECTION("all rows resized and new rows are Z_n / X_n") {
        Tableau t(1);
        t.Add();
        REQUIRE(A::Qubits(t) == 2);
        CheckSymplectic(t);
        CHECK(Str(A::Comp(t)[0]) == "+ZI");
        CHECK(Str(A::Had(t)[0]) == "+XI");
        CHECK(Str(A::Comp(t)[1]) == "+IZ");
        CHECK(Str(A::Had(t)[1]) == "+IX");
    }
    SECTION("existing state is preserved") {
        Tableau t(1);
        t.ApplyX(0);
        t.Add();
        CheckSymplectic(t);
        CHECK(t.MeasureZ(0) == true);
        CHECK(t.MeasureZ(1) == false);
    }
    SECTION("gates act on the new qubit") {
        Tableau t(1);
        t.Add();
        CheckSymplectic(t);
        t.ApplyX(1);
        t.ApplyCx(1, 0);
        CHECK(t.MeasureZ(0) == true);
        CHECK(t.MeasureZ(1) == true);
    }
    SECTION("repeated Add from zero qubits") {
        Tableau t(0);
        for (QubitIdx n = 1; n <= 4; ++n) {
            t.Add();
            REQUIRE(A::Qubits(t) == n);
            CheckSymplectic(t);
        }
    }
    SECTION("Add after entangling keeps correlations") {
        Tableau t(2);
        t.ApplyH(0);
        t.ApplyCx(0, 1);
        t.Add();
        CheckSymplectic(t);
        A::Seed(3);
        bool a = t.MeasureZ(0);
        CHECK(t.MeasureZ(1) == a);
        CHECK(t.MeasureZ(2) == false);
    }
}

TEST_CASE("Print", "[tableau][Print]") {
    SECTION("fresh 1-qubit tableau") {
        const Tableau t(1);
        CoutCapture cap;
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 1 | 0 | 0 |"
              "\nZ0: | 0 | 1 | 0 |"
              "\n");
    }
    SECTION("fresh 2-qubit tableau") {
        Tableau t(2);
        CoutCapture cap;
        t.Print();
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
        t.ApplyX(0);
        CoutCapture cap;
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 1 | 0 | 0 |"
              "\nZ0: | 0 | 1 | 1 |"
              "\n");
    }
    SECTION("rows swapped after H") {
        Tableau t(1);
        t.ApplyH(0);
        CoutCapture cap;
        t.Print();
        CHECK(cap.ss.str() ==
              "\n    | x | z | r |"
              "\nX0: | 0 | 1 | 0 |"
              "\nZ0: | 1 | 0 | 0 |"
              "\n");
    }
    SECTION("Print does not modify the tableau") {
        Tableau t(2);
        t.ApplyH(0);
        auto before = Str(A::Comp(t)[0]) + Str(A::Had(t)[0]);
        {
            CoutCapture cap;
            t.Print();
        }
        CHECK(Str(A::Comp(t)[0]) + Str(A::Had(t)[0]) == before);
    }
}

// ===========================================================================
// Private helpers
// ===========================================================================

TEST_CASE("SwapBit", "[tableau][private][SwapBit]") {
    Tableau t(1);
    SECTION("differing bits in two bitsets") {
        PauliString g = Pauli("XI");
        A::SwapBit(t, g.x_[0], g.z_[0]);
        CHECK(Str(g) == "+ZI");
        A::SwapBit(t, g.x_[0], g.z_[0]);
        CHECK(Str(g) == "+XI");
    }
    SECTION("both bits 0") {
        PauliString g = Pauli("IZ");
        A::SwapBit(t, g.x_[0], g.z_[0]);
        CHECK(Str(g) == "+IZ");
    }
    SECTION("both bits 1") {
        PauliString g = Pauli("YZ");
        A::SwapBit(t, g.x_[0], g.z_[0]);
        CHECK(Str(g) == "+YZ");
    }
    SECTION("positions within the same bitset") {
        boost::dynamic_bitset<> b(4);
        b[1] = true;
        A::SwapBit(t, b[1], b[3]);
        CHECK(b[0] == false);
        CHECK(b[1] == false);
        CHECK(b[2] == false);
        CHECK(b[3] == true);
    }
    SECTION("same position is a no-op") {
        boost::dynamic_bitset<> b(2);
        b[0] = true;
        A::SwapBit(t, b[0], b[0]);
        CHECK(b[0] == true);
        CHECK(b[1] == false);
    }
    SECTION("does not touch the sign or other bits") {
        PauliString g = Pauli("XZY", true);
        A::SwapBit(t, g.x_[1], g.z_[1]);
        CHECK(Str(g) == "-XXY");
    }
}

TEST_CASE("e: single-qubit table", "[tableau][private][e]") {
    Tableau t(1);
    struct Case {
        const char* a;
        const char* b;
        std::int64_t e;
    };
    auto c = GENERATE(values<Case>({
        {"I", "I", 0},
        {"I", "X", 0},
        {"I", "Y", 0},
        {"I", "Z", 0},
        {"X", "I", 0},
        {"X", "X", 0},
        {"X", "Y", 1},
        {"X", "Z", -1},
        {"Y", "I", 0},
        {"Y", "X", 3},
        {"Y", "Y", 4},
        {"Y", "Z", 1},
        {"Z", "I", 0},
        {"Z", "X", 1},
        {"Z", "Y", 3},
        {"Z", "Z", 0},
    }));
    CAPTURE(c.a, c.b);
    CHECK(A::E(t, Pauli(c.a), Pauli(c.b)) == c.e);
}

TEST_CASE("e: matches Pauli matrix multiplication", "[tableau][private][e]") {
    Tableau t(1);
    const std::complex<double> imaginary_unit(0, 1);
    for (char a : std::string("IXYZ")) {
        for (char b : std::string("IXYZ")) {
            CAPTURE(a, b);
            PauliString ga = Pauli(std::string(1, a));
            PauliString gb = Pauli(std::string(1, b));
            PauliString prod(1);
            prod.x_ = ga.x_ ^ gb.x_;
            prod.z_ = ga.z_ ^ gb.z_;
            char p = Str(prod)[1];

            std::int64_t e = ((A::E(t, ga, gb) % 4) + 4) % 4;
            std::complex<double> factor = std::pow(imaginary_unit, e);
            Mat lhs = Multiply(PauliMatrix(a), PauliMatrix(b));
            Mat rhs = PauliMatrix(p);
            for (std::size_t k = 0; k < 4; ++k) {
                CHECK(std::abs(lhs[k] - factor * rhs[k]) < 1e-12);
            }
        }
    }
}

TEST_CASE("e: is additive over tensor products", "[tableau][private][e]") {
    Tableau t(1);
    auto pair = GENERATE(table<std::string, std::string>({
        {"XY", "ZZ"},
        {"YYYY", "XZXZ"},
        {"IXYZ", "ZYXI"},
        {"ZZZ", "XXX"},
        {"YIY", "IYI"},
    }));
    auto [a, b] = pair;
    std::int64_t sum = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += A::E(t, Pauli(std::string(1, a[i])), Pauli(std::string(1, b[i])));
    }
    CAPTURE(a, b);
    CHECK(A::E(t, Pauli(a), Pauli(b)) == sum);
}

TEST_CASE("e: ignores sign bits", "[tableau][private][e]") {
    Tableau t(1);
    CHECK(A::E(t, Pauli("XZ", true), Pauli("ZY", false)) == A::E(t, Pauli("XZ"), Pauli("ZY")));
    CHECK(A::E(t, Pauli("Y", true), Pauli("Y", true)) == 4);
}

TEST_CASE("phase: sign of product for all signs and valid rules", "[tableau][private][phase]") {
    Tableau t(1);
    auto ab = GENERATE(table<std::string, std::string>({
        {"X", "Z"},
        {"Z", "X"},
        {"Y", "Y"},
        {"X", "X"},
        {"Z", "Y"},
        {"Y", "X"},
        {"XZ", "ZX"},
        {"YI", "IY"},
        {"XX", "ZZ"},
    }));
    auto ra = GENERATE(false, true);
    auto rb = GENERATE(false, true);
    std::uint8_t rule = GENERATE(0, 1, 2, 3);
    auto [a, b] = ab;
    CAPTURE(a, b, ra, rb, rule);

    PauliString ga = Pauli(a, ra), gb = Pauli(b, rb);
    std::int64_t total = A::E(t, ga, gb) + 2 * ra + 2 * rb + rule;
    std::int64_t mod = ((total % 4) + 4) % 4;

    if (mod % 2 == 1) return;  // Odd phases violate Phase()'s precondition.
    CHECK(A::Phase(t, ga, gb, rule) == (mod == 2));
}

TEST_CASE("phase: hand-picked values", "[tableau][private][phase]") {
    Tableau t(1);
    // Z*Z = I, no sign
    CHECK(A::Phase(t, Pauli("Z"), Pauli("Z"), 0) == false);
    // -Z * Z = -I
    CHECK(A::Phase(t, Pauli("Z", true), Pauli("Z"), 0) == true);
    // -Z * -Z = I
    CHECK(A::Phase(t, Pauli("Z", true), Pauli("Z", true), 0) == false);
    // Y*Y = I (e = 4)
    CHECK(A::Phase(t, Pauli("Y"), Pauli("Y"), 0) == false);
    // S rule: X*Z with rule 3 => e = -1 + 3 = 2 => negative
    CHECK(A::Phase(t, Pauli("X"), Pauli("Z"), 3) == true);
    // Z*X with rule 3 => 1 + 3 = 4 => positive
    CHECK(A::Phase(t, Pauli("Z"), Pauli("X"), 3) == false);
    // (XZ)(ZX) = (-iY)(iY) = +YY; bits are YY so e should be 0
    CHECK(A::Phase(t, Pauli("XZ"), Pauli("ZX"), 0) == false);
    // (XX)(ZZ) = (-iY)(-iY) = -YY
    CHECK(A::Phase(t, Pauli("XX"), Pauli("ZZ"), 0) == true);
}

TEST_CASE("Collapse", "[tableau][private][Collapse]") {
    SECTION("after H, Collapse(0, p) gives (-1)^p Z and X") {
        auto p = GENERATE(false, true);
        Tableau t(1);
        t.ApplyH(0);  // comp = X, had = Z
        A::Collapse(t, 0, p);
        CHECK(Str(A::Comp(t)[0]) == (p ? "-Z" : "+Z"));
        CHECK(Str(A::Had(t)[0]) == "+X");
        CHECK(t.MeasureZ(0) == p);
        CheckSymplectic(t);
    }
    SECTION("Y column picks up H sign and X^p sign") {
        auto p = GENERATE(false, true);
        Tableau t(1);
        A::Comp(t)[0] = Pauli("Y");
        A::Had(t)[0] = Pauli("Z");
        A::Collapse(t, 0, p);
        // H Y H = -Y, then X(-Y)X = +Y
        CHECK(Str(A::Comp(t)[0]) == (p ? "+Y" : "-Y"));
        CHECK(Str(A::Had(t)[0]) == "+X");
    }
    SECTION("only the pivot column is changed") {
        Tableau t(3);
        A::Comp(t)[0] = Pauli("XZY");
        A::Comp(t)[1] = Pauli("ZXI", true);
        A::Comp(t)[2] = Pauli("IYZ");
        A::Had(t)[0] = Pauli("ZIX");
        A::Had(t)[1] = Pauli("IZZ");
        A::Had(t)[2] = Pauli("YXI", true);
        A::Collapse(t, 1, false);
        CHECK(Str(A::Comp(t)[0]) == "+XXY");
        CHECK(Str(A::Comp(t)[1]) == "-ZZI");
        CHECK(Str(A::Comp(t)[2]) == "-IYZ");
        CHECK(Str(A::Had(t)[0]) == "+ZIX");
        CHECK(Str(A::Had(t)[1]) == "+IXZ");
        CHECK(Str(A::Had(t)[2]) == "-YZI");
    }
    SECTION("p = 1 flips sign of rows with Z on pivot after swap") {
        Tableau t(2);
        A::Comp(t)[0] = Pauli("IX");
        A::Comp(t)[1] = Pauli("IZ");
        A::Had(t)[0] = Pauli("XI");
        A::Had(t)[1] = Pauli("ZI");
        A::Collapse(t, 1, true);
        CHECK(Str(A::Comp(t)[0]) == "-IZ");
        CHECK(Str(A::Comp(t)[1]) == "+IX");
        CHECK(Str(A::Had(t)[0]) == "+XI");
        CHECK(Str(A::Had(t)[1]) == "+ZI");
    }
}

TEST_CASE("CleanUp", "[tableau][private][CleanUp]") {
    SECTION("physical state H0 H1 CX01: Z1 row X0X1 reduced to X0") {
        Tableau t(2);
        A::Comp(t)[0] = Pauli("XI");
        A::Comp(t)[1] = Pauli("XX");
        A::Had(t)[0] = Pauli("ZZ");
        A::Had(t)[1] = Pauli("IZ");
        CheckSymplectic(t);

        A::CleanUp(t, 0, A::Comp(t)[1]);
        CHECK(Str(A::Comp(t)[0]) == "+XX");
        CHECK(Str(A::Comp(t)[1]) == "+XI");
        CHECK(Str(A::Had(t)[0]) == "+IZ");
        CHECK(Str(A::Had(t)[1]) == "+ZZ");
        CheckSymplectic(t);
    }
    SECTION("column CNOT rule including sign updates") {
        Tableau t(2);
        A::Comp(t)[0] = Pauli("XZ");
        A::Comp(t)[1] = Pauli("XX");
        A::Had(t)[0] = Pauli("YZ", true);
        A::Had(t)[1] = Pauli("IY");
        A::CleanUp(t, 0, A::Comp(t)[1]);
        CHECK(Str(A::Comp(t)[0]) == "-YY");  // X0Z1 -> -Y0Y1
        CHECK(Str(A::Comp(t)[1]) == "+XI");
        CHECK(Str(A::Had(t)[0]) == "-XY");  // Y0Z1 -> X0Y1, sign kept
        CHECK(Str(A::Had(t)[1]) == "+ZY");  // Y1 -> Z0Y1
    }
    SECTION("only columns where gen has X are touched") {
        Tableau t(3);
        A::Comp(t)[0] = Pauli("XZZ");
        A::Comp(t)[1] = Pauli("XIX");
        A::Comp(t)[2] = Pauli("IXI");
        A::Had(t)[0] = Pauli("ZII");
        A::Had(t)[1] = Pauli("IZI");
        A::Had(t)[2] = Pauli("IIZ");
        A::CleanUp(t, 0, A::Comp(t)[1]);
        CHECK(Str(A::Comp(t)[0]) == "-YZY");
        CHECK(Str(A::Comp(t)[1]) == "+XII");
        CHECK(Str(A::Comp(t)[2]) == "+IXI");
        CHECK(Str(A::Had(t)[0]) == "+ZII");
        CHECK(Str(A::Had(t)[1]) == "+IZI");
        CHECK(Str(A::Had(t)[2]) == "+ZIZ");
    }
    SECTION("non-zero pivot leaves earlier columns alone") {
        Tableau t(3);
        A::Comp(t)[0] = Pauli("IXX");
        A::Comp(t)[1] = Pauli("XIZ");
        A::Comp(t)[2] = Pauli("ZII");
        A::Had(t)[0] = Pauli("III");
        A::Had(t)[1] = Pauli("III");
        A::Had(t)[2] = Pauli("III");
        A::CleanUp(t, 1, A::Comp(t)[0]);
        CHECK(Str(A::Comp(t)[0]) == "+IXI");
        CHECK(Str(A::Comp(t)[1]) == "+XZZ");
        CHECK(Str(A::Comp(t)[2]) == "+ZII");
    }
    SECTION("single X on pivot with no Z is a no-op") {
        Tableau t(2);
        t.ApplyH(0);
        std::vector<std::string> before;
        for (std::size_t i = 0; i < 2; ++i)
            before.push_back(Str(A::Comp(t)[i]) + Str(A::Had(t)[i]));
        A::CleanUp(t, 0, A::Comp(t)[0]);
        for (std::size_t i = 0; i < 2; ++i)
            CHECK(Str(A::Comp(t)[i]) + Str(A::Had(t)[i]) == before[i]);
    }
    SECTION("Y on pivot is rotated to X in the table (S-dagger on column)") {
        Tableau t(1);
        A::Comp(t)[0] = Pauli("Y");
        A::Had(t)[0] = Pauli("Z");
        A::CleanUp(t, 0, A::Comp(t)[0]);
        CHECK(Str(A::Comp(t)[0]) == "+X");
        CHECK(Str(A::Had(t)[0]) == "+Z");
    }
    SECTION("Z created on pivot by column ops is also removed") {
        Tableau t(2);
        A::Comp(t)[1] = Pauli("XY");  // after CNOT column op becomes +YZ
        A::Comp(t)[0] = Pauli("ZI");
        A::Had(t)[0] = Pauli("XX");
        A::Had(t)[1] = Pauli("IZ");
        A::CleanUp(t, 0, A::Comp(t)[1]);
        CHECK(A::Comp(t)[1].x_.count() == 1);
        CHECK(A::Comp(t)[1].x_[0] == true);
        CHECK(A::Comp(t)[1].z_[0] == false);
    }
}

TEST_CASE("RandomMeasurement", "[tableau][private][RandomMeasurement]") {
    SECTION("computational basis on |+>") {
        auto seen = Outcomes([](std::uint32_t s) {
            Tableau t(1);
            t.ApplyH(0);
            A::Seed(s);
            bool p = A::RandomMeasurement(t, 0, Tableau::Basis::kComputational);
            CHECK(A::Comp(t)[0].x_.none());
            CHECK(A::Comp(t)[0].sign_ == p);
            CHECK(t.MeasureZ(0) == p);
            CheckSymplectic(t);
            return p;
        });
        CHECK(seen.size() == 2);
    }
    SECTION("hadamard basis on |0>") {
        auto seen = Outcomes([](std::uint32_t s) {
            Tableau t(1);
            A::Seed(s);
            bool p = A::RandomMeasurement(t, 0, Tableau::Basis::kHadamard);
            CHECK(A::Had(t)[0].x_.none());
            CHECK(A::Had(t)[0].sign_ == p);
            CHECK(t.MeasureX(0) == p);
            CheckSymplectic(t);
            return p;
        });
        CHECK(seen.size() == 2);
    }
    SECTION("multi-X row: pivot is first X, other qubit untouched") {
        for (std::uint32_t s = 0; s < 16; ++s) {
            // State |++> written as H0 H1 CX01: Z1 row = X0X1.
            Tableau t(2);
            A::Comp(t)[0] = Pauli("XI");
            A::Comp(t)[1] = Pauli("XX");
            A::Had(t)[0] = Pauli("ZZ");
            A::Had(t)[1] = Pauli("IZ");
            A::Seed(s);
            bool p = A::RandomMeasurement(t, 1, Tableau::Basis::kComputational);
            CHECK(A::Comp(t)[1].x_.none());
            CHECK(A::Comp(t)[1].sign_ == p);
            CHECK(t.MeasureZ(1) == p);
            // Qubit 0 is still |+>.
            CHECK(t.MeasureX(0) == false);
            CheckSymplectic(t);
        }
    }
    SECTION("Y row on measured qubit") {
        for (std::uint32_t s = 0; s < 16; ++s) {
            Tableau t(1);
            t.ApplyH(0);
            t.ApplyS(0);  // had row = +Y
            A::Seed(s);
            bool p = A::RandomMeasurement(t, 0, Tableau::Basis::kHadamard);
            CHECK(A::Had(t)[0].x_.none());
            CHECK(A::Had(t)[0].sign_ == p);
            CheckSymplectic(t);
        }
    }
    SECTION("negative-sign row: returned outcome matches the collapsed row") {
        auto seen = Outcomes([](std::uint32_t s) {
            Tableau t(1);
            t.ApplyZ(0);
            t.ApplyH(0);  // comp row = -X
            A::Seed(s);
            bool p = A::RandomMeasurement(t, 0, Tableau::Basis::kComputational);
            CHECK(A::Comp(t)[0].x_.none());
            CHECK(A::Comp(t)[0].sign_ == p);
            return p;
        });
        CHECK(seen.size() == 2);
    }
    SECTION("same seed gives same outcome") {
        Tableau a(1), b(1);
        A::Seed(7);
        bool pa = A::RandomMeasurement(a, 0, Tableau::Basis::kHadamard);
        A::Seed(7);
        bool pb = A::RandomMeasurement(b, 0, Tableau::Basis::kHadamard);
        CHECK(pa == pb);
    }
}

}  // namespace stabaal
