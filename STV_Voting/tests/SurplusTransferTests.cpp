#include <vector>
#include <string>
#include <set>
#include <sstream>
#include <iostream>
#include <cmath>

std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>&, int);

// Minimal CSV splitter for our output
static std::vector<std::string> splitCsv(const std::string& line) {
    std::vector<std::string> cols; std::string cur; bool inQuotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char ch = line[i];
        if (ch == '"') { if (inQuotes && i + 1 < line.size() && line[i+1] == '"') { cur.push_back('"'); ++i; } else { inQuotes = !inQuotes; } }
        else if (ch == ',' && !inQuotes) { cols.push_back(cur); cur.clear(); }
        else cur.push_back(ch);
    }
    cols.push_back(cur); return cols;
}

static bool tryGetRow(const std::string& csv, int round, const std::string& candidate,
                      double& transferred, std::string& sources) {
    std::istringstream ss(csv); std::string line;
    while (std::getline(ss, line)) {
        if (line.empty() || line.rfind("Quota,",0)==0 || line.rfind("Round,",0)==0) continue;
        auto c = splitCsv(line); if (c.size() < 6) continue;
        int r = 0; try { r = std::stoi(c[0]); } catch (...) { continue; }
        auto cand = c[1]; if (!cand.empty() && cand.front()=='"' && cand.back()=='"') cand = cand.substr(1, cand.size()-2);
        if (r != round || cand != candidate) continue;
        transferred = 0.0; if (!c[4].empty()) try { transferred = std::stod(c[4]); } catch (...) {}
        sources = c[5]; if (!sources.empty() && sources.front()=='"' && sources.back()=='"') sources = sources.substr(1, sources.size()-2);
        return true;
    }
    return false;
}

static bool approx(double a, double b) { return std::fabs(a - b) <= 0.01; }
static int fail(const char* m){ std::cout << "FAIL: " << m << "\n"; return 1; }

// Seats=2
// Ballots: 2x A,B,C; 1x A; 1x B; 1x C
// Quota = ceil(5/3 to 2dp) = 1.67
// First phase: A elected with 3.00; surplus=1.33; perBallot=floor2(1.33/3)=0.44 using ALL A tickets.
//   - Two A,B,C -> B +0.88; one A-only has no-next -> split among remaining (B,C), split=floor2(0.44/2)=0.22 => B +0.22, C +0.22.
//   Expect Round 1: B Transferred=1.10 from A; C Transferred=0.22 from A.
// Later phase: B now ≥ quota; B surplus from last batch only (tickets from A’s transfer): 0.44,0.44,0.22
//   denom=3; perBallot=floor2(0.43/3)=0.14; two with next=C => +0.28; 0.22 no-next -> split to C (only remaining) => +0.14
//   Expect Round 2: C Transferred=0.42 from B.
int main() {
    std::vector<std::vector<std::string>> ballots = {
        {"A","B","C"}, {"A","B","C"}, {"A"}, {"B"}, {"C"}
    };

    // Capture program CSV
    std::streambuf* old = std::cout.rdbuf();
    std::ostringstream cap; std::cout.rdbuf(cap.rdbuf());

    auto winners = runMultiSeatElection(ballots, 2);

    std::cout.rdbuf(old);
    const std::string out = cap.str();

    if (winners.count("A") != 1 || winners.count("B") != 1 || winners.size() != 2)
        return fail("winners mismatch");

    // Round 1: proportional surplus (ALL A tickets)
    double tB=0, tC=0; std::string sB, sC;
    if (!tryGetRow(out, 1, "B", tB, sB)) return fail("missing R1 B");
    if (!tryGetRow(out, 1, "C", tC, sC)) return fail("missing R1 C");
    if (!approx(tB, 1.10)) return fail("R1 B transferred != 1.10");
    if (!approx(tC, 0.22)) return fail("R1 C transferred != 0.22");
    if (sB.find("A(1.10)") == std::string::npos) return fail("R1 B missing A(1.10)");
    if (sC.find("A(0.22)") == std::string::npos) return fail("R1 C missing A(0.22)");

    // Round 2: last-batch-only surplus from B
    double tC2=0; std::string sC2;
    if (!tryGetRow(out, 2, "C", tC2, sC2)) return fail("missing R2 C");
    if (!approx(tC2, 0.42)) return fail("R2 C transferred != 0.42");
    if (sC2.find("B(0.42)") == std::string::npos) return fail("R2 C missing B(0.42)");

    std::cout << "SurplusTransferTests: All tests passed.\n";
    return 0;
}