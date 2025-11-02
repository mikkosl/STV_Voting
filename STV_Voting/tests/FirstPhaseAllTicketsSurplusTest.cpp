#include <vector>
#include <string>
#include <set>
#include <sstream>
#include <iostream>
#include <cmath>

std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>&, int);

// Minimal CSV splitter
static std::vector<std::string> splitCsv(const std::string& line) {
    std::vector<std::string> cols; std::string cur; bool inQuotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char ch = line[i];
        if (ch == '"') {
            if (inQuotes && i + 1 < line.size() && line[i+1] == '"') { cur.push_back('"'); ++i; }
            else { inQuotes = !inQuotes; }
        } else if (ch == ',' && !inQuotes) {
            cols.push_back(cur); cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    cols.push_back(cur);
    return cols;
}

static bool tryGetRow(const std::string& csv, int round, const std::string& candidate,
                      double& transferredOut, std::string& sources)
{
    std::istringstream ss(csv);
    std::string line;
    while (std::getline(ss, line)) {
        if (line.empty()) continue;
        if (line.rfind("Quota,", 0) == 0) continue;
        if (line.rfind("Round,", 0) == 0) continue;

        auto cols = splitCsv(line);
        if (cols.size() < 6) continue;

        int r = 0;
        try { r = std::stoi(cols[0]); } catch (...) { continue; }

        std::string cand = cols[1];
        if (!cand.empty() && cand.front()=='"' && cand.back()=='"') cand = cand.substr(1, cand.size()-2);

        if (r != round || cand != candidate) continue;

        transferredOut = 0.0;
        if (!cols[4].empty()) {
            try { transferredOut = std::stod(cols[4]); } catch (...) { transferredOut = 0.0; }
        }

        sources = cols[5];
        if (!sources.empty() && sources.front()=='"' && sources.back()=='"')
            sources = sources.substr(1, sources.size()-2);

        return true;
    }
    return false;
}

static bool approx(double a, double b) { return std::fabs(a - b) <= 0.01; }
static int fail(const char* m){ std::cout << "FAIL: " << m << "\n"; return 1; }

// Scenario (seats=2, first phase):
// Ballots: 2x A,B,C; 1x A; 1x B; 1x C
// Quota = ceil(5/3 to 2dp) = 1.67
// A has 3.00 -> surplus 1.33; first-phase uses ALL A tickets (2 with next=B, 1 with no-next).
// perBallot = floor2(1.33/3) = 0.44
// Transfers: B +0.88 (two A,B,C) +0.22 (split of A-only across {B,C}) = 1.10; C +0.22
int main() {
    std::vector<std::vector<std::string>> ballots = {
        {"A","B","C"}, {"A","B","C"}, {"A"}, {"B"}, {"C"}
    };

    // Capture CSV
    std::streambuf* oldBuf = std::cout.rdbuf();
    std::ostringstream cap; std::cout.rdbuf(cap.rdbuf());

    auto winners = runMultiSeatElection(ballots, 2);

    std::cout.rdbuf(oldBuf);
    const std::string out = cap.str();

    // Winners should be A and B (not strictly required for this assertion, but a sanity check)
    if (winners.size() != 2 || winners.count("A") != 1 || winners.count("B") != 1)
        return fail("unexpected winners (sanity)");

    // Round 1 must reflect transfers computed using ALL A tickets
    double tB=0, tC=0; std::string sB, sC;
    if (!tryGetRow(out, 1, "B", tB, sB)) return fail("missing round 1 row for B");
    if (!tryGetRow(out, 1, "C", tC, sC)) return fail("missing round 1 row for C");

    if (!approx(tB, 1.10)) return fail("Round 1 B transferred != 1.10 (expected all A tickets used)");
    if (!approx(tC, 0.22)) return fail("Round 1 C transferred != 0.22 (expected split from A-only)");

    if (sB.find("A(1.10)") == std::string::npos) return fail("Round 1 B sources missing A(1.10)");
    if (sC.find("A(0.22)") == std::string::npos) return fail("Round 1 C sources missing A(0.22)");

    std::cout << "FirstPhaseAllTicketsSurplusTest: Passed.\n";
    return 0;
}