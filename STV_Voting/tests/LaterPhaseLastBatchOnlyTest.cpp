#include <vector>
#include <string>
#include <set>
#include <sstream>
#include <iostream>
#include <cmath>

std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>&, int);

struct Row {
    int round = -1;
    std::string cand;
    double transferred = 0.0;
    std::string sources;
};

static std::vector<std::string> splitCsv(const std::string& line) {
    std::vector<std::string> cols; std::string cur; bool inQuotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char ch = line[i];
        if (ch == '"') {
            if (inQuotes && i + 1 < line.size() && line[i+1] == '"') { cur.push_back('"'); ++i; }
            else { inQuotes = !inQuotes; }
        } else if (ch == ',' && !inQuotes) { cols.push_back(cur); cur.clear(); }
        else { cur.push_back(ch); }
    }
    cols.push_back(cur); return cols;
}

static std::vector<Row> parseRows(const std::string& csv) {
    std::vector<Row> rows;
    std::istringstream ss(csv); std::string line;
    while (std::getline(ss, line)) {
        if (line.empty() || line.rfind("Quota,",0)==0 || line.rfind("Round,",0)==0) continue;
        auto c = splitCsv(line); if (c.size() < 6) continue;
        Row r;
        try { r.round = std::stoi(c[0]); } catch (...) { continue; }
        r.cand = c[1];
        if (!r.cand.empty() && r.cand.front()=='"' && r.cand.back()=='"')
            r.cand = r.cand.substr(1, r.cand.size()-2);
        if (!c[4].empty()) { try { r.transferred = std::stod(c[4]); } catch (...) { r.transferred = 0.0; } }
        r.sources = c[5];
        if (!r.sources.empty() && r.sources.front()=='"' && r.sources.back()=='"')
            r.sources = r.sources.substr(1, r.sources.size()-2);
        rows.push_back(r);
    }
    return rows;
}

static bool approx(double a, double b) { return std::fabs(a - b) <= 0.01; }
static int fail(const char* m){ std::cout << "FAIL: " << m << "\n"; return 1; }

int main() {
    // Seats = 2; 9 ballots -> quota = 3.00
    // Later-phase: B crosses quota after D elimination; last-batch-only means only D->B (1.00) is used.
    // B surplus = 0.41; no-next after B => all to the only continuing C => C +0.41 from B.
    std::vector<std::vector<std::string>> ballots = {
        {"A","B"},  // A with next B
        {"A","C"},  // A with next C
        {"A"}, {"A"}, // A-only x2
        {"B"}, {"B"},
        {"C"}, {"C"},
        {"D","B"}   // D with next B
    };

    // Capture CSV
    std::streambuf* old = std::cout.rdbuf();
    std::ostringstream cap; std::cout.rdbuf(cap.rdbuf());
    auto winners = runMultiSeatElection(ballots, 2);
    std::cout.rdbuf(old);
    const std::string out = cap.str();

    if (winners.count("A") != 1 || winners.count("B") != 1 || winners.size() != 2)
        return fail("winners mismatch");

    // Parse all rows and find the round where C receives exactly 0.41 from B.
    const auto rows = parseRows(out);
    int roundWithB = -1;
    double cTransferred = 0.0;
    for (const auto& r : rows) {
        if (r.cand == "C" && r.sources.find("B(") != std::string::npos) {
            // Optional: ensure the amount from B is 0.41 in the Sources string
            if (r.sources.find("B(0.41)") != std::string::npos) {
                roundWithB = r.round;
                cTransferred = r.transferred;
                break;
            }
        }
    }
    if (roundWithB < 0) return fail("did not find B->C(0.41) transfer in any round");
    if (!approx(cTransferred, 0.41)) return fail("C transferred != 0.41 in the B-surplus round");

    // In that same round, no other candidate should list B as a source (last-batch-only, no-next -> only C gets it).
    for (const auto& r : rows) {
        if (r.round != roundWithB) continue;
        if (r.cand != "C" && r.sources.find("B(") != std::string::npos) {
            return fail("another candidate received transfer from B in B-surplus round");
        }
    }

    std::cout << "LaterPhaseLastBatchOnlyTest: Passed.\n";
    return 0;
}