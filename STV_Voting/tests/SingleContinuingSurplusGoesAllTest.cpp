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
            if (inQuotes && i + 1 < line.size() && line[i + 1] == '"') { cur.push_back('"'); ++i; }
            else { inQuotes = !inQuotes; }
        }
        else if (ch == ',' && !inQuotes) { cols.push_back(cur); cur.clear(); }
        else { cur.push_back(ch); }
    }
    cols.push_back(cur); return cols;
}

static std::vector<Row> parseRows(const std::string& csv) {
    std::vector<Row> rows;
    std::istringstream ss(csv); std::string line;
    while (std::getline(ss, line)) {
        if (line.empty() || line.rfind("Quota,", 0) == 0 || line.rfind("Round,", 0) == 0) continue;
        auto c = splitCsv(line); if (c.size() < 6) continue;
        Row r;
        try { r.round = std::stoi(c[0]); }
        catch (...) { continue; }
        r.cand = c[1];
        if (!r.cand.empty() && r.cand.front() == '"' && r.cand.back() == '"')
            r.cand = r.cand.substr(1, r.cand.size() - 2);
        if (!c[4].empty()) { try { r.transferred = std::stod(c[4]); } catch (...) { r.transferred = 0.0; } }
        r.sources = c[5];
        if (!r.sources.empty() && r.sources.front() == '"' && r.sources.back() == '"')
            r.sources = r.sources.substr(1, r.sources.size() - 2);
        rows.push_back(r);
    }
    return rows;
}

static bool approx(double a, double b) { return std::fabs(a - b) <= 0.01; }
static int fail(const char* m) { std::cout << "FAIL: " << m << "\n"; return 1; }

// Seats = 2; ballots = 6 => quota = ceil(600/3)/100 = 2.00
// A-only x5, B-only x1. A elected with 5.00, surplus = 3.00.
// First phase: ALL 5 A tickets are picked; each has no next.
// Continuing except A is {B} only, so per-ballot (3.00/5)=0.60 goes entirely to B on each ticket.
// Expect a round where B shows Transferred=3.00 with Sources including A(3.00); nobody else has A as source that round.
int main() {
    std::vector<std::vector<std::string>> ballots = {
        {"A"}, {"A"}, {"A"}, {"A"}, {"A"},
        {"B"}
    };

    // Capture CSV
    std::streambuf* old = std::cout.rdbuf();
    std::ostringstream cap; std::cout.rdbuf(cap.rdbuf());
    auto winners = runMultiSeatElection(ballots, 2);
    std::cout.rdbuf(old);

    if (winners.count("A") != 1 || winners.count("B") != 1 || winners.size() != 2)
        return fail("unexpected winners");

    const auto rows = parseRows(cap.str());

    // Find any round where B received exactly 3.00 from A
    int roundWithA = -1;
    double bTransferred = 0.0;
    for (const auto& r : rows) {
        if (r.cand == "B" && r.sources.find("A(3.00)") != std::string::npos) {
            roundWithA = r.round;
            bTransferred = r.transferred;
            break;
        }
    }
    if (roundWithA < 0) return fail("did not find A->B(3.00) transfer");
    if (!approx(bTransferred, 3.00)) return fail("B transferred != 3.00 in A-surplus round");

    // In the same round, no other candidate should list A as source (since only B is continuing)
    for (const auto& r : rows) {
        if (r.round != roundWithA) continue;
        if (r.cand != "B" && r.sources.find("A(") != std::string::npos)
            return fail("another candidate received transfer from A in A-surplus round");
    }

    std::cout << "SingleContinuingSurplusGoesAllTest: Passed.\n";
    return 0;
}