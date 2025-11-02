#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <sstream>
#include <iomanip>

std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>& ballots, int seats);

// Minimal CSV splitter (handles quoted candidate and quoted Sources at tail)
static std::vector<std::string> splitCsv(const std::string& line)
{
    std::vector<std::string> cols;
    std::string cur;
    bool inQuotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char ch = line[i];
        if (ch == '"' ) {
            if (inQuotes && i + 1 < line.size() && line[i+1] == '"') {
                cur.push_back('"'); ++i; // escaped quote
            } else {
                inQuotes = !inQuotes;
            }
        } else if (ch == ',' && !inQuotes) {
            cols.push_back(cur);
            cur.clear();
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

        // cols: 0=Round, 1=Candidate (quoted), 2=Votes, 3=Status, 4=Transferred, 5=Sources (quoted or empty)
        int r = 0;
        try { r = std::stoi(cols[0]); } catch (...) { continue; }

        // candidate is quoted, remove quotes if present
        std::string cand = cols[1];
        if (!cand.empty() && cand.front() == '"' && cand.back() == '"') {
            cand = cand.substr(1, cand.size()-2);
        }

        if (r != round || cand != candidate) continue;

        // Transferred can be empty
        transferredOut = 0.0;
        if (!cols[4].empty()) {
            try { transferredOut = std::stod(cols[4]); } catch (...) { transferredOut = 0.0; }
        }

        // Sources: remove quotes if present
        sources = cols[5];
        if (!sources.empty() && sources.front() == '"' && sources.back() == '"') {
            sources = sources.substr(1, sources.size()-2);
        }
        return true;
    }
    return false;
}

static int fail(const char* msg) { std::cout << "FAIL: " << msg << "\n"; return 1; }

// Scenario:
// Seats=2, ballots:
//  - A,B,C
//  - A,B,C
//  - A
//  - B
//  - C
// Quota = ceil(500/3)/100 = 1.67
// Round 1: A elected (3.00), surplus 1.33 distributed proportionally from ALL A tickets:
//   perBallot = floor2(1.33/3)=0.44
//   B gets 0.88 (two A,B,C) + 0.22 split from A-only => 1.10; C gets 0.22
//   Expect: Round 1: B Transferred=1.10 from A; C Transferred=0.22 from A
// Next phase: B now meets quota (2.10), elected; its surplus arises later => last-batch-only.
//   B surplus = 0.43. Last-batch at B are 0.44,0.44,0.22 fragments from A’s transfer.
//   denom=3, perBallot=0.14.
//   Two 0.44 have next=C => C +0.28; the 0.22 has no next => split among continuing (only C) => +0.14.
//   Expect: Round 2: C Transferred=0.42 from B.
int main()
{
    std::vector<std::vector<std::string>> ballots = {
        {"A","B","C"},
        {"A","B","C"},
        {"A"},
        {"B"},
        {"C"}
    };

    // Capture stdout
    std::streambuf* oldBuf = std::cout.rdbuf();
    std::ostringstream capture;
    std::cout.rdbuf(capture.rdbuf());

    auto winners = runMultiSeatElection(ballots, 2);

    // Restore stdout
    std::cout.rdbuf(oldBuf);
    const std::string out = capture.str();

    // Sanity winners
    if (winners.count("A") != 1 || winners.count("B") != 1 || winners.size() != 2) {
        std::cout << "Winners:";
        for (auto& w : winners) std::cout << " " << w;
        std::cout << " expected: A B\n";
        return fail("winners mismatch");
    }

    // Round 1 assertions (first phase proportional)
    {
        double tB = 0.0, tC = 0.0; std::string srcB, srcC;
        if (!tryGetRow(out, 1, "B", tB, srcB)) return fail("missing round 1 row for B");
        if (!tryGetRow(out, 1, "C", tC, srcC)) return fail("missing round 1 row for C");

        auto approx = [](double a, double b){ return std::fabs(a - b) <= 0.01; };
        if (!approx(tB, 1.10)) { std::cout << "Round1 B transferred=" << tB << " expected=1.10\n"; return fail("round1 B transferred mismatch"); }
        if (!approx(tC, 0.22)) { std::cout << "Round1 C transferred=" << tC << " expected=0.22\n"; return fail("round1 C transferred mismatch"); }

        if (srcB.find("A(1.10)") == std::string::npos)
            return fail("round1 B sources missing A(1.10)");
        if (srcC.find("A(0.22)") == std::string::npos)
            return fail("round1 C sources missing A(0.22)");
    }

    // Round 2 assertions (later phase last-batch only)
    {
        double tC = 0.0; std::string srcC;
        if (!tryGetRow(out, 2, "C", tC, srcC)) return fail("missing round 2 row for C");

        auto approx = [](double a, double b){ return std::fabs(a - b) <= 0.01; };
        if (!approx(tC, 0.42)) { std::cout << "Round2 C transferred=" << tC << " expected=0.42\n"; return fail("round2 C transferred mismatch"); }

        if (srcC.find("B(0.42)") == std::string::npos)
            return fail("round2 C sources missing B(0.42)");
    }

    std::cout << "SurplusTransferTests: All tests passed.\n";
    return 0;
}