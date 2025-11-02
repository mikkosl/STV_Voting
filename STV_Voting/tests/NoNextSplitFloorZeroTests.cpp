#include <vector>
#include <string>
#include <set>
#include <sstream>
#include <iostream>
#include <cmath>

// We will make perBallot = 0.01 and continuing recipients >= 3 so split=floor2(0.01/3)=0.00
std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>&, int);

static int fail(const char* m){ std::cout << "FAIL: " << m << "\n"; return 1; }
static std::vector<std::string> lines(const std::string& s) { std::vector<std::string> v; std::istringstream ss(s); std::string l; while (std::getline(ss,l)) v.push_back(l); return v; }

int main() {
    // Seats = 2; n = 100 -> quota = ceil(10000/3)/100 = 33.34
    // Winner A has 34 one-choice ballots: A (no next).
    // Others: B=33, C=33. A’s surplus = 0.66. denom = 34 (all A tickets have no next).
    // perBallot = floor2(0.66/34) = 0.01. continuing recipients = {B,C} (2) → split=floor2(0.01/2)=0.00.
    // Expect: No transfers to B or C from A’s surplus (split floors to 0), totals unchanged aside from A decreasing by 0.34 (34×0.01).
    std::vector<std::vector<std::string>> ballots;
    ballots.reserve(100);
    for (int i=0;i<34;++i) ballots.push_back({"A"});
    for (int i=0;i<33;++i) ballots.push_back({"B"});
    for (int i=0;i<33;++i) ballots.push_back({"C"});

    std::streambuf* old = std::cout.rdbuf();
    std::ostringstream cap; std::cout.rdbuf(cap.rdbuf());
    auto winners = runMultiSeatElection(ballots, 2);
    std::cout.rdbuf(old);
    const std::string out = cap.str();

    // Winners must be A and then B or C. But importantly, in the round after A’s surplus,
    // there should be no positive Transferred for B or C from A.
    bool sawAnyBCTransferFromA = false;
    for (const auto& ln : lines(out)) {
        if (ln.find(",\"B\",") != std::string::npos || ln.find(",\"C\",") != std::string::npos) {
            // Check Sources column contains A(…)
            if (ln.find(",\"") != std::string::npos && ln.find("A(") != std::string::npos) {
                // If there is an A(…) entry, the amount should be 0.00 due to split=0.00
                // Our CSV prints Transferred and Sources; we assert no A(positive)
                if (ln.find("A(0.00)") == std::string::npos) {
                    sawAnyBCTransferFromA = true;
                    break;
                }
            }
        }
    }
    if (sawAnyBCTransferFromA) return fail("unexpected positive split from A to B/C");

    std::cout << "NoNextSplitFloorZeroTests: Passed.\n";
    return 0;
}