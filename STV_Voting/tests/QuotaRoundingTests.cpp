#include <vector>
#include <string>
#include <set>
#include <sstream>
#include <iostream>

std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>&, int);

static int fail(const char* m){ std::cout << "FAIL: " << m << "\n"; return 1; }
static bool hasQuota(const std::string& out, const std::string& q) {
    return out.find("Quota," + q) != std::string::npos;
}

int main() {
    {
        // n=2, seats=2 -> quota = ceil(200/3)/100 = 0.67
        std::vector<std::vector<std::string>> ballots = { {"A"}, {"B"} };
        std::streambuf* old = std::cout.rdbuf();
        std::ostringstream cap; std::cout.rdbuf(cap.rdbuf());
        (void)runMultiSeatElection(ballots, 2);
        std::cout.rdbuf(old);
        if (!hasQuota(cap.str(), "0.67")) return fail("quota 0.67 mismatch");
    }
    {
        // n=5, seats=2 -> quota = ceil(500/3)/100 = 1.67
        std::vector<std::vector<std::string>> ballots = { {"A"},{"B"},{"C"},{"D"},{"E"} };
        std::streambuf* old = std::cout.rdbuf();
        std::ostringstream cap; std::cout.rdbuf(cap.rdbuf());
        (void)runMultiSeatElection(ballots, 2);
        std::cout.rdbuf(old);
        if (!hasQuota(cap.str(), "1.67")) return fail("quota 1.67 mismatch");
    }
    {
        // n=101, seats=2 -> quota = ceil(10100/3)/100 = ceil(3366.66)/100 = 33.67
        std::vector<std::vector<std::string>> ballots(101, std::vector<std::string>{"A"});
        std::streambuf* old = std::cout.rdbuf();
        std::ostringstream cap; std::cout.rdbuf(cap.rdbuf());
        (void)runMultiSeatElection(ballots, 2);
        std::cout.rdbuf(old);
        if (!hasQuota(cap.str(), "33.67")) return fail("quota 33.67 mismatch");
    }
    std::cout << "QuotaRoundingTests: All tests passed.\n";
    return 0;
}