#include <vector>
#include <string>
#include <set>
#include <iostream>

std::string runSingleSeatElection(const std::vector<std::vector<std::string>>& ballots);
std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>& ballots, int seats);

static int fail(const char* msg) { std::cout << "FAIL: " << msg << "\n"; return 1; }

int main()
{
    // Single-seat: A wins 2-1
    {
        std::vector<std::vector<std::string>> ballots = {
            {"A","B"},
            {"A"},
            {"B","A"}
        };
        auto winner = runSingleSeatElection(ballots);
        if (winner != "A") {
            std::cout << "Winner=" << winner << " expected=A\n";
            return fail("single-seat winner mismatch");
        }
    }

    // Multi-seat (tickets on): 2 seats; A=3, B=2, C=1 -> expect A,B
    {
        std::vector<std::vector<std::string>> ballots = {
            {"A"}, {"A"}, {"A"},
            {"B"}, {"B"},
            {"C"}
        };
        auto winners = runMultiSeatElection(ballots, 2);
        if (winners.count("A") != 1 || winners.count("B") != 1 || winners.size() != 2) {
            std::cout << "Winners:";
            for (auto& w : winners) std::cout << " " << w;
            std::cout << " expected: A B\n";
            return fail("multi-seat winners mismatch");
        }
    }

    std::cout << "All tests passed.\n";
    return 0;
}