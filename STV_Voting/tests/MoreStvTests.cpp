#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <cmath>

double calculateQuota(int validBallots, int seats);
std::string runSingleSeatElection(const std::vector<std::vector<std::string>>& ballots);
std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>& ballots, int seats);

static int fail(const char* msg) { std::cout << "FAIL: " << msg << "\n"; return 1; }

int main()
{
    // Quota rounding (ceil to 2 decimals)
    {
        if (std::abs(calculateQuota(5, 2) - 1.67) > 1e-9) return fail("quota 5 ballots, 2 seats != 1.67");
        if (std::abs(calculateQuota(10, 3) - 2.50) > 1e-9) return fail("quota 10 ballots, 3 seats != 2.50");
        if (std::abs(calculateQuota(1, 1) - 0.50) > 1e-9) return fail("quota 1 ballot, 1 seat != 0.50");
    }

    // Single-seat: majority straight away
    {
        std::vector<std::vector<std::string>> ballots = {
            {"B"}, {"B"}, {"A"}
        };
        auto winner = runSingleSeatElection(ballots);
        if (winner != "B") return fail("single-seat majority winner mismatch (expected B)");
    }

    // Single-seat: winner via elimination
    {
        std::vector<std::vector<std::string>> ballots = {
            {"A"},
            {"B","A"},
            {"B"}
        };
        auto winner = runSingleSeatElection(ballots);
        if (winner != "B") return fail("single-seat elimination winner mismatch (expected B)");
    }

    // Multi-seat: remaining == seats -> elect all remaining
    {
        // 2 seats, 2 candidates, each with one ballot
        std::vector<std::vector<std::string>> ballots = {
            {"A"}, {"B"}
        };
        auto winners = runMultiSeatElection(ballots, 2);
        if (winners.size() != 2 || winners.count("A") != 1 || winners.count("B") != 1)
            return fail("force-elect remaining != expected {A,B}");
    }

    // Multi-seat: immediate quota winners (tickets engine path)
    {
        // 5 ballots, 2 seats -> quota = ceil(500/3)/100 = 1.67
        // A=2, C=2 meet quota, B=1 does not.
        std::vector<std::vector<std::string>> ballots = {
            {"A"}, {"A"},
            {"B"},
            {"C"}, {"C"}
        };
        auto winners = runMultiSeatElection(ballots, 2);
        if (winners.size() != 2 || winners.count("A") != 1 || winners.count("C") != 1)
            return fail("multi-seat immediate quota != expected {A,C}");
    }

    // Edge: no ballots
    {
        std::vector<std::vector<std::string>> ballots = {};
        auto winner = runSingleSeatElection(ballots);
        if (!winner.empty()) return fail("single-seat expected no winner for empty ballots");

        auto winners = runMultiSeatElection(ballots, 2);
        if (!winners.empty()) return fail("multi-seat expected no winners for empty ballots");
    }

    std::cout << "MoreStvTests: All tests passed.\n";
    return 0;
}