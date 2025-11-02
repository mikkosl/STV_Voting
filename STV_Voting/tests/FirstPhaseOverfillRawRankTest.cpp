#include <vector>
#include <string>
#include <set>
#include <iostream>

std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>&, int);

static int fail(const char* msg) { std::cout << "FAIL: " << msg << "\n"; return 1; }

int main() {
    // Seats = 2; 9 ballots => quota = ceil(900/3)/100 = 3.00
    // First choices: A=3, B=3, C=3 (overfill: 3 eligible > 2 seats)
    // Raw 2nd-choice counts:
    //  - A-second = 4 (from B,A,C x3 and C,A,B x1)
    //  - B-second = 5 (from A,B,C x3 and C,B,A x2)
    //  - C-second = 0
    // Expected winners by raw ranks: B then A (C loses).
    std::vector<std::vector<std::string>> ballots = {
        {"A","B","C"}, {"A","B","C"}, {"A","B","C"},
        {"B","A","C"}, {"B","A","C"}, {"B","A","C"},
        {"C","A","B"}, {"C","B","A"}, {"C","B","A"}
    };

    auto winners = runMultiSeatElection(ballots, 2);

    if (winners.size() != 2) return fail("expected exactly 2 winners");
    if (!winners.count("B")) return fail("expected B to win by highest 2nd-choice");
    if (!winners.count("A")) return fail("expected A to win by next-highest 2nd-choice");
    if (winners.count("C"))  return fail("did not expect C to win");

    std::cout << "FirstPhaseOverfillRawRankTest: Passed.\n";
    return 0;
}