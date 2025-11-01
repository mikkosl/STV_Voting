#include <vector>
#include <string>
#include <set>
#include <iostream>

std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>& ballots, int seats);

static int fail(const char* msg) { std::cout << "FAIL: " << msg << "\n"; return 1; }

int main()
{
    // Test 1: A elected; A's last batch has no next -> split equally over continuing {B,C}.
    // 4 ballots, 2 seats. Quota = ceil(400/3)/100 = 1.34
    // A=2 ({"A","B"}, {"A"} [last with no next]), B=1 ({"B","C"}), C only appears as lower pref.
    // A surplus = 0.66; last batch denom=1; perBallot=0.66; split across {B,C} => +0.33 each (floored).
    // Expected winners: {A,B} deterministically (C stays at 0; B ends up ahead; C will be eliminated).
    {
        std::vector<std::vector<std::string>> ballots = {
            {"A","B"}, // A first, B ensures B is continuing and in allCandidates
            {"A"},     // A-only, last ticket at A has no next
            {"B","C"}  // B gets 1; C enters allCandidates but not first-choice
        };
        auto winners = runMultiSeatElection(ballots, 2);
        if (winners.count("A") != 1 || winners.count("B") != 1 || winners.size() != 2) {
            std::cout << "Winners:";
            for (auto& w : winners) std::cout << " " << w;
            std::cout << " expected: A B\n";
            return fail("last-batch no-next split (multi-continuing) mismatch");
        }
    }

    // Test 2: A and C both reach quota initially; when transferring A's surplus,
    // cont excludes elected {A,C}, so only B is continuing; A's last-batch no-next
    // surplus flows entirely to B (cont size = 1). Winners remain {A,C}.
    // 3 ballots, 2 seats. Quota = ceil(300/3)/100 = 1.00
    // A=2 ({"A","B"}, {"A"} [last with no next]), C=1 ({"C"}), B only present as lower pref.
    {
        std::vector<std::vector<std::string>> ballots = {
            {"A","B"}, // A first, B as lower pref (enters allCandidates)
            {"A"},     // A-only, last ticket at A has no next
            {"C"}      // C is also eligible in same phase
        };
        auto winners = runMultiSeatElection(ballots, 2);
        if (winners.count("A") != 1 || winners.count("C") != 1 || winners.size() != 2) {
            std::cout << "Winners:";
            for (auto& w : winners) std::cout << " " << w;
            std::cout << " expected: A C\n";
            return fail("last-batch no-next split (single continuing) mismatch");
        }
    }

    std::cout << "LastBatchNoNextTests: All tests passed.\n";
    return 0;
}