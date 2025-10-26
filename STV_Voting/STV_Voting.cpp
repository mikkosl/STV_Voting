// STV_Voting.cpp : Defines the entry point for the application.
//

#include "STV_Voting.h"

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <limits> // <-- add this
#include <random> // add this for drawing lots
#ifdef _WIN32
#include <conio.h>
#endif

using namespace std;

// Numeric epsilon for FP comparisons in tie logic
constexpr double kEps = 1e-9;

// Multi-seat Single Transferable Vote

// Step 1: Calculate the election quota
double calculateQuota(int validBallots, int seats);

// Step 2: Count first-choice votes
std::map<std::string, double> countFirstChoices(const std::vector<std::vector<std::string>>& ballots);

// Step 3: Elect candidates who meet the quota

// Step 4: Transfer surplus
void transferSurplus(std::string candidate, double surplus, const std::vector<std::vector<std::string>>& ballots, std::map<std::string, double>& voteCounts, const std::set<std::string>& elected);

// Step 5: Eliminate the lowest-vote candidate

// Step 6: Transfer votes from the eliminated candidate
void transferEliminatedVotes(std::string eliminated, const std::vector<std::vector<std::string>>& ballots, std::map<std::string, double>& voteCounts, const std::set<std::string>& elected);

// Step 7: Resolve ties per §11D
std::string resolveTieSTV(const std::vector<std::vector<std::string>>& ballots, const std::set<std::string>& tiedCandidates);

// Tie-breaking with round history (updated signature to accept current state)
std::string resolveTieSTVWithHistory(const std::vector<std::map<std::string, double>>& history,
                                     const std::vector<std::vector<std::string>>& ballots,
                                     const std::set<std::string>& tiedCandidates,
                                     const std::map<std::string, double>& voteCounts,
                                     const std::set<std::string>& elected);

// --- Implementations ---

// Add this helper function near the top of the file, before its first use:
static double floor2(double x) { return std::floor(x * 100.0) / 100.0; }

// Add this helper function near the other tie-break helpers (above its first use)
static std::string finalTiebreakPick(const std::set<std::string>& tiedCandidates)
{
    if (tiedCandidates.empty()) return {};

    // Draw lots uniformly among the tied candidates
    std::vector<std::string> pool(tiedCandidates.begin(), tiedCandidates.end());
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, pool.size() - 1);
    return pool[dist(gen)];
}

// Add this helper function above its first use (e.g., above resolveTieSTVWithHistory)
static std::string resolveTieNextContinuing(const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& tiedCandidates,
    const std::map<std::string, double>& voteCounts,
    const std::set<std::string>& elected)
{
    if (tiedCandidates.empty()) return {};

    auto isContinuing = [&](const std::string& c) -> bool {
        return elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end();
        };

    size_t maxLen = 0;
    std::vector<std::vector<std::string>> compressed;
    compressed.reserve(ballots.size());
    for (const auto& b : ballots) {
        std::vector<std::string> v;
        v.reserve(b.size());
        for (const auto& c : b) if (isContinuing(c)) v.push_back(c);
        maxLen = std::max(maxLen, v.size());
        compressed.push_back(std::move(v));
    }
    if (maxLen <= 1) return {};

    // For rank = 2..maxLen, pick the candidate with the highest count at that rank.
    for (size_t rank = 2; rank <= maxLen; ++rank) {
        std::map<std::string, int> rankCounts;
        for (const auto& c : tiedCandidates) rankCounts[c] = 0;

        for (const auto& v : compressed) {
            if (v.size() < rank) continue;
            const std::string& candAtRank = v[rank - 1];
            if (tiedCandidates.count(candAtRank)) rankCounts[candAtRank] += 1;
        }

        int best = -1;
        std::set<std::string> bestCands;
        for (const auto& [cand, cnt] : rankCounts) {
            if (cnt > best) {
                best = cnt;
                bestCands.clear();
                bestCands.insert(cand);
            }
            else if (cnt == best) {
                bestCands.insert(cand);
            }
        }

        if (bestCands.size() == 1) return *bestCands.begin();
        // else continue deeper
    }
    return {};
}

// Add this function above its first use (e.g., above resolveTieSTVWithHistoryForElimination)
// Original prefs, pick lowest (for elimination tiebreak)
static std::string resolveTieSTVLeast(const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& tiedCandidates)
{
    if (tiedCandidates.empty()) return {};
    size_t maxLen = 0;
    for (const auto& b : ballots) maxLen = std::max(maxLen, b.size());
    if (maxLen == 0) return finalTiebreakPick(tiedCandidates);

    for (size_t rank = 2; rank <= maxLen; ++rank) {
        std::map<std::string, int> rankCounts;
        for (const auto& c : tiedCandidates) rankCounts[c] = 0;

        for (const auto& b : ballots) {
            if (b.size() < rank) continue;
            const std::string& candAtRank = b[rank - 1];
            if (tiedCandidates.count(candAtRank)) rankCounts[candAtRank] += 1;
        }

        int best = std::numeric_limits<int>::max();
        std::set<std::string> bestCands;
        for (const auto& [cand, cnt] : rankCounts) {
            if (cnt < best) { best = cnt; bestCands.clear(); bestCands.insert(cand); }
            else if (cnt == best) { bestCands.insert(cand); }
        }
        if (bestCands.size() == 1) return *bestCands.begin();
    }
    return finalTiebreakPick(tiedCandidates);
}

// Helper: next-continuing-preference tiebreak (least preferred for elimination)
static std::string resolveTieNextContinuingLeast(const std::vector<std::vector<std::string>>& ballots,
                                                 const std::set<std::string>& tiedCandidates,
                                                 const std::map<std::string, double>& voteCounts,
                                                 const std::set<std::string>& elected)
{
    if (tiedCandidates.empty()) return {};

    auto isContinuing = [&](const std::string& c) -> bool {
        return elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end();
    };

    size_t maxLen = 0;
    std::vector<std::vector<std::string>> compressed;
    compressed.reserve(ballots.size());
    for (const auto& b : ballots) {
        std::vector<std::string> v;
        v.reserve(b.size());
        for (const auto& c : b) if (isContinuing(c)) v.push_back(c);
        maxLen = std::max(maxLen, v.size());
        compressed.push_back(std::move(v));
    }
    if (maxLen <= 1) return {};

    // For rank = 2..maxLen, pick the candidate with the lowest count at that rank.
    for (size_t rank = 2; rank <= maxLen; ++rank) {
        std::map<std::string, int> rankCounts;
        for (const auto& c : tiedCandidates) rankCounts[c] = 0;

        for (const auto& v : compressed) {
            if (v.size() < rank) continue;
            const std::string& candAtRank = v[rank - 1];
            if (tiedCandidates.count(candAtRank)) rankCounts[candAtRank] += 1;
        }

        int best = std::numeric_limits<int>::max();
        std::set<std::string> bestCands;
        for (const auto& [cand, cnt] : rankCounts) {
            if (cnt < best) {
                best = cnt;
                bestCands.clear();
                bestCands.insert(cand);
            } else if (cnt == best) {
                bestCands.insert(cand);
            }
        }

        if (bestCands.size() == 1) return *bestCands.begin();
        // else continue deeper
    }
    return {};
}

// §11D (elimination context): prefer lowest earlier total, then least next-continuing preferences
static std::string resolveTieSTVWithHistoryForElimination(const std::vector<std::map<std::string, double>>& history,
                                                          const std::vector<std::vector<std::string>>& ballots,
                                                          const std::set<std::string>& tiedCandidates,
                                                          const std::map<std::string, double>& voteCounts,
                                                          const std::set<std::string>& elected)
{
    if (tiedCandidates.empty()) return {};

    // Primary: use the most recent earlier phase where totals differ, preferring the lowest
    for (auto it = history.rbegin(); it != history.rend(); ++it) {
        const auto& totals = *it;

        double minV = std::numeric_limits<double>::infinity();
        for (const auto& cand : tiedCandidates) {
            auto f = totals.find(cand);
            const double v = (f == totals.end() ? 0.0 : f->second);
            if (v < minV) minV = v;
        }

        std::set<std::string> bestCands;
        for (const auto& cand : tiedCandidates) {
            auto f = totals.find(cand);
            const double v = (f == totals.end() ? 0.0 : f->second);
            if (std::abs(v - minV) <= kEps) bestCands.insert(cand);
        }
        if (bestCands.size() == 1) return *bestCands.begin();
    }

    // Secondary: least next continuing preferences (use continuing compression)
    if (auto chosen = resolveTieNextContinuingLeast(ballots, tiedCandidates, voteCounts, elected); !chosen.empty()) {
        return chosen;
    }

    return finalTiebreakPick(tiedCandidates);
}


// Print CSV header for election rounds
void printCsvHeader()
{
    std::cout << "Round";
    std::cout << ",Candidate";
    std::cout << ",Votes";
    std::cout << ",Status";
    std::cout << ",Transferred";
    std::cout << ",Sources";
    std::cout << "\n";
}

// Print a CSV row for each candidate in the round, including transfers and source breakdowns.
void printCsvRound(
    int round,
    const std::map<std::string, double>& voteCounts,
    const std::set<std::string>& elected,
    const std::set<std::string>& allCandidates,
    const std::map<std::string, double>& transferredAmounts,
    const std::map<std::string, std::vector<std::pair<std::string, double>>>& transferSources)
{
    for (const auto& cand : allCandidates) {
        std::cout << round;
        std::cout << "," << cand;

        auto it = voteCounts.find(cand);
        double votes = (it != voteCounts.end()) ? it->second : 0.0;
        std::cout << "," << std::fixed << std::setprecision(2) << votes;

        const bool isElected = elected.count(cand) != 0;
        std::string status;
        if (isElected) {
            status = "Elected";
        } else if (voteCounts.find(cand) == voteCounts.end()) {
            status = "Eliminated";
        } else {
            status = "Continuing";
        }
        std::cout << "," << status;

        // Transferred amount (total received this round)
        std::cout << ",";
        auto tIt = transferredAmounts.find(cand);
        if (tIt != transferredAmounts.end()) {
            std::cout << std::fixed << std::setprecision(2) << tIt->second;
        }

        // All sources with per-source amounts: e.g., Alice(1.20); Bob(0.80)
        std::cout << ",";
        auto sIt = transferSources.find(cand);
        if (sIt != transferSources.end()) {
            const auto& items = sIt->second;
            for (size_t i = 0; i < items.size(); ++i) {
                const auto& src = items[i];
                std::cout << src.first << "(" << std::fixed << std::setprecision(2) << src.second << ")";
                if (i + 1 < items.size()) std::cout << "; ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

// Add this function above its first use (e.g., above runMultiSeatElection)
// Computes the distribution of surplus votes for logging purposes.
// Returns a map of recipient candidate to amount received.
std::map<std::string, double> computeSurplusDistributionForLog(
    const std::string& candidate,
    double surplus,
    const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& elected,
    const std::map<std::string, double>& voteCounts)
{
    std::map<std::string, double> result;
    if (surplus <= 0.0) return result;

    auto isContinuing = [&](const std::string& c) -> bool {
        return elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end();
    };

    int transferable = 0;
    std::map<std::string, int> nextCounts;

    // Count next continuing preferences for ballots currently assigned to 'candidate'
    for (const auto& ballot : ballots) {
        std::string current;
        size_t curIndex = 0;
        for (; curIndex < ballot.size(); ++curIndex) {
            const auto& c = ballot[curIndex];
            if (voteCounts.find(c) != voteCounts.end()) { // mirrors transferSurplus
                current = c;
                break;
            }
        }
        if (current != candidate) continue;

        for (size_t j = curIndex + 1; j < ballot.size(); ++j) {
            const auto& nxt = ballot[j];
            if (isContinuing(nxt)) {
                ++transferable;
                nextCounts[nxt] += 1;
                break;
            }
        }
    }

    if (transferable == 0) return result;

    // Match transferSurplus: per-ballot value is floored to 2 decimals
    const double perBallot = floor2(surplus / static_cast<double>(transferable));

    // Aggregate by recipient using the same per-ballot value
    for (const auto& [rcpt, cnt] : nextCounts) {
        // Using operator<< later with fixed precision ensures clean 2-decimal output
        result[rcpt] = perBallot * static_cast<double>(cnt);
    }

    return result;
}

std::map<std::string, double> computeEliminationDistributionForLog(
    const std::string& eliminated,
    const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& elected,
    const std::map<std::string, double>& voteCounts)
{
    std::map<std::string, double> result;

    // A candidate is "continuing" if not elected and still present in voteCounts
    auto isContinuing = [&](const std::string& c) -> bool {
        return elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end();
    };

    // Collect next continuing recipients for ballots currently assigned to the eliminated candidate
    std::vector<std::string> nextRecipients;
    nextRecipients.reserve(ballots.size());

    for (const auto& ballot : ballots) {
        // Find current assignment: first candidate who is continuing OR the eliminated one
        std::string current;
        size_t curIndex = 0;
        for (; curIndex < ballot.size(); ++curIndex) {
            const auto& c = ballot[curIndex];
            if (isContinuing(c) || c == eliminated) { current = c; break; }
        }
        if (current != eliminated) continue;

        // Find next continuing preference after the eliminated candidate
        for (size_t j = curIndex + 1; j < ballot.size(); ++j) {
            const auto& c = ballot[j];
            if (isContinuing(c)) { nextRecipients.push_back(c); break; }
        }
        // If none, ballot exhausts — no transfer
    }

    const int transferable = static_cast<int>(nextRecipients.size());
    if (transferable == 0) return result;

    // Transfer each ballot at the eliminated candidate's current per-ballot value (2-decimal floor)
    const double eliminatedTotal = voteCounts.count(eliminated) ? voteCounts.at(eliminated) : 0.0;
    const double perBallot = floor2(eliminatedTotal / static_cast<double>(transferable));

    for (const auto& rcpt : nextRecipients) {
        result[rcpt] += perBallot;
    }

    return result;
}

// §11C(2): quota = validBallots / (seats + 1), rounded up to two decimals.
double calculateQuota(int validBallots, int seats)
{
    if (validBallots <= 0 || seats <= 0) return 0.0;

    // §11C(2): Regulatory threshold uplift to two decimals (ceil), not mathematical rounding.
    // quota = ceil( validBallots * 100 / (seats + 1) ) / 100
    long long numerator = static_cast<long long>(validBallots) * 100LL;
    long long denom = static_cast<long long>(seats) + 1LL;
    long long scaled = (numerator + denom - 1LL) / denom; // ceil(numerator/denom)

    return static_cast<double>(scaled) / 100.0;
}

// §11D: Tie-breaking by counting 2nd choices, then 3rd, etc. If still tied, use lot or deterministic lexical order.
// NOTE: If caller has a known earlier phase where candidates had different totals (per §11D),
// they should resolve using that history before calling this function.
std::string resolveTieSTV(const std::vector<std::vector<std::string>>& ballots,
                          const std::set<std::string>& tiedCandidates)
{
    if (tiedCandidates.empty()) return {};

    // Compute maximum ranking length present
    size_t maxLen = 0;
    for (const auto& b : ballots) maxLen = std::max(maxLen, b.size());
    if (maxLen == 0) {
        // No preferences at all; fall back to drawing lots
        return finalTiebreakPick(tiedCandidates);
    }

    // For rank = 2..maxLen, count occurrences at that exact rank.
    for (size_t rank = 2; rank <= maxLen; ++rank) {
        std::map<std::string, int> rankCounts;
        for (const auto& c : tiedCandidates) rankCounts[c] = 0;

        for (const auto& b : ballots) {
            if (b.size() < rank) continue;
            const std::string& candAtRank = b[rank - 1]; // 1-based rank -> index
            if (tiedCandidates.count(candAtRank)) {
                rankCounts[candAtRank] += 1;
            }
        }

        // Find the candidate(s) with highest count at this rank
        int best = -1;
        std::vector<std::string> bestCands;
        for (const auto& [cand, cnt] : rankCounts) {
            if (cnt > best) {
                best = cnt;
                bestCands.clear();
                bestCands.push_back(cand);
            } else if (cnt == best) {
                bestCands.push_back(cand);
            }
        }

        if (bestCands.size() == 1) {
            return bestCands.front();
        }

        // Continue to next rank if still tied
    }

    // If still tied after exhausting all ranks, resolve by drawing lots
    return finalTiebreakPick(tiedCandidates);
}

// New implementation placed near resolveTieSTV
std::string resolveTieSTVWithHistory(const std::vector<std::map<std::string, double>>& history,
                                     const std::vector<std::vector<std::string>>& ballots,
                                     const std::set<std::string>& tiedCandidates,
                                     const std::map<std::string, double>& voteCounts,
                                     const std::set<std::string>& elected)
{
    if (tiedCandidates.empty()) return {};

    // §11D primary: most recent earlier phase where totals differ (prefer highest)
    for (auto it = history.rbegin(); it != history.rend(); ++it) {
        const auto& totals = *it;

        double maxV = -std::numeric_limits<double>::infinity();
        for (const auto& cand : tiedCandidates) {
            auto f = totals.find(cand);
            const double v = (f == totals.end() ? 0.0 : f->second);
            if (v > maxV) maxV = v;
        }

        std::set<std::string> bestCands;
        for (const auto& cand : tiedCandidates) {
            auto f = totals.find(cand);
            const double v = (f == totals.end() ? 0.0 : f->second);
            if (std::abs(v - maxV) <= kEps) bestCands.insert(cand);
        }
        if (bestCands.size() == 1) return *bestCands.begin();
    }

    // §11D secondary: next continuing preferences (use MOST, not least)
    if (auto chosen = resolveTieNextContinuing(ballots, tiedCandidates, voteCounts, elected); !chosen.empty()) {
        return chosen;
    }

    return finalTiebreakPick(tiedCandidates);
}

std::map<std::string, double> countFirstChoices(const std::vector<std::vector<std::string>>& ballots)
{
    std::map<std::string, double> counts;

    // Pre-initialize all candidates observed in any position to 0.0
    for (const auto& b : ballots) {
        for (const auto& c : b) {
            if (!c.empty()) counts.emplace(c, 0.0);
        }
    }

    // Count first-choice for each ballot (if any)
    for (const auto& b : ballots) {
        if (!b.empty() && !b.front().empty()) {
            counts[b.front()] += 1.0;
        }
    }

    return counts;
}

void transferSurplus(std::string candidate,
                     double surplus,
                     const std::vector<std::vector<std::string>>& ballots,
                     std::map<std::string, double>& voteCounts,
                     const std::set<std::string>& elected)
{
    if (surplus <= 0.0) return;
    auto candIt = voteCounts.find(candidate);
    if (candIt == voteCounts.end()) return;

    auto isContinuing = [&](const std::string& c) -> bool {
        return elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end();
    };

    // Collect next continuing recipients for ballots currently assigned to 'candidate'
    std::vector<std::string> nextRecipients;
    nextRecipients.reserve(ballots.size());

    for (const auto& ballot : ballots) {
        // Find current assignment
        std::string current;
        size_t curIndex = 0;
        for (; curIndex < ballot.size(); ++curIndex) {
            const auto& c = ballot[curIndex];
            if (voteCounts.find(c) != voteCounts.end()) { current = c; break; }
        }
        if (current != candidate) continue;

        // Find next continuing preference
        for (size_t j = curIndex + 1; j < ballot.size(); ++j) {
            const auto& nxt = ballot[j];
            if (isContinuing(nxt)) { nextRecipients.push_back(nxt); break; }
        }
        // If none, ballot exhausts — do not redistribute
    }

    const int transferable = static_cast<int>(nextRecipients.size());
    if (transferable == 0) return;

    // Transfer value per transferable ballot, floored to 2 decimals
    const double transferValue = floor2(surplus / static_cast<double>(transferable));

    double totalTransferred = 0.0;
    for (const auto& rcpt : nextRecipients) {
        voteCounts[rcpt] += transferValue;
        totalTransferred += transferValue;
    }

    // Deduct transferred amount from elected candidate
    candIt->second -= totalTransferred;
    if (candIt->second < 0.0) candIt->second = 0.0;
}

void transferEliminatedVotes(std::string eliminated,
                             const std::vector<std::vector<std::string>>& ballots,
                             std::map<std::string, double>& voteCounts,
                             const std::set<std::string>& elected)
{
    if (voteCounts.find(eliminated) == voteCounts.end()) return;

    auto isContinuing = [&](const std::string& c) -> bool {
        return elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end();
    };

    // Gather next continuing recipients for ballots currently assigned to the eliminated candidate
    std::vector<std::string> nextRecipients;
    nextRecipients.reserve(ballots.size());

    for (const auto& ballot : ballots) {
        // Find current assignment: first continuing OR the eliminated one
        std::string current;
        size_t curIndex = 0;
        for (; curIndex < ballot.size(); ++curIndex) {
            const auto& c = ballot[curIndex];
            if (isContinuing(c) || c == eliminated) { current = c; break; }
        }
        if (current != eliminated) continue;

        // Find next continuing preference to receive transfer
        for (size_t j = curIndex + 1; j < ballot.size(); ++j) {
            const auto& c = ballot[j];
            if (isContinuing(c)) { nextRecipients.push_back(c); break; }
        }
        // If none, ballot exhausts
    }

    const int transferable = static_cast<int>(nextRecipients.size());
    if (transferable == 0) return;

    // Per-ballot value = eliminated total divided by number of transferable ballots, floored to 2 decimals
    const double eliminatedTotal = voteCounts[eliminated];
    const double perBallot = floor2(eliminatedTotal / static_cast<double>(transferable));

    double totalTransferred = 0.0;
    for (const auto& rcpt : nextRecipients) {
        voteCounts[rcpt] += perBallot;
        totalTransferred += perBallot;
    }

    // Note: caller erases the eliminated candidate from voteCounts.
    // Any tiny rounding remainder stays unallocated (consistent with surplus handling).
}

std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>& ballots, int seats)
{
    std::set<std::string> elected;
    if (seats <= 0 || ballots.empty()) return elected;

    // Collect all candidates from ballots
    std::set<std::string> allCandidates;
    for (const auto& b : ballots)
        for (const auto& c : b) allCandidates.insert(c);
    if (allCandidates.empty()) return elected;

    // Quota and initial first-choice tallies
    const double quota = calculateQuota(static_cast<int>(ballots.size()), seats);
    std::map<std::string, double> voteCounts = countFirstChoices(ballots);
    for (const auto& c : allCandidates)
        if (!voteCounts.count(c)) voteCounts[c] = 0.0;

    // Print quota before CSV printouts
    std::cout << "Quota," << std::fixed << std::setprecision(2) << quota << "\n";

    // Track round history for §11D tie-breaks
    std::vector<std::map<std::string, double>> history;
    history.push_back(voteCounts);

    // CSV header and initial snapshot (Round 0)
    printCsvHeader();
    int round = 0; 
    {
        std::map<std::string, double> noneAmt;
        std::map<std::string, std::vector<std::pair<std::string, double>>> noneSrc;
        printCsvRound(round++, voteCounts, elected, allCandidates, noneAmt, noneSrc);
    }

    // Helper: elect everyone meeting quota (but never exceed seats)
    auto electByQuota = [&]() -> std::vector<std::string>
    {
        std::vector<std::string> eligible;
        for (const auto& [cand, votes] : voteCounts) {
            if (!elected.count(cand) && (votes + kEps) >= quota) {
                eligible.push_back(cand);
            }
        }
        if (eligible.empty()) return {};

        const int seatsLeft = seats - static_cast<int>(elected.size());
        std::vector<std::string> newly;

        // If we have room for all, elect them all.
        if (static_cast<int>(eligible.size()) <= seatsLeft) {
            for (const auto& cand : eligible) {
                elected.insert(cand);
                newly.push_back(cand);
            }
            return newly;
        }

        // More eligible than seats left — select deterministically using vote totals and §11D tie-break.
        std::set<std::string> remainingEligible(eligible.begin(), eligible.end());
        for (int i = 0; i < seatsLeft; ++i) {
            double bestVotes = -std::numeric_limits<double>::infinity();
            for (const auto& c : remainingEligible) {
                bestVotes = std::max(bestVotes, voteCounts[c]);
            }

            std::set<std::string> tied;
            for (const auto& c : remainingEligible) {
                if (std::abs(voteCounts[c] - bestVotes) <= kEps) tied.insert(c);
            }

            std::string chosen = (tied.size() == 1)
                ? *tied.begin()
                : resolveTieSTVWithHistory(history, ballots, tied, voteCounts, elected);

            elected.insert(chosen);
            newly.push_back(chosen);
            remainingEligible.erase(chosen);
        }
        return newly;
    };

    // Helper: if remaining candidates == seats left, elect them all
    auto forceElectIfOnlyAsManyLeft = [&]()
    {
        const int seatsLeft = seats - static_cast<int>(elected.size());
        int remaining = 0;
        for (const auto& [cand, _] : voteCounts)
            if (!elected.count(cand)) ++remaining;

        if (remaining == seatsLeft) {
            for (const auto& [cand, _] : voteCounts)
                if (!elected.count(cand)) elected.insert(cand);
        }
    };

    // Main loop
    while (static_cast<int>(elected.size()) < seats) {
        // 1) Elect all who reach quota
        auto newly = electByQuota();
        if (!newly.empty()) {
            // Process surpluses in descending vote order
            std::sort(newly.begin(), newly.end(), [&](const std::string& a, const std::string& b) {
                if (std::abs(voteCounts[a] - voteCounts[b]) > kEps) return voteCounts[a] > voteCounts[b];
                return a < b;
            });

            // Accumulate a per-round transfer log (amounts received by recipients, and all sources)
            std::map<std::string, double> transferredAmounts;
            std::map<std::string, std::vector<std::pair<std::string, double>>> sourceBreakdown;

            for (const auto& cand : newly) {
                const double surplus = voteCounts[cand] - quota;
                if (surplus > 0.0) {
                    auto dist = computeSurplusDistributionForLog(cand, surplus, ballots, elected, voteCounts);
                    for (std::map<std::string, double>::const_iterator it = dist.begin(); it != dist.end(); ++it) {
                        const std::string& rcpt = it->first;
                        double amt = it->second;
                        transferredAmounts[rcpt] += amt;
                        sourceBreakdown[rcpt].push_back(std::make_pair(cand, amt));
                    }
                    transferSurplus(cand, surplus, ballots, voteCounts, elected);
                }
            }

            // Print a CSV round showing the transfers (if any)
            printCsvRound(round++, voteCounts, elected, allCandidates, transferredAmounts, sourceBreakdown);

            // Record post-surplus state for §11D history-based tie-breaks
            history.push_back(voteCounts);

            forceElectIfOnlyAsManyLeft();
            if (static_cast<int>(elected.size()) >= seats) break;

            // Loop to allow newly transferred votes to cause further elections
            continue;
        }

        // 2) No one meets quota -> eliminate lowest and transfer their votes
        bool hasCandidate = false;
        double minVotes = 0.0;
        for (const auto& [cand, votes] : voteCounts) {
            if (elected.count(cand)) continue;
            if (!hasCandidate || (votes + kEps) < minVotes) {
                hasCandidate = true;
                minVotes = votes;
            }
        }
        if (!hasCandidate) break; // no candidates left

        // Find all with the minimum vote (within epsilon)
        std::set<std::string> tied;
        for (const auto& [cand, votes] : voteCounts) {
            if (elected.count(cand)) continue;
            if (std::abs(votes - minVotes) <= kEps) tied.insert(cand);
        }

        // Tie case:
        std::string toEliminate = (tied.size() == 1)
            ? *tied.begin()
            : resolveTieSTVWithHistoryForElimination(history, ballots, tied, voteCounts, elected);

        // Compute elimination distribution for logging before mutating voteCounts
        auto elimDist = computeEliminationDistributionForLog(toEliminate, ballots, elected, voteCounts);

        // Transfer eliminated candidate's ballots to next continuing preferences
        transferEliminatedVotes(toEliminate, ballots, voteCounts, elected);

        // Remove eliminated candidate from the tally
        voteCounts.erase(toEliminate);

        // Prepare and print CSV round for elimination transfers
        {
            std::map<std::string, std::vector<std::pair<std::string, double>>> sources;
            for (const auto& kv : elimDist) {
                const std::string& rcpt = kv.first;
                double amt = kv.second;
                sources[rcpt].push_back(std::make_pair(toEliminate, amt));
            }
            printCsvRound(round++, voteCounts, elected, allCandidates, elimDist, sources);
        }

        // Record the state for future §11D resolutions
        history.push_back(voteCounts);

        forceElectIfOnlyAsManyLeft();
        if (static_cast<int>(elected.size()) >= seats) break;
        if (voteCounts.empty()) break;
    }

    // Fallback: if seats remain (e.g., due to placeholder transfers), fill by top remaining tallies
    if (static_cast<int>(elected.size()) < seats) {
        int seatsLeft = seats - static_cast<int>(elected.size());

        // Build remaining continuing candidates
        std::set<std::string> remaining;
        for (const auto& [cand, v] : voteCounts)
            if (!elected.count(cand)) remaining.insert(cand);

        // Pick winners one-by-one using §11D with history
        for (int i = 0; i < seatsLeft && !remaining.empty(); ++i) {
            double bestVotes = -std::numeric_limits<double>::infinity();
            for (const auto& c : remaining)
                bestVotes = std::max(bestVotes, voteCounts[c]);

            std::set<std::string> tied;
            for (const auto& c : remaining)
                if (std::abs(voteCounts[c] - bestVotes) <= kEps) tied.insert(c);

            std::string chosen = (tied.size() == 1)
                ? *tied.begin()
                : resolveTieSTVWithHistory(history, ballots, tied, voteCounts, elected);

            elected.insert(chosen);
            remaining.erase(chosen);
        }

        // Final snapshot (no transfers in fallback)
        std::map<std::string, double> noneAmt;
        std::map<std::string, std::vector<std::pair<std::string, double>>> noneSrc;
        printCsvRound(round++, voteCounts, elected, allCandidates, noneAmt, noneSrc);
    }

    return elected;
}

// Helper to trim leading/trailing whitespace
static std::string trim(const std::string& s)
{
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

// Display a row-numbered candidate list (1-based).
void displayCandidateList(const std::vector<std::string>& candidates)
{
    if (candidates.empty()) {
        std::cout << "No candidates.\n";
        return;
    }
    const int width = static_cast<int>(std::to_string(candidates.size()).size());
    std::cout << "Candidates:\n";
    for (size_t i = 0; i < candidates.size(); ++i) {
        std::cout << std::setw(width) << (i + 1) << ") " << candidates[i] << "\n";
    }
}

// Parse a single ranked ballot from a comma-separated string of row numbers.
std::vector<std::string> parseBallotFromRowNumbers(const std::vector<std::string>& candidates,
                                                   const std::string& line)
{
    std::vector<std::string> ballot;
    if (candidates.empty()) return ballot;

    std::set<size_t> used; // to prevent duplicates
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) {
        const std::string t = trim(token);
        if (t.empty()) continue;

        try {
            size_t idxPos = 0;
            long long val = std::stoll(t, &idxPos, 10);
            if (idxPos != t.size()) continue; // trailing characters
            if (val <= 0) continue;           // require 1-based positive
            size_t oneBased = static_cast<size_t>(val);
            if (oneBased == 0 || oneBased > candidates.size()) {
                std::cout << "Ignoring out-of-range index: " << t << "\n";
                continue;
            }
            size_t zeroBased = oneBased - 1;
            if (used.insert(zeroBased).second) {
                ballot.push_back(candidates[zeroBased]);
            } else {
                std::cout << "Ignoring duplicate index: " << t << "\n";
            }
        } catch (const std::exception&) {
            std::cout << "Ignoring invalid token: " << t << "\n";
        }
    }
    return ballot;
}

// Input multiple ballots, one per line as comma-separated row numbers.
std::vector<std::vector<std::string>> inputBallotsByRowNumbers(const std::vector<std::string>& candidates)
{
    std::vector<std::vector<std::string>> ballots;
    displayCandidateList(candidates);

    std::cout << "Enter ballots as comma-separated row numbers in preference order (e.g., 1,3,2).\n";
    std::cout << "Press Enter on an empty line to finish.\n";

    int i = 1;
    while (true) {
        std::cout << i << ": ";
        std::string line;
        if (!std::getline(std::cin, line)) break; // EOF
        const std::string s = trim(line);
        if (s.empty()) break;

        auto ballot = parseBallotFromRowNumbers(candidates, s);
        if (ballot.empty()) {
            std::cout << "(no valid selections, ballot ignored)\n";
            ++i;
            continue;
        }
        ballots.push_back(std::move(ballot));
        ++i;
    }

    std::cout << "Captured " << ballots.size() << " ballot(s).\n";
    return ballots;
}

std::vector<std::string> inputCandidateNames()
{
    std::vector<std::string> candidates;
    std::set<std::string> seen;

    std::cout << "Enter candidate names, one per line. Press Enter on an empty line to finish.\n";
    for (;;)
    {
        std::string line;
        if (!std::getline(std::cin, line)) break; // EOF
        auto name = trim(line);
        if (name.empty()) break; // end input
        if (seen.insert(name).second) {
            candidates.push_back(std::move(name));
        } else {
            std::cout << "(duplicate ignored)\n";
        }
    }

    std::cout << "Captured " << candidates.size() << " unique candidate(s).\n";
    return candidates;
}

int inputNumberOfSeats()
{
    for (;;)
    {
        std::cout << "Enter number of seats (positive integer above 1): ";
        std::string line;
        if (!std::getline(std::cin, line)) {
            std::cout << "Input stream closed.\n";
            return 0; // signal no input
        }

        const auto s = trim(line);
        if (s.empty()) {
            std::cout << "Please enter a value.\n";
            continue;
        }

        try {
            size_t idx = 0;
            long long v = std::stoll(s, &idx, 10);
            if (idx != s.size()) {
                std::cout << "Invalid characters detected. Try again.\n";
                continue;
            }
            if (v <= 1 || v > 1000000) {
                std::cout << "Seats must be between 2 and 1,000,000. Try again.\n";
                continue;
            }
            return static_cast<int>(v);
        } catch (const std::exception&) {
            std::cout << "Invalid number. Try again.\n";
        }
    }
}

int main()
{
    std::vector<std::vector<std::string>> multiSeatBallots; 
    std::vector<std::string> candidates;
    int seats;

    candidates = inputCandidateNames();
    seats = inputNumberOfSeats();
    multiSeatBallots = inputBallotsByRowNumbers(candidates);

    auto winners = runMultiSeatElection(multiSeatBallots, seats);

    std::cout << "\nElected (" << winners.size() << "): ";
    bool first = true;
    for (const auto& w : winners) {
        if (!first) std::cout << ", ";
        std::cout << w;
        first = false;
    }

    // Prompt until we get a clear Y or N
    for (;;)
    {
        std::cout << "\n\nRun a new election (Y/N)?\n";
        char again = 'N';
#ifdef _WIN32
        again = _getch(); // single key, no Enter needed
#else
        std::string line;
        if (!std::getline(std::cin, line)) return 0;
        // take first non-space char if present
        again = 0;
        for (char ch : line) {
            if (!std::isspace(static_cast<unsigned char>(ch))) { again = ch; break; }
        }
#endif
        if (std::toupper(static_cast<unsigned char>(again)) == 'Y') {
            std::cout << "\n\n";
            break; // run another election (outer loop continues)
        }
        else if (std::toupper(static_cast<unsigned char>(again)) == 'N') {
            return 0; // exit program
        }
        else {
            std::cout << "Unrecognized input.\n";
        }
    }
}

