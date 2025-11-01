// STV_Voting.cpp : Defines the entry point for the application.
//

#include "STV_Voting.h"

#include <cstdlib> // add near other includes
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
#include <cctype> // for std::toupper, std::isspace
#include <new> // for std::bad_alloc
#ifdef _WIN32
#include <conio.h>
#endif

#ifndef STV_EXPERIMENTAL_TICKETS
#define STV_EXPERIMENTAL_TICKETS 1
#endif

using namespace std;

// Numeric epsilon for FP comparisons in tie logic
constexpr double kEps = 1e-9;

#if defined(STV_EXPERIMENTAL_TICKETS) && STV_EXPERIMENTAL_TICKETS

// Move this struct definition to the top of the #if block, before any use of Ticket
struct Ticket {
    std::vector<std::string> prefs;
    size_t pos = 0;          // index of current continuing choice
    double weight = 1.0;     // rounded at each transfer with floor2
    int batchId = 0;         // incremented when rerouted; used for "last batch" logic
};

// Add a single file-scope batch id shared by helper functions and election routine
static int gBatchId = 1;

#endif // STV_EXPERIMENTAL_TICKETS

// Multi-seat Single Transferable Vote

// Step 1: Calculate the election quota
double calculateQuota(int validBallots, int seats);

// Step 2: Count first-choice votes
std::map<std::string, double> countFirstChoices(const std::vector<std::vector<std::string>>& ballots);

// Step 3: Elect candidates who meet the quota

// Step 4: Transfer surplus
void transferSurplus(std::string candidate, double surplus, double quota,
                     const std::vector<std::vector<std::string>>& ballots,
                     std::map<std::string, double>& voteCounts,
                     const std::set<std::string>& elected);

// Step 5: Eliminate the lowest-vote candidate

// Step 6: Transfer votes from the eliminated candidate
void transferEliminatedVotes(std::string eliminated, const std::vector<std::vector<std::string>>& ballots, std::map<std::string, double>& voteCounts, const std::set<std::string>& elected);

// Step 7: Resolve ties per §11D
// Add this prototype above first use (e.g., above runSingleSeatElection)
static std::string resolveTieByRawRanks(
    const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& tiedCandidates);
// Tie-breaking with round history (updated signature to accept current state)
std::string resolveTieSTVWithHistory(const std::vector<std::map<std::string, double>>& history,
                                     const std::vector<std::vector<std::string>>& ballots,
                                     const std::set<std::string>& tiedCandidates,
                                     const std::map<std::string, double>& voteCounts,
                                     const std::set<std::string>& elected);

// Add this prototype above runSingleSeatElection
static std::string resolveTieSTVWithHistoryForElimination(
    const std::vector<std::map<std::string, double>>& history,
    const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& tiedCandidates,
    const std::map<std::string, double>& voteCounts,
    const std::set<std::string>& elected);

#if defined(STV_EXPERIMENTAL_TICKETS) && STV_EXPERIMENTAL_TICKETS
std::set<std::string> runMultiSeatElection_Tickets(const std::vector<std::vector<std::string>>&, int);

// Place this prototype near the other prototypes inside the STV_EXPERIMENTAL_TICKETS block
static std::map<std::string, double> recomputeCountsFromTickets(
    const std::vector<Ticket>&,
    const std::set<std::string>&,
    const std::set<std::string>&,
    const std::set<std::string>&);

#endif

// --- Implementations ---

// Add this helper function near the top of the file, before its first use:
static double floor2(double x) { return std::floor(x * 100.0) / 100.0; }

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

// Helper to clamp displayed counts for candidates marked as Elected to the quota/threshold
static std::map<std::string, double> makeDisplayCountsClampedToQuota(
    const std::map<std::string, double>& voteCounts,
    const std::set<std::string>& electedForThisRow,
    double quota)
{
    std::map<std::string, double> display = voteCounts;
    for (const auto& c : electedForThisRow) {
        auto it = display.find(c);
        if (it != display.end() && it->second > quota) it->second = quota;
    }
    return display;
}

// Add above printCsvRound
static std::string csvQuote(const std::string& s)
{
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');
    for (char ch : s) {
        if (ch == '"') out.push_back('"'); // escape double-quote by doubling
        out.push_back(ch);
    }
    out.push_back('"');
    return out;
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
        std::cout << "," << csvQuote(cand);

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

        std::cout << ",";
        auto tIt = transferredAmounts.find(cand);
        if (tIt != transferredAmounts.end()) {
            std::cout << std::fixed << std::setprecision(2) << tIt->second;
        }

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

// Add these implementations above main() (e.g., after transferEliminatedVotes and before runMultiSeatElection)

// Recompute per-round totals for single-seat election.
// Rule: ballot goes to top continuing preference; if none remains, split 1.00
// equally among all continuing candidates, floored to 2 decimals per share.
static std::map<std::string, double> computeSingleSeatRoundTotals(
    const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& continuing)
{
    std::map<std::string, double> totals;
    for (const auto& c : continuing) totals[c] = 0.0;

    if (continuing.empty()) return totals;

    for (const auto& ballot : ballots) {
        // Find the highest-ranked continuing candidate on this ballot
        std::string top;
        for (const auto& c : ballot) {
            if (continuing.count(c)) { top = c; break; }
        }

        if (!top.empty()) {
            totals[top] += 1.0;
        } else {
            // No next preference -> split equally across all continuing candidates
            const size_t m = continuing.size();
            if (m == 0) continue;
            const double share = floor2(1.0 / static_cast<double>(m));
            if (share <= 0.0) continue; // nothing to add if share floors to 0.00
            for (const auto& c : continuing) {
                totals[c] += share;
            }
        }
    }

    return totals;
}

// Compute display-only transfers for single-seat elimination.
// Ballots whose current top is the eliminated candidate:
//  - if they have a next continuing preference, transfer 1.00 to that candidate
//  - else split 1.00 equally among all remaining continuing candidates (floored per share)
static std::map<std::string, double> computeSingleSeatEliminationDistributionForLog(
    const std::string& eliminated,
    const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& continuingBefore,
    const std::set<std::string>& continuingAfter)
{
    std::map<std::string, double> result;
    if (continuingAfter.empty()) return result;

    for (const auto& ballot : ballots) {
        // Find top continuing BEFORE elimination
        size_t elimIndex = std::numeric_limits<size_t>::max();
        std::string topBefore;
        for (size_t i = 0; i < ballot.size(); ++i) {
            if (continuingBefore.count(ballot[i])) {
                topBefore = ballot[i];
                elimIndex = i;
                break;
            }
        }
        if (topBefore != eliminated) continue;

        // Find next continuing AFTER elimination
        bool transferred = false;
        for (size_t j = elimIndex + 1; j < ballot.size(); ++j) {
            const auto& c = ballot[j];
            if (continuingAfter.count(c)) {
                result[c] += 1.0;
                transferred = true;
                break;
            }
        }
        if (!transferred) {
            // No next preference => split equally across all continuingAfter (2-decimal floor per share)
            const double share = floor2(1.0 / static_cast<double>(continuingAfter.size()));
            if (share > 0.0) {
                for (const auto& c : continuingAfter) result[c] += share;
            }
        }
    }

    return result;
}

// Single-seat election per rules §11B(1..8)
std::string runSingleSeatElection(const std::vector<std::vector<std::string>>& ballots)
{
    // Collect candidates
    std::set<std::string> continuing;
    for (const auto& b : ballots) for (const auto& c : b) if (!c.empty()) continuing.insert(c);
    if (continuing.empty()) return {};

    // Keep all candidates for consistent CSV output
    const std::set<std::string> allCandidates = continuing;

    // History of round totals for §11D tie-breaking
    std::vector<std::map<std::string, double>> history;

    // Print the majority threshold and CSV header once
    bool printedMajorityThreshold = false;
    bool printedHeader = false;

    int round = 0;

    // Round 0: initial snapshot
    {
        auto totals0 = computeSingleSeatRoundTotals(ballots, continuing);
        history.push_back(totals0);

        double sum0 = 0.0;
        for (const auto& kv : totals0) sum0 += kv.second;
        const double need0 = std::nextafter(0.5 * sum0, std::numeric_limits<double>::infinity()); // strictly greater than half

        std::string majorityCand0;
        double best0 = -std::numeric_limits<double>::infinity();
        for (const auto& [cand, v] : totals0) {
            if (v > best0) { best0 = v; majorityCand0 = cand; }
        }

        if (!printedMajorityThreshold) {
            std::cout << "MajorityThreshold," << std::fixed << std::setprecision(2) << need0 << "\n";
            printedMajorityThreshold = true;
        }
        if (!printedHeader) {
            printCsvHeader();
            printedHeader = true;
        }

        std::set<std::string> electedForDisplay0;
        if (best0 > need0) electedForDisplay0.insert(majorityCand0);

        std::map<std::string, double> noneAmt;
        std::map<std::string, std::vector<std::pair<std::string, double>>> noneSrc;
        // Remove clamping: show raw totals
        auto display0 = totals0;
        printCsvRound(round++, display0, electedForDisplay0, allCandidates, noneAmt, noneSrc);

        if (true)
        {

        } if (best0 > need0) return majorityCand0;
    }

    while (!continuing.empty()) {
        // 2+5+7: recompute totals each round following the rules
        auto totals = computeSingleSeatRoundTotals(ballots, continuing);
        history.push_back(totals);

        // Majority threshold = over half of the votes counted this round
        double sum = 0.0;
        for (const auto& kv : totals) sum += kv.second;
        const double need = std::nextafter(0.5 * sum, std::numeric_limits<double>::infinity()); // strictly greater than half

        // 3/6: check majority
        std::string majorityCand;
        double best = -std::numeric_limits<double>::infinity();
        for (const auto& [cand, v] : totals) {
            if (v > best) { best = v; majorityCand = cand; }
        }

        // Display winners for this round if any (for CSV only)
        std::set<std::string> electedForDisplay;
        if (best > need) {
            electedForDisplay.insert(majorityCand);
        }

        if (best > need) return majorityCand;

        // 4/7: eliminate the lowest-vote candidate (resolve ties per §11D elimination rules)
        double minV = std::numeric_limits<double>::infinity();
        for (const auto& [cand, v] : totals) minV = std::min(minV, v);

        std::set<std::string> tied;
        for (const auto& [cand, v] : totals) if (std::abs(v - minV) <= kEps) tied.insert(cand);

        // when computing toEliminate
        std::string toEliminate;
        if (tied.size() == 1) {
            toEliminate = *tied.begin();
        } else {
            // §11D: if no surplus transfers this phase, use raw 2nd/3rd/...
            toEliminate = resolveTieByRawRanks(ballots, tied);
            if (toEliminate.empty()) {
                // Fix: pass 'totals' and 'continuing' as voteCounts and elected
                toEliminate = resolveTieSTVWithHistoryForElimination(history, ballots, tied, totals, continuing);
            }
        }

        // Prepare transfer log for this elimination (display-only)
        std::set<std::string> continuingAfter = continuing;
        continuingAfter.erase(toEliminate);

        auto elimDist = computeSingleSeatEliminationDistributionForLog(
            toEliminate, ballots, continuing, continuingAfter);

        // Recompute totals AFTER elimination for the next round's state
        auto totalsAfter = computeSingleSeatRoundTotals(ballots, continuingAfter);

        // Combine elimination transfers with re-split delta so Transferred matches Votes delta
        std::map<std::string, double> transferredCombined = elimDist;
        std::map<std::string, std::vector<std::pair<std::string, double>>> sources;

        // Seed sources with elimination amounts
        for (const auto& kv : elimDist) {
            const std::string& rcpt = kv.first;
            double amt = kv.second;
            sources[rcpt].push_back(std::make_pair(toEliminate, amt));
        }

        // Add re-split delta (difference between after and before totals not explained by elimDist)
        for (const auto& kv : totalsAfter) {
            const std::string& cand = kv.first;
            const double afterV = kv.second;
            const double beforeV = (totals.count(cand) ? totals.at(cand) : 0.0);
            const double loggedElim = (elimDist.count(cand) ? elimDist.at(cand) : 0.0);
            const double resplit = afterV - beforeV - loggedElim;
            if (resplit > 1e-9) {
                transferredCombined[cand] += resplit;
                sources[cand].push_back(std::make_pair("Resplit", floor2(resplit)));
            }
        }

        // Optional display-only winner marking if majority achieved after elimination
        double sumAfter = 0.0;
        for (const auto& kv : totalsAfter) sumAfter += kv.second;
        const double needAfter = std::nextafter(0.5 * sumAfter, std::numeric_limits<double>::infinity());
        std::string majorityAfter;
        double bestAfter = -std::numeric_limits<double>::infinity();
        for (const auto& [cand, v] : totalsAfter) if (v > bestAfter) { bestAfter = v; majorityAfter = cand; }
        if (bestAfter > needAfter) electedForDisplay.insert(majorityAfter);

        // Show real totals this round (clamp starts next round)
        printCsvRound(round++, totalsAfter, electedForDisplay, allCandidates, transferredCombined, sources);

        // Advance state to AFTER elimination
        continuing = std::move(continuingAfter);
        history.push_back(totalsAfter);

        // If exactly one remains, they win
        if (continuing.size() == 1) return *continuing.begin();
        // If majority reached after elimination, declare winner
        if (bestAfter > needAfter) return majorityAfter;
    }

    return {};
}

// Add above first use (near other tie-break helpers)
static std::string resolveTieByRawRanks(
    const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& tiedCandidates)
{
    if (tiedCandidates.empty()) return {};
    size_t maxLen = 0;
    for (const auto& b : ballots) maxLen = std::max(maxLen, b.size());

    for (size_t rank = 2; rank <= maxLen; ++rank) {
        std::map<std::string, int> cnt;
        for (const auto& c : tiedCandidates) cnt[c] = 0;

        for (const auto& b : ballots) {
            if (b.size() < rank) continue;
            const std::string& at = b[rank - 1];
            if (tiedCandidates.count(at)) cnt[at] += 1;
        }

        int best = -1;
        std::set<std::string> bestCands;
        for (const auto& [c, v] : cnt) {
            if (v > best) { best = v; bestCands.clear(); bestCands.insert(c); }
            else if (v == best) { bestCands.insert(c); }
        }
        if (bestCands.size() == 1) return *bestCands.begin();
    }
    return {};
}

// Replace existing finalTiebreakPick with this version
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

    // Secondary: least next continuing preferences (use MOST, not least)
    if (auto chosen = resolveTieNextContinuing(ballots, tiedCandidates, voteCounts, elected); !chosen.empty()) {
        return chosen;
    }

    return finalTiebreakPick(tiedCandidates);
}

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
    int noNextCount = 0;
    std::map<std::string, int> nextCounts;

    // Count next continuing preferences and "no next" for ballots currently on 'candidate'
    for (const auto& ballot : ballots) {
        std::string current;
        size_t curIndex = 0;
        for (; curIndex < ballot.size(); ++curIndex) {
            const auto& c = ballot[curIndex];
            if (voteCounts.find(c) != voteCounts.end()) {
                current = c;
                break;
            }
        }
        if (current != candidate) continue;

        bool foundNext = false;
        for (size_t j = curIndex + 1; j < ballot.size(); ++j) {
            const auto& nxt = ballot[j];
            if (isContinuing(nxt)) {
                ++transferable;
                nextCounts[nxt] += 1;
                foundNext = true;
                break;
            }
        }
        if (!foundNext) ++noNextCount;
    }

    const int denom = transferable + noNextCount;
    if (denom == 0) return result;

    // Match transferSurplus: per-ballot value uses (transferable + noNextCount)
    const double perBallot = floor2(surplus / static_cast<double>(denom));

    // 1) Next recipients get full perBallot times their count
    for (const auto& [rcpt, cnt] : nextCounts) {
        result[rcpt] += perBallot * static_cast<double>(cnt);
    }

    // 2) Each “no next” ballot splits perBallot equally over all continuing candidates
    if (noNextCount > 0) {
        std::vector<std::string> cont;
        for (const auto& kv : voteCounts) {
            const auto& c = kv.first;
            if (elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end()) {
                cont.push_back(c);
            }
        }
        if (!cont.empty()) {
            const double split = floor2(perBallot / static_cast<double>(cont.size()));
            if (split > 0.0) {
                for (const auto& c : cont) {
                    result[c] += split * static_cast<double>(noNextCount);
                }
            }
        }
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

// New implementation
std::string resolveTieSTVWithHistory(const std::vector<std::map<std::string, double>>& history,
                                     const std::vector<std::vector<std::string>>& ballots,
                                     const std::set<std::string>& tiedCandidates,
                                     const std::map<std::string, double>& voteCounts,
                                     const std::set<std::string>& elected)
{
    if (tiedCandidates.empty()) return {};

    // §11D primary: most recent earlier phase where totals differ, preferring the highest
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

// Replace transferSurplus with this version
void transferSurplus(std::string candidate,
                     double surplus,
                     double quota,
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

    std::vector<std::string> nextRecipients;
    nextRecipients.reserve(ballots.size());

    for (const auto& ballot : ballots) {
        std::string current;
        size_t curIndex = 0;
        for (; curIndex < ballot.size(); ++curIndex) {
            const auto& c = ballot[curIndex];
            if (voteCounts.find(c) != voteCounts.end()) { current = c; break; }
        }
        if (current != candidate) continue;

        for (size_t j = curIndex + 1; j < ballot.size(); ++j) {
            const auto& nxt = ballot[j];
            if (isContinuing(nxt)) { nextRecipients.push_back(nxt); break; }
        }
    }

    const int transferable = static_cast<int>(nextRecipients.size());

    int noNextCount = 0;
    for (const auto& ballot : ballots) {
        std::string current; size_t curIndex = 0;
        for (; curIndex < ballot.size(); ++curIndex) {
            const auto& c = ballot[curIndex];
            if (voteCounts.find(c) != voteCounts.end()) { current = c; break; }
        }
        if (current != candidate) continue;

        bool hasNext = false;
        for (size_t j = curIndex + 1; j < ballot.size(); ++j) {
            const auto& nxt = ballot[j];
            if (isContinuing(nxt)) { hasNext = true; break; }
        }
        if (!hasNext) ++noNextCount;
    }

    const int denom = transferable + noNextCount;
    if (denom == 0) return;

    const double perBallot = floor2(surplus / static_cast<double>(denom));
    if (perBallot <= 0.0) return;

    // 1) next recipients get full perBallot times their count
    std::map<std::string, int> nextCounts;
    for (const auto& nxt : nextRecipients) {
        nextCounts[nxt]++;
    }

    double movedTotal = 0.0;
    for (const auto& [rcpt, cnt] : nextCounts) {
        const double add = perBallot * static_cast<double>(cnt);
        voteCounts[rcpt] += add;
        movedTotal += add;
    }

    // 2) each “no next” ballot splits perBallot equally over all continuing
    if (noNextCount > 0) {
        std::vector<std::string> cont;
        cont.reserve(voteCounts.size());
        for (const auto& kv : voteCounts) {
            if (elected.count(kv.first) == 0 && voteCounts.find(kv.first) != voteCounts.end())
                cont.push_back(kv.first);
        }
        // Remove the source candidate from cont if present
        cont.erase(std::remove(cont.begin(), cont.end(), candidate), cont.end());

        if (!cont.empty()) {
            const double split = floor2(perBallot / static_cast<double>(cont.size()));
            if (split > 0.0) {
                for (int k = 0; k < noNextCount; ++k) {
                    for (const auto& c : cont) {
                        voteCounts[c] += split;
                    }
                }
                movedTotal += split * static_cast<double>(cont.size()) * static_cast<double>(noNextCount);
            }
        }
        // If cont is empty, these ballots’ perBallot exhausts (no addition, no movedTotal increase)
    }

    // Deduct exactly what was moved (post-flooring) from the winner to conserve totals.
    candIt->second = std::max(0.0, floor2(candIt->second - movedTotal));
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

    for (const auto& rcpt : nextRecipients) {
        voteCounts[rcpt] += perBallot;
    }

    return;
}

std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>& ballots, int seats)
{
#if defined(STV_EXPERIMENTAL_TICKETS) && STV_EXPERIMENTAL_TICKETS
    return runMultiSeatElection_Tickets(ballots, seats);
#else
    std::set<std::string> elected;
    std::set<std::string> pendingSurplus; // winners elected via elimination, awaiting surplus transfer

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
        auto totals0 = voteCounts;
        std::set<std::string> electedForDisplay0;
        std::map<std::string, double> noneAmt;
        std::map<std::string, std::vector<std::pair<std::string, double>>> noneSrc;
        // Remove clamping: show raw totals
        auto display0 = totals0;
        printCsvRound(round++, display0, electedForDisplay0, allCandidates, noneAmt, noneSrc);
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

            std::string chosen;
            if (tied.size() == 1) {
                chosen = *tied.begin();
            } else {
                // §11D: use raw 2nd/3rd/... preferences first
                chosen = resolveTieByRawRanks(ballots, tied);
                if (chosen.empty()) {
                    chosen = resolveTieSTVWithHistory(history, ballots, tied, voteCounts, elected);
                }
            }

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
        // Snapshot who was elected at the start of this round (display-only)
        const std::set<std::string> electedAtStart = elected;

        // 0) Transfer any pending surpluses from winners elected in a prior elimination round
        if (!pendingSurplus.empty()) {
            std::vector<std::string> toProcess;
            toProcess.reserve(pendingSurplus.size());
            for (const auto& c : pendingSurplus) {
                if ((voteCounts[c] - quota) > kEps) toProcess.push_back(c);
            }

            if (!toProcess.empty()) {
                // Process in descending vote order, then by name
                std::sort(toProcess.begin(), toProcess.end(), [&](const std::string& a, const std::string& b) {
                    if (std::abs(voteCounts[a] - voteCounts[b]) > kEps) return voteCounts[a] > voteCounts[b];
                    return a < b;
                });

                // Aggregate per-round log
                std::map<std::string, double> transferredAmounts;
                std::map<std::string, std::vector<std::pair<std::string, double>>> sourceBreakdown;

                for (const auto& cand : toProcess) {
                    const double surplus = voteCounts[cand] - quota;
                    if (surplus <= 0.0) continue;

                    auto dist = computeSurplusDistributionForLog(cand, surplus, ballots, elected, voteCounts);
                    for (auto it = dist.begin(); it != dist.end(); ++it) {
                        const std::string& rcpt = it->first;
                        double amt = it->second;
                        transferredAmounts[rcpt] += amt;
                        sourceBreakdown[rcpt].push_back(std::make_pair(cand, amt));
                    }
                    transferSurplus(cand, surplus, quota, ballots, voteCounts, elected);
                }

                // Print a round for these surplus transfers (display-only elected snapshot)
                {
                    // Mark anyone who now reached quota as Elected (display-only)
                    std::set<std::string> electedForDisplay = elected;
                    for (const auto& kv : voteCounts) {
                        if (!elected.count(kv.first) && (kv.second + kEps) >= quota) {
                            electedForDisplay.insert(kv.first);
                        }
                    }

                    // Clamp only those already elected at the start of the round
                    // Clamp all elected to the quota for display
                    // auto displayCounts = makeDisplayCountsClampedToQuota(voteCounts, electedForDisplay, quota);
                    // Clamp only candidates elected at the start of the round; show real totals for newly crossing
                    auto displayCounts = makeDisplayCountsClampedToQuota(voteCounts, electedAtStart, quota);

                    printCsvRound(round++, displayCounts, electedForDisplay, allCandidates,
                                  transferredAmounts, sourceBreakdown);
                }

                // Record state for §11D history
                history.push_back(voteCounts);

                // Remove processed from pendingSurplus
                for (const auto& c : toProcess) pendingSurplus.erase(c);

                forceElectIfOnlyAsManyLeft();
                if (static_cast<int>(elected.size()) >= seats) break;

                // After handling pending surplus, restart loop to evaluate new elections
                continue;
            }

            // If nothing had positive surplus, clear the pending set
            pendingSurplus.clear();
        }

        // 1) Elect all who reach quota
        auto newly = electByQuota();
        if (!newly.empty()) {
            // Process surpluses in descending vote order
            std::sort(newly.begin(), newly.end(), [&](const std::string& a, const std::string& b) {
                if (std::abs(voteCounts[a] - voteCounts[b]) > kEps) return voteCounts[a] > voteCounts[b];
                return a < b;
            });

            std::map<std::string, double> transferredAmounts;
            std::map<std::string, std::vector<std::pair<std::string, double>>> sourceBreakdown;

            // in runMultiSeatElection, where surplus is processed
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
                }
                transferSurplus(cand, surplus, quota, ballots, voteCounts, elected);
            }

            // Print a CSV round showing the transfers (if any)
            // Clamp only those already elected at the start of the round; show real totals for newly crossing
            auto displayAfterElections = makeDisplayCountsClampedToQuota(voteCounts, electedAtStart, quota);
            printCsvRound(round++, displayAfterElections, elected, allCandidates, transferredAmounts, sourceBreakdown);

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

        std::string toEliminate;
        if (tied.size() == 1) {
            toEliminate = *tied.begin();
        } else {
            // §11D: if no surplus transfers this phase, use raw 2nd/3rd/...
            toEliminate = resolveTieByRawRanks(ballots, tied);
            if (toEliminate.empty()) {
                // Fix: add declaration and initialization for 'history'
                std::vector<std::map<std::string, double>> history;
                history.push_back(voteCounts);

                toEliminate = resolveTieSTVWithHistoryForElimination(history, ballots, tied, voteCounts, elected);
            }
        }

        // Compute elimination distribution for logging before mutating voteCounts
        auto elimDist = computeEliminationDistributionForLog(toEliminate, ballots, elected, voteCounts);

        // Transfer eliminated candidate's ballots to next continuing preferences
        transferEliminatedVotes(toEliminate, ballots, voteCounts, elected);

        // Remove eliminated candidate from the tally
        voteCounts.erase(toEliminate);

        // Immediately lock candidates who now meet quota due to this elimination transfer.
        // Do not transfer their surplus here; that will occur at the start of the next loop iteration.
        {
            std::vector<std::string> eligibleNow;
            for (const auto& [cand, v] : voteCounts) {
                if (!elected.count(cand) && (v + kEps) >= quota) {
                    eligibleNow.push_back(cand);
                }
            }

            if (!eligibleNow.empty()) {
                int seatsLeft = seats - static_cast<int>(elected.size());
                if (static_cast<int>(eligibleNow.size()) <= seatsLeft) {
                    for (const auto& c : eligibleNow) {
                        elected.insert(c);
                        if ((voteCounts[c] - quota) > kEps) pendingSurplus.insert(c);
                    }
                } else {
                    // More eligible than seats left — select deterministically
                    std::set<std::string> remainingEligible(eligibleNow.begin(), eligibleNow.end());
                    for (int i = 0; i < seatsLeft && !remainingEligible.empty(); ++i) {
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
                        if ((voteCounts[chosen] - quota) > kEps) pendingSurplus.insert(chosen);
                        remainingEligible.erase(chosen);
                    }
                }
            }
        }

        // Prepare and print CSV round for elimination transfers
        {
            std::map<std::string, std::vector<std::pair<std::string, double>>> sources;
            for (const auto& kv : elimDist) {
                const std::string& rcpt = kv.first;
                double amt = kv.second;
                sources[rcpt].push_back(std::make_pair(toEliminate, amt));
            }
            // Prepare and print CSV round for elimination transfers
            // Clamp only those already elected at the start of the round; show real totals for newly crossing
            auto displayCounts = makeDisplayCountsClampedToQuota(voteCounts, electedAtStart, quota);
            printCsvRound(round++, displayCounts, elected, allCandidates, elimDist, sources);
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
        auto displayFinal = makeDisplayCountsClampedToQuota(voteCounts, elected, quota);
        printCsvRound(round++, displayFinal, elected, allCandidates, noneAmt, noneSrc);
    }

    return elected;
#endif
}

// Add this implementation near the end of the file, after the #endif for STV_EXPERIMENTAL_TICKETS

#if defined(STV_EXPERIMENTAL_TICKETS) && STV_EXPERIMENTAL_TICKETS

// Add this function prototype above its first use in runMultiSeatElection_Tickets (inside the #if STV_EXPERIMENTAL_TICKETS block):

static std::pair<std::map<std::string, double>,
    std::map<std::string, std::vector<std::pair<std::string, double>>>>
    transferSurplusLastBatch(const std::string& cand,
        double quota,
        std::vector<Ticket>& tickets,
        const std::set<std::string>& allCandidates,
        const std::set<std::string>& elected,
        const std::set<std::string>& eliminated);

// Add this function prototype above its first use in runMultiSeatElection_Tickets (inside the #if STV_EXPERIMENTAL_TICKETS block):

static std::pair<std::map<std::string, double>,
    std::map<std::string, std::vector<std::pair<std::string, double>>>>
    eliminateAndTransferFull(const std::string& eliminatedCand,
        std::vector<Ticket>& tickets,
        const std::set<std::string>& elected,
        const std::set<std::string>& eliminated);

// Add this function above its first use in runMultiSeatElection_Tickets (inside the #if STV_EXPERIMENTAL_TICKETS block):

static std::map<std::string, double> recomputeCountsFromTickets(
    const std::vector<Ticket>& tickets,
    const std::set<std::string>& allCandidates,
    const std::set<std::string>& elected,
    const std::set<std::string>& eliminated)
{
    std::map<std::string, double> counts;
    // Initialize all non-eliminated candidates (including elected) so they appear in display.
    for (const auto& c : allCandidates) {
        if (eliminated.count(c) == 0) counts[c] = 0.0;
    }

    // Sum weights for all non-eliminated candidates (do NOT exclude elected here).
    for (const auto& t : tickets) {
        if (t.weight <= 0.0 || t.pos >= t.prefs.size()) continue;
        const std::string& c = t.prefs[t.pos];
        if (eliminated.count(c) == 0) {
            counts[c] += t.weight;
        }
    }
    return counts;
}

std::set<std::string> runMultiSeatElection_Tickets(const std::vector<std::vector<std::string>>& ballots, int seats)
{
    // Basic ticket-based STV implementation (experimental)
    std::set<std::string> elected;
    std::set<std::string> eliminated;
    std::set<std::string> allCandidates;
    std::vector<Ticket> tickets;

    // Collect all candidates
    for (const auto& b : ballots)
        for (const auto& c : b)
            if (!c.empty()) allCandidates.insert(c);

    // Initialize tickets
    for (const auto& b : ballots) {
        if (!b.empty()) {
            Ticket t;
            t.prefs = b;
            t.pos = 0;
            t.weight = 1.0;
            t.batchId = gBatchId++;
            tickets.push_back(t);
        }
    }

    if (allCandidates.empty() || seats <= 0) return elected;

    // Calculate quota
    const double quota = calculateQuota(static_cast<int>(ballots.size()), seats);
    std::cout << "Quota," << std::fixed << std::setprecision(2) << quota << "\n";
    printCsvHeader();
    int round = 0;

    // Initial counts
    auto voteCounts = recomputeCountsFromTickets(tickets, allCandidates, elected, eliminated);
    {
        auto totals0 = voteCounts;
        std::set<std::string> electedForDisplay0;
        std::map<std::string, double> noneAmt;
        std::map<std::string, std::vector<std::pair<std::string, double>>> noneSrc;
        // Remove clamping: show raw totals
        auto display0 = totals0;
        printCsvRound(round++, display0, electedForDisplay0, allCandidates, noneAmt, noneSrc);
    }

    // Main loop
    while (static_cast<int>(elected.size()) < seats) {
        // 1) Elect all who reach quota
        std::vector<std::string> eligible;
        for (const auto& [cand, votes] : voteCounts) {
            if (!elected.count(cand) && (votes + kEps) >= quota) {
                eligible.push_back(cand);
            }
        }
        if (!eligible.empty()) {
            // Sort by votes descending, then name
            std::sort(eligible.begin(), eligible.end(), [&](const std::string& a, const std::string& b) {
                if (std::abs(voteCounts[a] - voteCounts[b]) > kEps) return voteCounts[a] > voteCounts[b];
                return a < b;
                });

            std::map<std::string, double> transferredAmounts;
            std::map<std::string, std::vector<std::pair<std::string, double>>> sourceBreakdown;

            for (const auto& cand : eligible) {
                elected.insert(cand);
                double surplus = voteCounts[cand] - quota;
                if (surplus > 0.0) {
                    auto [transferred, sources] = transferSurplusLastBatch(cand, quota, tickets, allCandidates, elected, eliminated);
                    for (const auto& [rcpt, amt] : transferred) {
                        transferredAmounts[rcpt] += amt;
                    }
                    for (const auto& [rcpt, srcs] : sources) {
                        for (const auto& src : srcs) {
                            sourceBreakdown[rcpt].push_back(src);
                        }
                    }
                }
            }

            voteCounts = recomputeCountsFromTickets(tickets, allCandidates, elected, eliminated);
            auto displayCounts = makeDisplayCountsClampedToQuota(voteCounts, elected, quota);
            printCsvRound(round++, displayCounts, elected, allCandidates, transferredAmounts, sourceBreakdown);

            // If enough elected, break
            if (static_cast<int>(elected.size()) >= seats) break;
            continue;
        }

        // 2) Eliminate lowest
        double minVotes = std::numeric_limits<double>::infinity();
        for (const auto& [cand, votes] : voteCounts) {
            if (elected.count(cand)) continue;
            if (!minVotes || (votes + kEps) < minVotes) {
                minVotes = votes;
            }
        }
        if (!minVotes) break; // no candidates left

        // Find all with the minimum vote (within epsilon)
        std::set<std::string> tied;
        for (const auto& [cand, votes] : voteCounts) {
            if (elected.count(cand)) continue;
            if (std::abs(votes - minVotes) <= kEps) tied.insert(cand);
        }

        std::string toEliminate;
        if (tied.size() == 1) {
            toEliminate = *tied.begin();
        }
        else {
            // §11D: if no surplus transfers this phase, use raw 2nd/3rd/...
            toEliminate = resolveTieByRawRanks(ballots, tied);
            if (toEliminate.empty()) {
                // Fix: add declaration and initialization for 'history'
                std::vector<std::map<std::string, double>> history;
                history.push_back(voteCounts);

                toEliminate = resolveTieSTVWithHistoryForElimination(history, ballots, tied, voteCounts, elected);
            }
        }

        auto [transferred, sources] = eliminateAndTransferFull(toEliminate, tickets, elected, eliminated);
        voteCounts = recomputeCountsFromTickets(tickets, allCandidates, elected, eliminated);
        auto displayCounts = makeDisplayCountsClampedToQuota(voteCounts, elected, quota);
        printCsvRound(round++, displayCounts, elected, allCandidates, transferred, sources);

        // If only as many left as seats, elect all
        int seatsLeft = seats - static_cast<int>(elected.size());
        int remaining = 0;
        for (const auto& c : allCandidates)
            if (!elected.count(c) && !eliminated.count(c)) ++remaining;
        if (remaining == seatsLeft) {
            for (const auto& c : allCandidates)
                if (!elected.count(c) && !eliminated.count(c)) elected.insert(c);
            break;
        }
    }

    return elected;
}

#endif // STV_EXPERIMENTAL_TICKETS

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

// Validating parser: collects errors so the caller can re-prompt
struct ParseBallotResult {
    std::vector<std::string> ballot;
    bool hadOutOfRange = false;
    bool hadDuplicate = false;
    bool hadInvalid = false;
    bool hasError() const { return hadOutOfRange || hadDuplicate || hadInvalid; }
};

static ParseBallotResult parseBallotFromRowNumbersWithValidation(
    const std::vector<std::string>& candidates,
    const std::string& line)
{
    ParseBallotResult res;
    if (candidates.empty()) return res;

    std::set<size_t> used;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) {
        const std::string t = trim(token);
        if (t.empty()) continue;

        try {
            size_t idxPos = 0;
            long long val = std::stoll(t, &idxPos, 10);
            if (idxPos != t.size() || val <= 0) {
                std::cout << "Invalid token: " << t << "\n";
                res.hadInvalid = true;
                continue;
            }
            size_t oneBased = static_cast<size_t>(val);
            if (oneBased == 0 || oneBased > candidates.size()) {
                std::cout << "Out-of-range index: " << t << "\n";
                res.hadOutOfRange = true;
                continue;
            }
            size_t zeroBased = oneBased - 1;
            if (used.insert(zeroBased).second) {
                res.ballot.push_back(candidates[zeroBased]);
            } else {
                std::cout << "Duplicate index: " << t << "\n";
                res.hadDuplicate = true;
            }
        } catch (const std::exception&) {
            std::cout << "Invalid token: " << t << "\n";
            res.hadInvalid = true;
        }
    }
    return res;
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

        ParseBallotResult parsed = parseBallotFromRowNumbersWithValidation(candidates, s);

        if (parsed.hasError()) {
            std::cout << "(please re-enter ballot " << i << ": ";
            bool first = true;
            if (parsed.hadOutOfRange) { std::cout << (first ? "" : ", ") << "out-of-range index"; first = false; }
            if (parsed.hadDuplicate)  { std::cout << (first ? "" : ", ") << "duplicate index";   first = false; }
            if (parsed.hadInvalid)    { std::cout << (first ? "" : ", ") << "invalid token";     first = false; }
            std::cout << " detected)\n";
            continue; // re-prompt same ballot number
        }

        if (parsed.ballot.empty()) {
            std::cout << "(no valid selections, ballot ignored)\n";
            ++i;
            continue;
        }

        ballots.push_back(std::move(parsed.ballot));
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
        std::cout << "Enter number of seats (positive integer): ";
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
            if (v <= 0 || v > 1000000) {
                std::cout << "Seats must be between 1 and 1,000,000. Try again.\n";
                continue;
            }
            return static_cast<int>(v);
        } catch (const std::exception&) {
            std::cout << "Invalid number. Try again.\n";
        }
    }
}


#if defined(STV_EXPERIMENTAL_TICKETS) && STV_EXPERIMENTAL_TICKETS

static bool isContinuingCandidate(const std::string& c,
    const std::set<std::string>& elected,
    const std::set<std::string>& eliminated)
{
    return elected.count(c) == 0 && eliminated.count(c) == 0;
}

static void advanceToNextContinuing(Ticket& t,
    const std::set<std::string>& elected,
    const std::set<std::string>& eliminated)
{
    while (t.pos < t.prefs.size()) {
        const std::string& c = t.prefs[t.pos];
        if (isContinuingCandidate(c, elected, eliminated)) return;
        ++t.pos;
    }
}

static std::map<std::string, std::vector<int>> buildAssignedMap(
    const std::vector<Ticket>& tickets,
    const std::set<std::string>& elected,
    const std::set<std::string>& eliminated)
{
    std::map<std::string, std::vector<int>> assigned;
    for (int i = 0; i < static_cast<int>(tickets.size()); ++i) {
        const auto& t = tickets[i];
        if (t.weight <= 0.0 || t.pos >= t.prefs.size()) continue;
        const std::string& c = t.prefs[t.pos];
        if (isContinuingCandidate(c, elected, eliminated)) {
            assigned[c].push_back(i);
        }
    }
    return assigned;
}

// Transfer elected candidate's surplus using only the "last equal-valued batch" that caused the surplus.
// Returns per-recipient amounts and source breakdown for CSV.
static std::pair<std::map<std::string, double>,
    std::map<std::string, std::vector<std::pair<std::string, double>>>>
    transferSurplusLastBatch(const std::string& cand,
        double quota,
        std::vector<Ticket>& tickets,
        const std::set<std::string>& allCandidates,
        const std::set<std::string>& elected,
        const std::set<std::string>& eliminated)
{
    std::map<std::string, double> transferred;
    std::map<std::string, std::vector<std::pair<std::string, double>>> sources;

    // Compute current total for cand
    double total = 0.0;
    for (const auto& t : tickets) {
        if (t.weight <= 0.0 || t.pos >= t.prefs.size()) continue;
        if (t.prefs[t.pos] == cand) total += t.weight;
    }
    double surplus = total - quota;
    if (surplus <= kEps) return { transferred, sources };

    // Identify last equal-valued batch among tickets currently at cand
    int maxBatch = std::numeric_limits<int>::min();
    for (const auto& t : tickets) {
        if (t.weight <= 0.0 || t.pos >= t.prefs.size()) continue;
        if (t.prefs[t.pos] == cand) maxBatch = std::max(maxBatch, t.batchId);
    }
    if (maxBatch == std::numeric_limits<int>::min()) return { transferred, sources };

    std::vector<int> lastBatchIdx;
    lastBatchIdx.reserve(tickets.size());
    for (int i = 0; i < static_cast<int>(tickets.size()); ++i) {
        const auto& t = tickets[i];
        if (t.weight <= 0.0 || t.pos >= t.prefs.size()) continue;
        if (t.prefs[t.pos] == cand && t.batchId == maxBatch) {
            lastBatchIdx.push_back(i);
        }
    }
    if (lastBatchIdx.empty()) return { transferred, sources };

    // Determine next recipients (or no-next) for last batch
    int transferable = 0;
    int noNextCount = 0;
    std::vector<int> nextPos(lastBatchIdx.size(), -1); // next continuing position index in prefs
    for (size_t k = 0; k < lastBatchIdx.size(); ++k) {
        const int idx = lastBatchIdx[k];
        const auto& t = tickets[idx];
        size_t p = t.pos + 1;
        int foundPos = -1;
        for (; p < t.prefs.size(); ++p) {
            const std::string& nxt = t.prefs[p];
            if (isContinuingCandidate(nxt, elected, eliminated)) {
                foundPos = static_cast<int>(p);
                break;
            }
        }
        if (foundPos >= 0) { transferable++; nextPos[k] = foundPos; }
        else { noNextCount++; }
    }

    int denom = transferable + noNextCount;
    if (denom == 0) return { transferred, sources };

    const double perBallot = floor2(surplus / static_cast<double>(denom));
    if (perBallot <= 0.0) return { transferred, sources };

    const int thisBatchId = gBatchId++;
    // For "no next" case, precompute continuing recipients (exclude elected/eliminated and also exclude cand)
    std::vector<std::string> cont;
    cont.reserve(allCandidates.size());
    for (const auto& c : allCandidates) {
        if (isContinuingCandidate(c, elected, eliminated) && c != cand) cont.push_back(c);
    }

    // Apply transfers: split each ticket's perBallot out of the candidate's ticket into new ticket(s).
    for (size_t k = 0; k < lastBatchIdx.size(); ++k) {
        const int idx = lastBatchIdx[k];
        Ticket& t = tickets[idx];

        // Reduce original ticket by perBallot (floor to 2 decimals for the ticket weight)
        t.weight = floor2(t.weight - perBallot);
        if (t.weight < 0.0) t.weight = 0.0; // safety clamp

        if (nextPos[k] >= 0) {
            // Create a new ticket fragment for the recipient
            Ticket fragment = t;
            fragment.pos = static_cast<size_t>(nextPos[k]);
            fragment.weight = perBallot;
            fragment.batchId = thisBatchId;
            tickets.push_back(fragment);

            const std::string& rcpt = fragment.prefs[fragment.pos];
            transferred[rcpt] += perBallot;
            sources[rcpt].push_back(std::make_pair(cand, perBallot));
        }
        else {
            // No next -> split equally across all continuing candidates
            if (!cont.empty()) {
                const double split = floor2(perBallot / static_cast<double>(cont.size()));
                if (split > 0.0) {
                    for (const auto& rcpt : cont) {
                        Ticket fragment = t;
                        // Move fragment to rcpt (find its position in prefs or synthesize)
                        size_t p = 0;
                        bool found = false;
                        for (; p < fragment.prefs.size(); ++p) {
                            if (fragment.prefs[p] == rcpt) { found = true; break; }
                        }
                        if (!found) {
                            fragment.prefs.push_back(rcpt);
                            p = fragment.prefs.size() - 1;
                        }
                        fragment.pos = p;
                        fragment.weight = split;
                        fragment.batchId = thisBatchId;
                        tickets.push_back(fragment);

                        transferred[rcpt] += split;
                        sources[rcpt].push_back(std::make_pair(cand, split));
                    }
                }
            }
            // If cont is empty, perBallot exhausts (no addition, no movedTotal increase)
        }
    }

    return { transferred, sources };
}

// Eliminate candidate and transfer each ticket at full current weight to the next continuing preference.
// Returns per-recipient amounts and source breakdown for CSV.
static std::pair<std::map<std::string, double>,
    std::map<std::string, std::vector<std::pair<std::string, double>>>>
    eliminateAndTransferFull(const std::string& eliminatedCand,
        std::vector<Ticket>& tickets,
        const std::set<std::string>& elected,
        const std::set<std::string>& eliminated)
{
    std::map<std::string, double> transferred;
    std::map<std::string, std::vector<std::pair<std::string, double>>> sources;

    const int thisBatchId = gBatchId++;

    for (auto& t : tickets) {
        if (t.weight <= 0.0 || t.pos >= t.prefs.size()) continue;
        if (t.prefs[t.pos] != eliminatedCand) continue;

        // Find next continuing recipient
        size_t p = t.pos + 1;
        int nextFoundPos = -1;
        for (; p < t.prefs.size(); ++p) {
            const std::string& nxt = t.prefs[p];
            if (isContinuingCandidate(nxt, elected, eliminated)) { nextFoundPos = static_cast<int>(p); break; }
        }
        if (nextFoundPos >= 0) {
            const std::string& rcpt = t.prefs[static_cast<size_t>(nextFoundPos)];
            transferred[rcpt] += t.weight;
            sources[rcpt].push_back(std::make_pair(eliminatedCand, t.weight));

            // Move ticket in place to recipient with full weight
            t.pos = static_cast<size_t>(nextFoundPos);
            t.batchId = thisBatchId;
        }
        else {
            // Exhausted
            t.weight = 0.0;
        }
    }

    return { transferred, sources };
}

#endif // STV_EXPERIMENTAL_TICKETS

#ifndef STV_NO_MAIN
int main()
{
    try {
        for (;;)
        {
            std::vector<std::vector<std::string>> singleSeatBallots;
            std::vector<std::vector<std::string>> multiSeatBallots;
            std::vector<std::string> candidates;
            int seats;

            candidates = inputCandidateNames();
            seats = inputNumberOfSeats();
            
            if (seats == 1) {
                singleSeatBallots = inputBallotsByRowNumbers(candidates);

                auto winner = runSingleSeatElection(singleSeatBallots);
                if (!winner.empty()) {
                    std::cout << "\nElected (1): " << winner;
                } else {
                    std::cout << "\nNo winner could be determined.";
                }
            }
            else {
                multiSeatBallots = inputBallotsByRowNumbers(candidates);

                auto winners = runMultiSeatElection(multiSeatBallots, seats);

                std::cout << "\nElected (" << winners.size() << "): ";
                bool first = true;
                for (const auto& w : winners) {
                    if (!first) std::cout << ", ";
                    std::cout << w;
                    first = false;
                }
            }

            // Prompt until we get a clear Y or N
            bool runAgain = false;
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
                    bool runAgainInner = true; // to satisfy compiler warnings
                    runAgain = runAgainInner;
                    break; // rerun the entire workflow
                }
                else if (std::toupper(static_cast<unsigned char>(again)) == 'N') {
                    return 0; // exit program
                }
                else {
                    std::cout << "Unrecognized input.\n";
                }
            }

            if (!runAgain) break; // safety
        }
    } catch (const std::bad_alloc& ex) {
        std::cerr << "Error: memory allocation failed (std::bad_alloc). " << ex.what() << "\n";
        return EXIT_FAILURE;
    } catch (const std::exception& ex) {
        std::cerr << "Unhandled exception: " << ex.what() << "\n";
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unhandled unknown exception.\n";
        return EXIT_FAILURE;
    }
}
#endif

