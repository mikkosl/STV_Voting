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
#include <sstream> // added
#ifdef _WIN32
#include <conio.h>
#endif

using namespace std;

// Multi-seat Single Transferable Vote

// Step 1: Calculate the election quota
double calculateQuota(int validBallots, int seats);

// Step 2: Count first-choice votes
std::map<std::string, double> countFirstChoices(const std::vector<std::vector<std::string>>& ballots);

// Step 3: Elect candidates who meet the quota
std::set<std::string> electCandidates(const std::map<std::string, double>& voteCounts, double quota);

// Step 4: Transfer surplus
void transferSurplus(std::string candidate, double surplus, const std::vector<std::vector<std::string>>& ballots, std::map<std::string, double>& voteCounts, const std::set<std::string>& elected);

// Step 5: Eliminate the lowest-vote candidate
std::string eliminateLowestCandidate(const std::map<std::string, double>& voteCounts, const std::set<std::string>& elected);

// Step 6: Transfer votes from the eliminated candidate
void transferEliminatedVotes(std::string eliminated, const std::vector<std::vector<std::string>>& ballots, std::map<std::string, double>& voteCounts, const std::set<std::string>& elected);

// Step 7: Resolve ties per §11D
std::string resolveTieSTV(const std::vector<std::vector<std::string>>& ballots, const std::set<std::string>& tiedCandidates);

// Tie-breaking with round history
std::string resolveTieSTVWithHistory(const std::vector<std::map<std::string, double>>& history,
                                     const std::vector<std::vector<std::string>>& ballots,
                                     const std::set<std::string>& tiedCandidates);

// Main loop: repeat until all seats are filled
std::set<std::string> runMultiSeatElection(const std::vector<std::vector<std::string>>& ballots, int seats);

std::map<std::string, double> computeSurplusDistributionForLog(
    const std::string& candidate,
    double surplus,
    const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& elected,
    const std::map<std::string, double>& voteCounts);

// --- Implementations ---

// Print CSV header for election rounds
void printCsvHeader()
{
    std::cout << "Round";
    std::cout << ",Candidate";
    std::cout << ",Votes";
    std::cout << ",Status";
    std::cout << ",Transferred";
    std::cout << ",Source";
    std::cout << "\n";
}

// Print a CSV row for each candidate in the round, including transfers and status.
void printCsvRound(
    int round,
    const std::map<std::string, double>& voteCounts,
    const std::set<std::string>& elected,
    const std::set<std::string>& allCandidates,
    const std::map<std::string, double>& transferredAmounts,
    const std::map<std::string, std::string>& recipients)
{
    for (const auto& cand : allCandidates) {
        std::cout << round;
        std::cout << "," << cand;
        auto it = voteCounts.find(cand);
        double votes = (it != voteCounts.end()) ? it->second : 0.0;
        std::cout << "," << std::fixed << std::setprecision(2) << votes;

        // Existing elected flag (kept for status calculation only; not printed)
        const bool isElected = elected.count(cand) != 0;
        // std::cout << "," << (isElected ? "Y" : ""); // removed from CSV

        // Status column
        std::string status;
        if (isElected) {
            status = "Elected";
        } else if (voteCounts.find(cand) == voteCounts.end()) {
            status = "Eliminated";
        } else {
            status = "Continuing";
        }
        std::cout << "," << status;

        // Transfer info
        auto tIt = transferredAmounts.find(cand);
        std::cout << "," << (tIt != transferredAmounts.end() ? std::to_string(tIt->second) : "");
        auto rIt = recipients.find(cand);
        std::cout << "," << (rIt != recipients.end() ? rIt->second : "");
        std::cout << "\n";
    }
}


std::map<std::string, double> computeEliminationDistributionForLog(
    const std::string& eliminated,
    const std::vector<std::vector<std::string>>& ballots,
    const std::set<std::string>& elected,
    const std::map<std::string, double>& voteCounts)
{
    std::map<std::string, double> result;

    auto isContinuing = [&](const std::string& c) -> bool {
        return elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end();
    };

    for (const auto& ballot : ballots) {
        // Find current assignment (first continuing candidate on this ballot)
        std::string current;
        size_t curIndex = 0;
        for (; curIndex < ballot.size(); ++curIndex) {
            const auto& c = ballot[curIndex];
            if (isContinuing(c)) { current = c; break; }
        }

        // Only transfer ballots currently assigned to the eliminated candidate
        if (current != eliminated) continue;

        // Find next continuing preference after the eliminated candidate on this ballot
        std::string next;
        for (size_t j = curIndex + 1; j < ballot.size(); ++j) {
            const auto& c = ballot[j];
            if (isContinuing(c)) { next = c; break; }
        }

        // Transfer one full vote to the next continuing candidate if present
        if (!next.empty()) {
            result[next] += 1.0;
        }
        // else: ballot exhausts; no transfer
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

// §11D: Tie-breaking by counting 2nd choices, then 3rd, etc. If still tied, use deterministic lexical order.
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
        // No preferences at all; fall back to deterministic lexical order
        return *tiedCandidates.begin();
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

    // If still tied after exhausting all ranks, resolve deterministically (lexical order).
    return *tiedCandidates.begin();
}

// New implementation placed near resolveTieSTV
std::string resolveTieSTVWithHistory(const std::vector<std::map<std::string, double>>& history,
                                     const std::vector<std::vector<std::string>>& ballots,
                                     const std::set<std::string>& tiedCandidates)
{
    if (tiedCandidates.empty()) return {};

    // §11D first: use the most recent earlier phase where totals differ
    for (auto it = history.rbegin(); it != history.rend(); ++it) {
        const auto& totals = *it;

        double best = -std::numeric_limits<double>::infinity();
        std::vector<std::string> bestCands;
        for (const auto& cand : tiedCandidates) {
            auto f = totals.find(cand);
            const double v = (f == totals.end() ? 0.0 : f->second);
            if (v > best) {
                best = v;
                bestCands.clear();
                bestCands.push_back(cand);
            } else if (v == best) {
                bestCands.push_back(cand);
            }
        }
        if (bestCands.size() == 1) return bestCands.front();
        // else continue to earlier phase
    }

    // Fallback to your existing rank-based tie-breaker
    return resolveTieSTV(ballots, tiedCandidates);
}

// Placeholder implementations to keep linkage intact; core STV loop to be implemented.
std::map<std::string, double> countFirstChoices(const std::vector<std::vector<std::string>>& ballots)
{
    std::map<std::string, double> counts;
    for (const auto& b : ballots) {
        if (!b.empty()) counts[b.front()] += 1.0;
    }
    return counts;
}

std::set<std::string> electCandidates(const std::map<std::string, double>& voteCounts, double quota)
{
    std::set<std::string> elected;
    for (const auto& [cand, votes] : voteCounts) {
        if (votes >= quota) elected.insert(cand);
    }
    return elected;
}

void transferSurplus(std::string candidate,
                     double surplus,
                     const std::vector<std::vector<std::string>>& ballots,
                     std::map<std::string, double>& voteCounts,
                     const std::set<std::string>& elected)
{
    // §11 surplus transfer (Inclusive, unweighted per transferable ballot):
    // Distribute exactly the surplus across ballots currently assigned to the elected
    // candidate that have a next continuing preference. Each such ballot transfers an equal
    // fraction: transferValue = surplus / transferableBallotsCount.
    //
    // This ensures the total transferred equals the surplus (subject to tiny FP rounding),
    // and the elected candidate's tally is reduced by that same amount.

    if (surplus <= 0.0) return;

    auto candIt = voteCounts.find(candidate);
    if (candIt == voteCounts.end()) return;

    auto isContinuing = [&](const std::string& c) -> bool {
        // Continuing = not elected and still present in voteCounts (i.e., not already eliminated)
        return elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end();
    };

    // Count only ballots currently assigned to 'candidate' that also have a next continuing preference.
    int transferable = 0;
    std::map<std::string, int> nextCounts;

    for (const auto& ballot : ballots) {
        // Current assignment = first candidate still in the tally map (could be elected)
        std::string current;
        size_t curIndex = 0;
        for (; curIndex < ballot.size(); ++curIndex) {
            const auto& c = ballot[curIndex];
            if (voteCounts.find(c) != voteCounts.end()) { // still in the tally map
                current = c;
                break;
            }
        }
        if (current != candidate) continue; // not assigned to the elected candidate

        // Find next continuing (not elected, still in count) after the elected candidate
        for (size_t j = curIndex + 1; j < ballot.size(); ++j) {
            const auto& nxt = ballot[j];
            if (isContinuing(nxt)) {
                ++transferable;
                nextCounts[nxt] += 1;
                break;
            }
        }
        // If none found, ballot exhausts; it does not participate in surplus transfer.
    }

    if (transferable == 0) {
        // Nothing to transfer under §11 (no next preferences on assigned ballots).
        return;
    }

    // Each transferable ballot moves this fraction of a vote
    const double transferValue = surplus / static_cast<double>(transferable);

    // Apply transfers; track total to conserve mass
    double totalTransferred = 0.0;
    for (const auto& [rcpt, cnt] : nextCounts) {
        const double amt = transferValue * static_cast<double>(cnt);
        voteCounts[rcpt] += amt;
        totalTransferred += amt;
    }

    // Deduct the amount moved from the elected candidate
    candIt->second -= totalTransferred;
    if (candIt->second < 0.0) candIt->second = 0.0; // guard against tiny rounding negatives
}

std::string eliminateLowestCandidate(const std::map<std::string, double>& voteCounts, const std::set<std::string>& elected)
{
    std::string lowest;
    double minVotes = std::numeric_limits<double>::infinity();
    for (const auto& [cand, votes] : voteCounts) {
        if (elected.count(cand)) continue;
        if (votes < minVotes || (votes == minVotes && cand < lowest)) {
            minVotes = votes;
            lowest = cand;
        }
    }
    return lowest;
}

void transferEliminatedVotes(std::string eliminated,
                             const std::vector<std::vector<std::string>>& ballots,
                             std::map<std::string, double>& voteCounts,
                             const std::set<std::string>& elected)
{
    // If the eliminated candidate is not in the current count, nothing to do.
    if (voteCounts.find(eliminated) == voteCounts.end()) return;

    auto isContinuing = [&](const std::string& c) -> bool {
        // Continuing = not elected and still present in voteCounts (i.e., not already eliminated)
        return elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end();
    };

    for (const auto& ballot : ballots) {
        // Find current assignment (first continuing candidate on this ballot)
        std::string current;
        size_t curIndex = 0;
        for (; curIndex < ballot.size(); ++curIndex) {
            const auto& c = ballot[curIndex];
            if (isContinuing(c)) { current = c; break; }
        }

        // Only transfer ballots currently assigned to the eliminated candidate
        if (current != eliminated) continue;

        // Find next continuing preference after the eliminated candidate on this ballot
        std::string next;
        for (size_t j = curIndex + 1; j < ballot.size(); ++j) {
            const auto& c = ballot[j];
            if (isContinuing(c)) { next = c; break; }
        }

        // Transfer one full vote to the next continuing candidate if present
        if (!next.empty()) {
            voteCounts[next] += 1.0;
        }
        // else: ballot exhausts; no transfer
    }
    // Caller removes the eliminated candidate from voteCounts after this.
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
        std::map<std::string, std::string> noneSrc;
        printCsvRound(round++, voteCounts, elected, allCandidates, noneAmt, noneSrc);
    }

    // Helper: elect everyone meeting quota (but never exceed seats)
    auto electByQuota = [&]() -> std::vector<std::string>
    {
        std::vector<std::string> eligible;
        for (const auto& [cand, votes] : voteCounts) {
            if (!elected.count(cand) && votes >= quota) {
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
                if (voteCounts[c] == bestVotes) tied.insert(c);
            }

            std::string chosen = (tied.size() == 1)
                ? *tied.begin()
                : resolveTieSTVWithHistory(history, ballots, tied);

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
                if (voteCounts[a] != voteCounts[b]) return voteCounts[a] > voteCounts[b];
                return a < b;
            });

            // Accumulate a per-round transfer log (amounts received by recipients, and source)
            std::map<std::string, double> transferredAmounts;
            std::map<std::string, std::string> recipients;

            for (const auto& cand : newly) {
                const double surplus = voteCounts[cand] - quota;
                if (surplus > 0.0) {
                    // Compute distribution for logging, then perform the actual transfer.
                    auto dist = computeSurplusDistributionForLog(cand, surplus, ballots, elected, voteCounts);
                    for (std::map<std::string, double>::const_iterator it = dist.begin(); it != dist.end(); ++it) {
                        const std::string& rcpt = it->first;
                        double amt = it->second;
                        transferredAmounts[rcpt] += amt;
                        recipients[rcpt] = cand; // source is the elected candidate
                    }
                    transferSurplus(cand, surplus, ballots, voteCounts, elected);
                }
            }

            // Print a CSV round showing the transfers (if any)
            printCsvRound(round++, voteCounts, elected, allCandidates, transferredAmounts, recipients);

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
            if (!hasCandidate || votes < minVotes) {
                hasCandidate = true;
                minVotes = votes;
            }
        }
        if (!hasCandidate) break; // no candidates left

        // Find all with the minimum vote
        std::set<std::string> tied;
        for (const auto& [cand, votes] : voteCounts) {
            if (elected.count(cand)) continue;
            if (votes == minVotes) tied.insert(cand);
        }

        // Tie case:
        std::string toEliminate = (tied.size() == 1)
            ? *tied.begin()
            : resolveTieSTVWithHistory(history, ballots, tied);

        // Compute elimination distribution for logging before mutating voteCounts
        auto elimDist = computeEliminationDistributionForLog(toEliminate, ballots, elected, voteCounts);

        // Transfer eliminated candidate's ballots to next continuing preferences
        transferEliminatedVotes(toEliminate, ballots, voteCounts, elected);

        // Remove eliminated candidate from the tally
        voteCounts.erase(toEliminate);

        // Prepare and print CSV round for elimination transfers
        {
            std::map<std::string, std::string> recipients;
            for (const auto& [rcpt, _] : elimDist) {
                recipients[rcpt] = toEliminate; // source is the eliminated candidate
            }
            printCsvRound(round++, voteCounts, elected, allCandidates, elimDist, recipients);
        }

        // Record the state for future §11D resolutions
        history.push_back(voteCounts);

        forceElectIfOnlyAsManyLeft();
        if (static_cast<int>(elected.size()) >= seats) break;
        if (voteCounts.empty()) break;
    }

    // Fallback: if seats remain (e.g., due to placeholder transfers), fill by top remaining tallies
    if (static_cast<int>(elected.size()) < seats) {
        std::vector<std::pair<std::string, double>> remaining;
        for (const auto& [cand, v] : voteCounts)
            if (!elected.count(cand)) remaining.emplace_back(cand, v);

        std::sort(remaining.begin(), remaining.end(), [](const auto& l, const auto& r) {
            if (l.second != r.second) return l.second > r.second;
            return l.first < r.first;
        });

        for (const auto& [cand, _] : remaining) {
            if (static_cast<int>(elected.size()) >= seats) break;
            elected.insert(cand);
        }

        // Final snapshot (no transfers in fallback)
        std::map<std::string, double> noneAmt;
        std::map<std::string, std::string> noneSrc;
        printCsvRound(round++, voteCounts, elected, allCandidates, noneAmt, noneSrc);
    }

    return elected;
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
    // Mirrors §11-compliant logic in transferSurplus:
    // Only ballots assigned to 'candidate' that have a next continuing preference are transferable.
    std::map<std::string, double> result;
    if (surplus <= 0.0) return result;

    auto isContinuing = [&](const std::string& c) -> bool {
        return elected.count(c) == 0 && voteCounts.find(c) != voteCounts.end();
    };

    int transferable = 0;
    std::map<std::string, int> nextCounts;

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

    const double transferValue = surplus / static_cast<double>(transferable);

    for (const auto& [rcpt, cnt] : nextCounts) {
        result[rcpt] = transferValue * static_cast<double>(cnt);
    }

    return result;
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

#ifdef _WIN32
    std::cout << "\nPress any key to exit...\n";
    (void)_getch();      // does not require Enter and does not echo
#else
    std::cout << "\nPress Enter to exit...\n";
    std::cin.get();      // portable fallback (requires Enter)
#endif
    return 0;
}
