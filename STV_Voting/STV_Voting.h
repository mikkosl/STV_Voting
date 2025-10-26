// STV_Voting.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>
#include <string>
#include <vector>

// Read candidate names from stdin, one per line, until an empty line.
// - Trims surrounding whitespace.
// - Ignores duplicate names (case-sensitive).
// - Preserves insertion order for unique names.
std::vector<std::string> inputCandidateNames();

// Read number of seats (positive integer) from stdin with basic validation.
// Re-prompts until a valid value is entered. Returns the parsed value.
int inputNumberOfSeats();

// Display a row-numbered candidate list (1-based).
void displayCandidateList(const std::vector<std::string>& candidates);

// Input multiple ballots, one per line as comma-separated row numbers.
// - Empty line ends input.
// - Displays candidates once and prompts for ballots.
// - Skips empty/invalid ballots.
std::vector<std::vector<std::string>> inputBallotsByRowNumbers(const std::vector<std::string>& candidates);
