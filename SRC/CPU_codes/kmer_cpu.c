#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <chrono>

// Constants for scoring
const int MATCH_SCORE = 1;
const int GAP_COST = -1;
const int MISMATCH_COST = -1;
const int KMER_LENGTH = 15; // Length of k-mers

// Struct for alignment result
struct AlignmentResult {
    int score;
    std::string aligned_seq1;
    std::string aligned_seq2;
};

// Function Declarations
std::vector<std::string> extract_kmers(const std::string& seq, int k);
void build_kmer_index(const std::vector<std::string>& sequences, int k, std::unordered_map<std::string, std::vector<int>>& kmer_index);
AlignmentResult smith_waterman(const std::string& seq1, const std::string& seq2);
AlignmentResult search_and_align_with_kmers(const std::string& query,
                                            const std::unordered_map<std::string, std::vector<int>>& kmer_index,
                                            const std::vector<std::string>& reference_sequences);
std::string alignment_markup(const std::string& seq1, const std::string& seq2);
double calculate_percent_identity(const std::string& seq1, const std::string& seq2);
std::vector<std::string> load_sequences_from_file(const std::string& filename);
void process_queries_sequentially(const std::vector<std::string>& queries,
                                  const std::unordered_map<std::string, std::vector<int>>& kmer_index,
                                  const std::vector<std::string>& reference_sequences);

// Function Definitions

// Extract k-mers using a sliding window
std::vector<std::string> extract_kmers(const std::string& seq, int k) {
    std::vector<std::string> kmers;
    if (seq.size() < k) return kmers;
    for (size_t i = 0; i <= seq.size() - k; ++i) {
        kmers.push_back(seq.substr(i, k));
    }
    return kmers;
}

// Build k-mer index from reference sequences
void build_kmer_index(const std::vector<std::string>& sequences, int k, std::unordered_map<std::string, std::vector<int>>& kmer_index) {
    for (int i = 0; i < sequences.size(); ++i) {
        std::vector<std::string> kmers = extract_kmers(sequences[i], k);
        for (const auto& kmer : kmers) {
            kmer_index[kmer].push_back(i);
        }
    }
}

// Smith-Waterman alignment
AlignmentResult smith_waterman(const std::string& seq1, const std::string& seq2) {
    int m = seq1.size();
    int n = seq2.size();
    std::vector<std::vector<int>> score_matrix(m + 1, std::vector<int>(n + 1, 0));

    int max_score = 0;
    int max_i = -1, max_j = -1;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? MATCH_SCORE : MISMATCH_COST);
            int delete_ = score_matrix[i - 1][j] + GAP_COST;
            int insert = score_matrix[i][j - 1] + GAP_COST;
            score_matrix[i][j] = std::max({0, match, delete_, insert});

            if (score_matrix[i][j] > max_score) {
                max_score = score_matrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    std::string aligned_seq1 = "", aligned_seq2 = "";
    int i = max_i, j = max_j;
    while (i > 0 && j > 0 && score_matrix[i][j] > 0) {
        if (score_matrix[i][j] == score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? MATCH_SCORE : MISMATCH_COST)) {
            aligned_seq1 = seq1[i - 1] + aligned_seq1;
            aligned_seq2 = seq2[j - 1] + aligned_seq2;
            --i;
            --j;
        } else if (score_matrix[i][j] == score_matrix[i - 1][j] + GAP_COST) {
            aligned_seq1 = seq1[i - 1] + aligned_seq1;
            aligned_seq2 = "-" + aligned_seq2;
            --i;
        } else {
            aligned_seq1 = "-" + aligned_seq1;
            aligned_seq2 = seq2[j - 1] + aligned_seq2;
            --j;
        }
    }

    return {max_score, aligned_seq1, aligned_seq2};
}

// Search using k-mer index with cache and perform alignment
AlignmentResult search_and_align_with_kmers(const std::string& query,
                                            const std::unordered_map<std::string, std::vector<int>>& kmer_index,
                                            const std::vector<std::string>& reference_sequences) {
    std::unordered_map<int, AlignmentResult> results;
    std::vector<std::string> query_kmers = extract_kmers(query, KMER_LENGTH);

    for (const auto& kmer : query_kmers) {
        if (kmer_index.find(kmer) != kmer_index.end()) {
            for (int idx : kmer_index.at(kmer)) {
                AlignmentResult result = smith_waterman(query, reference_sequences[idx]);
                if (results.find(idx) == results.end() || result.score > results[idx].score) {
                    results[idx] = result;
                }
            }
        }
    }

    // Find the best alignment
    int best_score = 0;
    AlignmentResult best_result;
    for (const auto& entry : results) {
        if (entry.second.score > best_score) {
            best_score = entry.second.score;
            best_result = entry.second;
        }
    }

    return best_result;
}

// Generate alignment markup line
std::string alignment_markup(const std::string& seq1, const std::string& seq2) {
    std::string markup;
    for (size_t i = 0; i < seq1.length(); ++i) {
        if (i < seq2.length() && seq1[i] == seq2[i] && seq1[i] != '-') {
            markup += "|"; // Add a match indicator for aligned bases
        } else {
            markup += " "; // Add a space for mismatches or gaps
        }
    }
    return markup;
}

// Calculate percent identity
double calculate_percent_identity(const std::string& seq1, const std::string& seq2) {
    int matches = 0;
    int aligned_length = 0;

    for (size_t i = 0; i < seq1.length() && i < seq2.length(); ++i) {
        if (seq1[i] == '-' || seq2[i] == '-') {
            aligned_length++;
        } else {
            aligned_length++;
            if (seq1[i] == seq2[i]) {
                matches++;
            }
        }
    }

    return (aligned_length > 0) ? (static_cast<double>(matches) / aligned_length) * 100.0 : 0.0;
}

// Load sequences from a file
std::vector<std::string> load_sequences_from_file(const std::string& filename) {
    std::vector<std::string> sequences;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            sequences.push_back(line);
        }
    }
    file.close();
    return sequences;
}

// Process queries sequentially
void process_queries_sequentially(const std::vector<std::string>& queries,
                                  const std::unordered_map<std::string, std::vector<int>>& kmer_index,
                                  const std::vector<std::string>& reference_sequences) {
    for (const auto& query : queries) {
        AlignmentResult result = search_and_align_with_kmers(query, kmer_index, reference_sequences);

        // Print alignment result
        std::cout << "Query Sequence: " << query << "\n";
        std::cout << "Smith-Waterman alignment score: " << result.score << "\n";
        std::cout << "Aligned Query Sequence:      " << result.aligned_seq1 << "\n";
        std::cout << "                             " << alignment_markup(result.aligned_seq1, result.aligned_seq2) << "\n";
        std::cout << "Aligned Reference Sequence:  " << result.aligned_seq2 << "\n";
        double percent_identity = calculate_percent_identity(result.aligned_seq1, result.aligned_seq2);
        std::cout << "Percent Identity: " << percent_identity << "%\n";
        std::cout << "Alignment Length: " << result.aligned_seq1.length() << "\n\n";
    }
}

// Main function
int main() {
    std::vector<std::string> reference_sequences = load_sequences_from_file("../Database/database.txt");
    std::vector<std::string> queries = load_sequences_from_file("../Database/query.txt");

    if (reference_sequences.empty() || queries.empty()) {
        std::cerr << "Error: No sequences loaded. Please check your input files.\n";
        return 1;
    }

    // Build k-mer index
    std::unordered_map<std::string, std::vector<int>> kmer_index;
    build_kmer_index(reference_sequences, KMER_LENGTH, kmer_index);

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Process queries sequentially
    process_queries_sequentially(queries, kmer_index, reference_sequences);

    // Stop timing
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;

    // Calculate and print throughput
    double throughput = queries.size() / duration.count();
    std::cout << "Total Execution Time: " << duration.count() << " seconds\n";
    std::cout << "Throughput: " << throughput << " queries/second\n";

    return 0;
}
