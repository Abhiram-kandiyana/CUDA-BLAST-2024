#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <thread>


// Constants for scoring
const int MATCH_SCORE = 1;
const int GAP_COST = -1;
const int MISMATCH_COST = -1;
// Constants for Alignment tracebacking
const int NONE = 0;
const int DIAGONAL = 1;
const int UP = 2;
const int LEFT = 3;


// Enum for choosing alignment algorithm
enum AlignmentAlgorithm { NEEDLEMAN_WUNSCH, SMITH_WATERMAN };
struct AlignmentResult {
   int score;
   std::string aligned_seq1;
   std::string aligned_seq2;
};


// Function to perform Needleman-Wunsch alignment
int needleman_wunsch(const std::string& seq1, const std::string& seq2) {
   int m = seq1.size();
   int n = seq2.size();
   std::vector<std::vector<int>> score_matrix(m + 1, std::vector<int>(n + 1, 0));


   // Initialize score matrix with gap penalties
   for (int i = 0; i <= m; ++i) score_matrix[i][0] = i * GAP_COST;
   for (int j = 0; j <= n; ++j) score_matrix[0][j] = j * GAP_COST;


   // Fill score matrix
   for (int i = 1; i <= m; ++i) {
       for (int j = 1; j <= n; ++j) {
           int match = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? MATCH_SCORE : MISMATCH_COST);
           int delete_ = score_matrix[i - 1][j] + GAP_COST;
           int insert = score_matrix[i][j - 1] + GAP_COST;
       }
   }
   return score_matrix[m][n];  // Return global alignment score
}


// Function to perform Smith-Waterman alignment
AlignmentResult smith_waterman(const std::string& seq1, const std::string& seq2) {
   int m = seq1.size();
   int n = seq2.size();
   std::vector<std::vector<int>> score_matrix(m + 1, std::vector<int>(n + 1, 0));
   std::vector<std::vector<int>> traceback_matrix(m + 1, std::vector<int>(n + 1, -1));


   int max_score = 0;  // Keep track of the maximum score for local alignment
   int max_i=-1;
   int max_j=-1;




   // Fill score matrix
   for (int i = 1; i <= m; ++i) {
       for (int j = 1; j <= n; ++j) {
           int match = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? MATCH_SCORE : MISMATCH_COST);
           int delete_ = score_matrix[i - 1][j] + GAP_COST;
           int insert = score_matrix[i][j - 1] + GAP_COST;
           int direction = NONE;
           score_matrix[i][j] = std::max({0, match, delete_, insert});  // Reset negative scores to zero
           if (score_matrix[i][j] > max_score) {
               max_score = score_matrix[i][j]; // Update max score
               max_i = i; //Update max i
               max_j = j; //Update max j
           }
           //Fill Traceback matrix
           if (score_matrix[i][j] == match) direction = DIAGONAL;
           else if (score_matrix[i][j] == delete_) direction = UP;
           else if (score_matrix[i][j] == insert) direction = LEFT;
           traceback_matrix[i][j] = direction;


       }
   }
   // Traceback to get the alignment seq
   int i = max_i;
   int j = max_j;
   std::string aligned_seq1 = "";
   std::string aligned_seq2 = "";
   while (score_matrix[i][j] > 0) {
       if (traceback_matrix[i][j] == DIAGONAL) {
         aligned_seq1 = seq1[i - 1] + aligned_seq1;
         aligned_seq2 = seq2[j - 1] + aligned_seq2;
         i--;
         j--;
       }
       else if (traceback_matrix[i][j] == UP) {
         aligned_seq1 = seq1[i - 1] + aligned_seq1;
         aligned_seq2 = "-" + aligned_seq2;
         i--;
       }
       else if (traceback_matrix[i][j] == LEFT) {
         aligned_seq1 = "-" + aligned_seq1;
         aligned_seq2 = seq2[j - 1] + aligned_seq2;
         j--;
       }
       else{
         printf("Error while running traceback");
         break;
       }
   }
   return {max_score,aligned_seq1,aligned_seq2};  // Return highest local alignment score
}


// Hierarchical Index and Cache setup
std::unordered_map<char, std::vector<std::string>> sequence_index;
std::unordered_map<std::string, AlignmentResult> cache;


// Build a simple index based on the first character of sequences
void build_index(const std::vector<std::string>& sequences) {
   for (const auto& seq : sequences) {
       sequence_index[seq[0]].push_back(seq);
   }
}


// Search in index and use cache for faster access
AlignmentResult search_and_align(const std::string& query, AlignmentAlgorithm algorithm) {
   // Check cache
   if (cache.find(query) != cache.end()) {
       return cache[query];
   }


   // Search in index based on the first character
   char key = query[0];
   int max_score = 0;
   std::string final_aligned_seq1 ="";
   std::string final_aligned_seq2 ="";
   AlignmentResult result;
   AlignmentResult final_result;


   if (sequence_index.find(key) != sequence_index.end()) {
       for (const auto& seq : sequence_index[key]) {
           int score;
           if (algorithm == NEEDLEMAN_WUNSCH) {
               score = needleman_wunsch(query, seq);
           } else {  // SMITH_WATERMAN
               result = smith_waterman(query, seq);
           }
           if (result.score > max_score) {
               max_score = result.score;
               final_result = result;


           }
           // max_score = std::max(max_score, score);
       }
   }


   // Store in cache
   cache[query] = final_result;


   return final_result;
   // return {max_score,final_aligned_seq1,final_aligned_seq2};
}


std::string alignment_markup(const std::string& seq1, const std::string& seq2) {
   std::string markup;
   for (size_t i = 0; i < seq1.length(); ++i) {
       if (seq1[i] == seq2[i]) {
           if (seq1[i] != '-') {  // Exclude gaps
               markup += "|";
           } else {
               markup += " ";
           }
       } else {
           markup += " ";
       }
   }
   return markup;
}


double calculate_percent_identity(const std::string& seq1, const std::string& seq2) {
   int matches = 0;
   int aligned_length = 0;
   for (size_t i = 0; i < seq1.length(); ++i) {
       if (seq1[i] != '-' && seq2[i] != '-') {
           aligned_length++;
           if (seq1[i] == seq2[i]) {
               matches++;
           }
       }
   }
   return (aligned_length > 0) ? (static_cast<double>(matches) / aligned_length) * 100.0 : 0.0;
}




// Function to simulate adaptive batching (for simplicity, it processes a batch size of 5)
void process_batches(const std::vector<std::string>& queries, AlignmentAlgorithm algorithm) {
   const int batch_size = 5;
   for (size_t i = 0; i < queries.size(); i += batch_size) {
       std::cout << "Processing batch " << (i / batch_size + 1) << "...\n";
       for (size_t j = i; j < std::min(i + batch_size, queries.size()); ++j) {
           AlignmentResult final_result = search_and_align(queries[j], algorithm);
           std::cout << "Query Sequence: " << queries[j] << "\n";
           std::cout << "Smith-Waterman alignment score: " << final_result.score << "\n";
           std::cout << "Aligned Query Sequence:      " << final_result.aligned_seq1 << "\n";
           std::cout << "                             " << alignment_markup(final_result.aligned_seq1, final_result.aligned_seq2) << "\n";
           std::cout << "Aligned Reference Sequence:  " << final_result.aligned_seq2 << "\n";
           double percent_identity = calculate_percent_identity(final_result.aligned_seq1, final_result.aligned_seq2);
           std::cout << "Percent Identity: " << percent_identity << "%\n";
           std::cout << "Alignment Length: " << final_result.aligned_seq1.length() << "\n\n";
       }
       std::this_thread::sleep_for(std::chrono::milliseconds(500));  // Simulate processing delay
   }
}


// Function to load sequences from a file
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


// Main function
int main() {
   // Load reference sequences from file
   std::vector<std::string> reference_sequences = load_sequences_from_file("../Database/database.txt");
   // Load query sequences from file
   std::vector<std::string> queries = load_sequences_from_file("../Database/query.txt");


   // Build the hierarchical index
   build_index(reference_sequences);


   // Choose alignment algorithm (Needleman-Wunsch or Smith-Waterman)
   AlignmentAlgorithm algorithm = SMITH_WATERMAN;


   // Process the queries in adaptive batches
   process_batches(queries, algorithm);


   return 0;
}
