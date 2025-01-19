#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <cuda.h>

// Constants for scoring
const int MATCH_SCORE = 1;
const int GAP_COST = -1;
const int MISMATCH_COST = -1;
const int KMER_LENGTH = 15; // Length of k-mers
const int N_MAX=150;

// Global variables to accumulate timings
static double total_gpu_time = 0.0;        // Accumulate total GPU kernel time
static double total_kmer_search_time = 0.0; // Accumulate total time picking candidate sequences via k-mers
static double total_calculate_percent_identity_time = 0.0;


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
    if (seq.size() < (size_t)k) return kmers;
    for (size_t i = 0; i <= seq.size() - k; ++i) {
        kmers.push_back(seq.substr(i, k));
    }
    return kmers;
}

// Build k-mer index from reference sequences
// Build k-mer index from reference sequences
void build_kmer_index(const std::vector<std::string>& sequences, int k, std::unordered_map<std::string, std::vector<int>>& kmer_index) {
    for (int i = 0; i < (int)sequences.size(); ++i) {
        std::vector<std::string> kmers = extract_kmers(sequences[i], k);
        for (const auto& kmer : kmers) {
            kmer_index[kmer].push_back(i);
        }
    }
}

// Kernel Function
__global__ void smith_waterman_kernel_batch(
    const char* seq1,
    const char* seq2,
    int m,
    int num_refs,
    int* score_matrices,
    int* max_scores,
    int* max_positions,
    char* aligned_seq1_array,
    char* aligned_seq2_array,
    int* alignment_lengths
)
{
    int n = N_MAX;

    int ref_string_idx = blockIdx.x;
    if(ref_string_idx >= num_refs) 
        return;

    __shared__ char s_ref_seq[N_MAX];

    int tid = threadIdx.x;

    const char* ref_seq = seq2 + ref_string_idx * N_MAX;
    int* score_matrix = score_matrices + ref_string_idx * (m+1) * (N_MAX + 1);
    int* max_score = max_scores + ref_string_idx;
    int* max_pos = max_positions + ref_string_idx * 2;

    // Pointers for alignment results
    char* aligned_seq1 = aligned_seq1_array + ref_string_idx * (m + N_MAX);  // Max possible alignment length
    char* aligned_seq2 = aligned_seq2_array + ref_string_idx * (m + N_MAX);
    int* alignment_length = alignment_lengths + ref_string_idx;

    // Load the reference sequence into shared memory
    for (int i = tid; i < N_MAX; i += blockDim.x) {
        s_ref_seq[i] = ref_seq[i];
    }
    __syncthreads();  // Ensure all threads have loaded the reference sequence

    for (int k = 2; k <= m + N_MAX; ++k) {
        int start_i = max(1, k - N_MAX);
        int end_i = min(m, k - 1);
        int num_elements = end_i - start_i + 1;

        for (int idx = tid; idx < num_elements; idx += blockDim.x) {
            int i = start_i + idx;
            int j = k - i;

            // Compute match/mismatch score
            int match = (seq1[i - 1] == ref_seq[j - 1]) ? MATCH_SCORE : MISMATCH_COST;

            // Calculate scores from neighboring cells
            int diag_score = score_matrix[(i - 1) * (N_MAX + 1) + (j - 1)] + match;
            int up_score = score_matrix[(i - 1) * (N_MAX + 1) + j] + GAP_COST;
            int left_score = score_matrix[i * (N_MAX + 1) + (j - 1)] + GAP_COST;

            int cell_score = max(0, max(diag_score, max(up_score, left_score)));

            score_matrix[i * (N_MAX + 1) + j] = cell_score;

             // Atomically update max score and position
            int old_max = atomicMax(max_score, cell_score);
            if (cell_score > old_max) {
                    atomicExch(&max_pos[0], i);
                    atomicExch(&max_pos[1], j);
                }
            }
            __syncthreads();
        }

        if(tid == 0)
        {
            int i = max_pos[0];
            int j = max_pos[1];
            int pos = 0;

            // Perform traceback until score is zero
            while (i > 0 && j > 0 && score_matrix[i * (n + 1) + j] > 0) {
                int current_score = score_matrix[i * (n + 1) + j];
                int diag_score = score_matrix[(i - 1) * (n + 1) + (j - 1)];
                int up_score = score_matrix[(i - 1) * (n + 1) + j];
                int left_score = score_matrix[i * (n + 1) + (j - 1)];

                if (current_score == diag_score + (seq1[i - 1] == s_ref_seq[j - 1] ? MATCH_SCORE : MISMATCH_COST)) {
                    aligned_seq1[pos] = seq1[i - 1];
                    aligned_seq2[pos] = s_ref_seq[j - 1];
                    --i;
                    --j;
                } else if (current_score == up_score + GAP_COST) {
                    aligned_seq1[pos] = seq1[i - 1];
                    aligned_seq2[pos] = '-';
                    --i;
                } else if (current_score == left_score + GAP_COST) {
                    aligned_seq1[pos] = '-';
                    aligned_seq2[pos] = s_ref_seq[j - 1];
                    --j;
                } else {
                    printf("Should not reach here in Smith-Waterman. Something is wrong");
                    break;
                }
                ++pos;
            }
            alignment_length[0] = pos;
            }
}

// Smith-Waterman alignment
std::vector<AlignmentResult> smith_waterman_batch(const std::string& h_seq1, const std::vector<std::string>& h_seq2_list) {
    int m = (int)h_seq1.size();
    int num_refs = (int)h_seq2_list.size();

    char* d_seq1;
    cudaMalloc(&d_seq1, m * sizeof(char));
    cudaMemcpy(d_seq1, h_seq1.data(), m * sizeof(char), cudaMemcpyHostToDevice);

    char* d_seq2;
    int total_ref_len = N_MAX * num_refs;

    cudaMalloc(&d_seq2, total_ref_len * sizeof(char));

    for(int i=0;i < num_refs;i++)
    {
        cudaMemcpy(d_seq2+i*N_MAX,h_seq2_list[i].data(),N_MAX*sizeof(char),cudaMemcpyHostToDevice); 
    }

    int* d_score_matrices;
    int score_matrix_size = (m + 1) * (N_MAX + 1);
    cudaMalloc(&d_score_matrices, num_refs * score_matrix_size * sizeof(int));

    int* d_max_scores;
    int* d_max_positions;
    cudaMalloc(&d_max_scores, num_refs * sizeof(int));
    cudaMalloc(&d_max_positions, num_refs * 2 * sizeof(int));


    // Initialize score matrices and max scores/positions
    cudaMemset(d_score_matrices, 0, num_refs * score_matrix_size * sizeof(int));
    cudaMemset(d_max_scores, 0, num_refs * sizeof(int));
    cudaMemset(d_max_positions, 0, num_refs * 2 * sizeof(int));

    //Device memory for alignment lengths
    char* d_aligned_seq1_arr;
    char* d_aligned_seq2_arr;
    int* d_alignment_lengths;
    int max_alignment_length = m + N_MAX;

    cudaMalloc(&d_aligned_seq1_arr, num_refs * max_alignment_length * sizeof(char));
    cudaMalloc(&d_aligned_seq2_arr, num_refs * max_alignment_length * sizeof(char));
    cudaMalloc(&d_alignment_lengths, num_refs * sizeof(int));

    cudaMemset(d_alignment_lengths, 0, num_refs * sizeof(int));

    int threadsPerBlock = N_MAX;
    int blocksPerGrid = num_refs;

    // Measure GPU kernel time using CUDA events
    cudaEvent_t start_ev, stop_ev;
    cudaEventCreate(&start_ev);
    cudaEventCreate(&stop_ev);

    cudaEventRecord(start_ev);

    smith_waterman_kernel_batch<<<blocksPerGrid, threadsPerBlock>>>(
        d_seq1,
        d_seq2,
        m,
        num_refs,
        d_score_matrices,
        d_max_scores,
        d_max_positions,
        d_aligned_seq1_arr,
        d_aligned_seq2_arr,
        d_alignment_lengths
    );

    cudaEventRecord(stop_ev);
    cudaEventSynchronize(stop_ev);

    float gpu_ms = 0.0f;
    cudaEventElapsedTime(&gpu_ms, start_ev, stop_ev);
    total_gpu_time += (gpu_ms / 1000.0);

    cudaEventDestroy(start_ev);
    cudaEventDestroy(stop_ev);
    // cudaDeviceSynchronize();

    std::vector<int> h_max_scores(num_refs);
    std::vector<int> h_alignment_lengths(num_refs);
    std::vector<AlignmentResult> results(num_refs);

    cudaMemcpy(h_max_scores.data(), d_max_scores, num_refs * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_alignment_lengths.data(), d_alignment_lengths, num_refs * sizeof(int), cudaMemcpyDeviceToHost);

    for (int i = 0; i < num_refs; ++i) {
        int alignment_length = h_alignment_lengths[i];
        char* h_aligned_seq1 = (char*)malloc(alignment_length * sizeof(char));
        char* h_aligned_seq2 = (char*)malloc(alignment_length * sizeof(char));
    
        cudaMemcpy(h_aligned_seq1, d_aligned_seq1_arr + i * max_alignment_length, alignment_length * sizeof(char), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_aligned_seq2, d_aligned_seq2_arr + i * max_alignment_length, alignment_length * sizeof(char), cudaMemcpyDeviceToHost);
    
        // Build the aligned sequences (reverse since they were constructed backwards)
        std::string aligned_seq1(h_aligned_seq1, alignment_length);
        std::string aligned_seq2(h_aligned_seq2, alignment_length);
        std::reverse(aligned_seq1.begin(), aligned_seq1.end());
        std::reverse(aligned_seq2.begin(), aligned_seq2.end());
    
        // Store the result
        results[i] = {h_max_scores[i], aligned_seq1, aligned_seq2};
    
        free(h_aligned_seq1);
        free(h_aligned_seq2);
    }

    cudaFree(d_seq1);
    cudaFree(d_seq2);
    cudaFree(d_score_matrices);
    cudaFree(d_max_scores);
    cudaFree(d_max_positions);
    cudaFree(d_aligned_seq1_arr);
    cudaFree(d_aligned_seq2_arr);
    cudaFree(d_alignment_lengths);

    return results;
}

// Search using k-mer index with cache and perform alignment
AlignmentResult search_and_align_with_kmers(const std::string& query,
                                            const std::unordered_map<std::string, std::vector<int>>& kmer_index,
                                            const std::vector<std::string>& reference_sequences) {
    
    // Start timing Candidate Selection
    auto kmer_search_start = std::chrono::high_resolution_clock::now();

    std::vector<std::string> query_kmers = extract_kmers(query, KMER_LENGTH);

    std::vector<int> candidate_indices;

    for (const auto& kmer : query_kmers) {
        auto it = kmer_index.find(kmer);
        if (it != kmer_index.end()) {
            const auto& refs = it->second;
            candidate_indices.insert(candidate_indices.end(), refs.begin(), refs.end());     
        }
    }

    // Removing duplicate indices
    std::sort(candidate_indices.begin(), candidate_indices.end());
    candidate_indices.erase(std::unique(candidate_indices.begin(), candidate_indices.end()), candidate_indices.end());

    // Candidate picking done
    auto kmer_search_end = std::chrono::high_resolution_clock::now();
    total_kmer_search_time += std::chrono::duration<double>(kmer_search_end - kmer_search_start).count();

    // If no candidates are found, return an empty AlignmentResult
    if (candidate_indices.empty()) {
        return {0, "", ""};
    }

    // Extract the candidate reference sequences from the full database
    std::vector<std::string> candidate_refs;
    candidate_refs.reserve(candidate_indices.size());
    for (int idx : candidate_indices) {
        candidate_refs.push_back(reference_sequences[idx]);
    }

    std::vector<AlignmentResult> results = smith_waterman_batch(query, candidate_refs);

    // Find the best alignment
    int best_score = 0;
    AlignmentResult best_result;
    for (const auto& res : results) {
        if (res.score > best_score) {
            best_score = res.score;
            best_result = res;
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
        auto start_time = std::chrono::high_resolution_clock::now();
        double percent_identity = calculate_percent_identity(result.aligned_seq1, result.aligned_seq2);
        auto end_time = std::chrono::high_resolution_clock::now();
        total_calculate_percent_identity_time += std::chrono::duration<double>(end_time - start_time).count();

        std::cout << "Percent Identity: " << percent_identity << "%\n";
        std::cout << "Alignment Length: " << result.aligned_seq1.length() << "\n\n";
    }
}

// Main function
int main() {
    std::vector<std::string> reference_sequences = load_sequences_from_file("proj_ref.txt");
    std::vector<std::string> queries = load_sequences_from_file("proj_query.txt");

    if (reference_sequences.empty() || queries.empty()) {
        std::cerr << "Error: No sequences loaded. Please check your input files.\n";
        return 1;
    }

    // Time building the k-mer index
    auto kmer_start_time = std::chrono::high_resolution_clock::now();
    // Build k-mer index
    std::unordered_map<std::string, std::vector<int>> kmer_index;
    build_kmer_index(reference_sequences, KMER_LENGTH, kmer_index);
    auto kmer_end_time = std::chrono::high_resolution_clock::now();
    double kmer_build_time = std::chrono::duration<double>(kmer_end_time - kmer_start_time).count();
    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Process queries sequentially
    process_queries_sequentially(queries, kmer_index, reference_sequences);

    // Stop timing
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;

    // Calculate and print throughput
    double throughput = queries.size() / duration.count();
    std::cout << "Total Execution Time (All Queries): " << duration.count() << " seconds\n";
    std::cout << "Throughput: " << throughput << " queries/second\n";
    std::cout << "Queries BLASTED: " << queries.size() << "\n";

    // Print the measured times
    std::cout << "Time to build k-mer index: " << kmer_build_time << " s\n";
    std::cout << "Total GPU kernel time: " << total_gpu_time << " s\n";
    std::cout << "Total time picking candidate references via k-mers: " << total_kmer_search_time << " s\n";
    std::cout << "Total time calculating percent identity: " << total_calculate_percent_identity_time << " s\n";



    return 0;
}
