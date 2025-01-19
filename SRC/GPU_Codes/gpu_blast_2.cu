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
const int N_MAX=150;
const int TILE_SIZE = 150;
// Constants for Alignment tracebacking
// const int NONE = 0;
// const int DIAGONAL = 1;
// const int UP = 2;
// const int LEFT = 3;
__constant__ char* ConstQuery;

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
            int match = seq1[i - 1] == ref_seq[j - 1] ? MATCH_SCORE : MISMATCH_COST;

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
    int m = h_seq1.size();
    int num_refs = h_seq2_list.size();

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
    cudaDeviceSynchronize();

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
    cudaFree(d_alignment_lengths);

    return results;
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

    char key = query[0];

    AlignmentResult final_result = {0,"",""};

    if (sequence_index.find(key) != sequence_index.end()) {

        std::vector<std::string>& ref_strings_matched = sequence_index[key];


        std::vector<AlignmentResult> results = smith_waterman_batch(query, ref_strings_matched);

        int max_score = -1;

        for (const auto& result : results) {
            if (result.score > max_score) {
                max_score = result.score;
                final_result = result;
            }
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
    std::vector<std::string> reference_sequences = load_sequences_from_file("proj_ref.txt");
    // Load query sequences from file
    std::vector<std::string> queries = load_sequences_from_file("proj_query.txt");

    // Build the hierarchical index
    build_index(reference_sequences);

    // Choose alignment algorithm (Needleman-Wunsch or Smith-Waterman)
    AlignmentAlgorithm algorithm = SMITH_WATERMAN;

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Process the queries in adaptive batches
    process_batches(queries, algorithm);

    // Stop timing
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;

    // Calculate and print throughput
    double throughput = queries.size() / duration.count();
    std::cout << "Total Execution Time: " << duration.count() << " seconds\n";
    std::cout << "Throughput: " << throughput << " queries/second\n";
    std::cout << "Queries BLASTED: " << queries.size() << "\n";

    return 0;
}


