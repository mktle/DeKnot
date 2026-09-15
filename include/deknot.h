#ifndef DEKNOT_H
#define DEKNOT_H

#define DK_VERSION "0.1"

#ifdef __cplusplus
extern "C" {
#endif

/***********************
 * Assembly parameters *
 ***********************/
typedef struct {
    // s value when selecting syncmers
    int s;
    // initial number of syncmers in each k-syncmer
    int init_k_prime;
    // number of k' growth iterations before the assembler will be forced to stop
    int max_iterations;
    // drop unitigs below this frequency during compaction
    int compact_thresh;
    // set to non-zero value to return homopolymer compressed unitigs
    int return_hpc_utg;
    // directory to output intermediate graph files (NULL if none)
    const char* out_dir;
    // non-zero value to print logging statements
    int verbose;
    // (bool) whether errors should be dropped based on strand evidence (e.g., ONT)
    int check_strand;
} dk_opt_t;

/*********************************
 * External graph representation *
 *********************************/
typedef struct {
    // id of source unitig (i.e., the index of the unitig sequence in the dk_graph_t char** unitigs member)
    int from_utg;
    int to_utg;
    char from_dir;
    char to_dir;
} dk_edge_t;

typedef struct {
    char** unitigs;
    dk_edge_t* edges;
    int n_unitigs;
    int n_edges;
} dk_graph_t;

typedef struct {
    // expanded (non-HPC) haplotype walk sequences
    char** walks;
    int n_walks;
    // per input read (same order as reads_arr): 1/2 = haplotype, 0 = homozygous/unassigned
    int* read_assignments;
} dk_walks_t;

/************************
 * Function definitions *
 ************************/

/**
 * Initializes and sets default assembly parameters *
 * 
 * @param opt           opt object to be initialized
 */
void dk_opt_init(dk_opt_t* opt);

/**
 * Sets assembly 
 *
 * @param preset        one of {hifi, ont} //TODO: implement
 * @param opt           opt object to be set
 * 
 * @return 0 on success, -1 for unknown preset
 */
int dk_set_opt(const char* preset, dk_opt_t* opt);

/**
 * Performs homopolymer expansion on compressed unitigs
 * 
 * @param opt           assembly parameters
 * @param reads_arr     original input read sequences
 * @param read_names    names corresponding to the reads in reads_arr (optional and can be set to NULL)
 * @param n_reads       number of reads in reads_arr
 * 
 * @return              final assembly graph represented as a list of unitigs and a list of edges
 */
dk_graph_t* dk_assemble_reads(dk_opt_t* opt, char** reads_arr, char** read_names, int n_reads);

/**
 * Assembles reads and returns the phased haplotype walks directly
 *
 * @param opt           assembly parameters
 * @param reads_arr     input read sequences
 * @param read_names    names corresponding to the reads (optional, can be NULL); must be unique when provided, since they are used to map haplotype assignments back to input order
 * @param n_reads       number of reads in reads_arr
 *
 * @return              phased walks and per-read haplotype assignments (NULL on failure); free with dk_walks_free
 */
dk_walks_t* dk_assemble_walks(dk_opt_t* opt, char** reads_arr, char** read_names, int n_reads);

/**
 * Frees a dk_walks_t returned by dk_assemble_walks
 *
 * @param to_free       walks object to be freed
 */
void dk_walks_free(dk_walks_t* to_free);

/**
 * Performs homopolymer expansion on compressed unitigs
 * 
 * @param opt               assembly parameters
 * @param reads_arr         original input read sequences (must not be homopolymer compressed)
 * @param hpc_seqs_arr      homopolymer compressed unitigs to expand
 * @param n_unitigs         number of sequences to expand (i.e., hpc_seqs size)
 * 
 * @return homopolymer expanded unitigs, where len(ret_seqs) == n_unitigs == len(hpc_seqs)
 */
char** dk_expand_sequences(dk_opt_t* opt, char** reads_arr, int n_reads, char** hpc_seqs_arr, int n_unitigs);

/**
 * Writes a given graph in .gfa format
 *
 * @param out_file      path of .gfa file to be written
 * @param graph         graph object to be written
 */
void dk_write_gfa(const char* file_path, dk_graph_t* graph);

/**
 * Frees the dynamically allocated arrays within a dk_graph_t
 * 
 * @param to_free       graph to be freed
 */
void dk_graph_free(dk_graph_t* to_free);

/*
 * Helper function for debugging/checkpointing
 */
bool check_for_checkpoint(int curr_k_prime);

#ifdef __cplusplus
}
#endif
#endif
