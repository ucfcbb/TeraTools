#ifndef R_SA_LCP_FASTA_H
#define R_SA_LCP_FASTA_H

#include"../thirdparty/kseq.h"
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#include<chrono>
#include<omp.h>
#include<atomic>
#include<mutex>
#include<thread>
#include<condition_variable>

KSEQ_INIT(int, read)

kseq_t* open_fasta(const std::string& fasta_file, FILE** fp) {
    *fp = fopen(fasta_file.c_str(), "r");
    if (!*fp) {
        if (errno == ENOENT) {
            std::cerr << "File not found (" + fasta_file + ")" << std::endl;
            exit(1);
        } else {
            std::cerr << "Problem opening file (" + fasta_file + "), " + strerror(errno) << std::endl;
            exit(1);
        }
        return nullptr;
    }
    return kseq_init(fileno(*fp));
}

// Stores sequence information for multi-threaded processing
struct SeqInfo {
    char* seq_content;
    char* seq_name;
    int64_t seq_len;
    int64_t seq_name_len;
    int64_t seq_comment_len;
};

template<typename ProcessFunc>
void process_sequences(kseq_t* seq, uint16_t threads, ProcessFunc process_func) {
    #pragma omp parallel
    {
        SeqInfo seq_info;
        while (true) {
            #pragma omp critical(read_seq)
            {
                seq_info.seq_len = kseq_read(seq);
                if (seq_info.seq_len >= 0) {
                    if (threads > 1) {
                        seq_info.seq_content = strdup(seq->seq.s);
                        seq_info.seq_name = strdup(seq->name.s);  // Copy name for multi-threaded
                    } 
                    // avoid extra copy for single thread
                    else {
                        seq_info.seq_content = seq->seq.s;
                        seq_info.seq_name = seq->name.s;
                    }
                    seq_info.seq_name_len = seq->name.l;
                    seq_info.seq_comment_len = seq->comment.l;
                }
            }
            if (seq_info.seq_len <= 0) {
                break;
            }

            process_func(seq_info);
        }

        if (threads > 1) {
            free(seq_info.seq_content);
            free(seq_info.seq_name);
        }
    }
}
#endif //#ifndef R_SA_LCP_FASTA_H