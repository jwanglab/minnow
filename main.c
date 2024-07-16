#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "incl/minimap2/minimap.h"
#include "incl/klib/kseq.h"
#include "incl/klib/kvec.h"
#include "incl/klib/khash.h"
#include "paf.h"
#include <getopt.h>
#include <string.h>

#ifndef M_E
#define M_E 2.718281828459045
#endif

KSEQ_INIT(gzFile, gzread);

// collated expression data
typedef struct expr {
  double ct;
} expr;

KHASH_MAP_INIT_STR(transcript2expr, expr);

void version() {
  printf("minnow version 0.1\n");
}

void usage() {
  printf("Usage: minnow [options]\n");
  printf("Options:\n");
  printf("  -q: FASTA/Q[.gz] file with reads\n");
  printf("  -r: Reference FASTA/Q[.gz]\n");
  printf("  -p: PAF (instead of -q and -r)\n");
  printf("  -t: Threads (default: 1)\n");
  printf("  -c: 'careful': more accurate but slower\n");
  printf("  -v, --verbose: verbose\n");
  printf("  -h, --help: show this\n");
  printf("  --version: show version information\n");
}

static struct option long_options[] = {
// if these are the same as a single-character option, put that character in the 4th field instead of 0
  { "paf",     required_argument, 0, 0 },
  { "verbose", no_argument,       0, 'v' },
  { "help",    no_argument,       0, 'h' },
  { "version", no_argument,       0, 0 },
  { "careful", no_argument,       0, 'c' },
  { 0, 0, 0, 0}
};

int main(int argc, char *argv[]) {
  mm_idxopt_t iopt;
  mm_mapopt_t mopt;

  char* read_fasta = NULL;
  char* ref_fasta = NULL;
  char* paf_file = NULL;
  char* preset = "map-ont";
  int n_threads = 1;
  int verbose = 0;
  int careful = 0;

  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "q:r:p:t:cvh", long_options, &long_idx)) != -1) {
    switch (opt) {
      case 'q':
        read_fasta = optarg;
        break;
      case 'r':
        ref_fasta = optarg;
        break;
      case 'p':
        paf_file = optarg;
        break;
      case 't':
        n_threads = atoi(optarg);
        break;
      case 'c':
        careful = 1;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        usage();
        return 0;
        break;
      case '?':
        if (optopt == 'p' || optopt == 'q' || optopt == 'r' || optopt == 't')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      case 0:
        // as long as all the long arguments have characters too, I don't think this section will be used
        if (long_idx == 0) verbose = 1; // --verbose
        else if (long_idx == 1) {usage(); return 0;} // --help
        else if (long_idx == 2) {version(); return 0;} // --version
        else if (long_idx == 3) careful = 1; // --careful
        break;
      default:
        usage();
        return 1;
    }
  }

  if(paf_file == NULL && (read_fasta == NULL || ref_fasta == NULL)) {
    fprintf(stderr, "Either -p PAF OR (-q reads AND -r reference FASTA/Q[.gz]) is required\n");
    usage();
    return 1;
  }

  mm_verbose = 2; // disable message output to stderr
  mm_set_opt(0, &iopt, &mopt); // initialize with defaults
  mm_set_opt(preset, &iopt, &mopt); // then add ont presets
  if(careful)
    mopt.flag |= MM_F_CIGAR; // perform alignment
  mopt.flag |= MM_F_NO_PRINT_2ND; // skip all secondary alignments

  // set optimized presets for sensitive FFPE alignment
  // -k12 -w1 -n2 -m20
  iopt.k = 12;
  iopt.w = 1;
  mopt.min_cnt = 2;
  mopt.min_chain_score = 20;

  // open query file for reading
  gzFile f;
  kseq_t *ks;

  // generic hashing variables
  khint_t bin; // hash bin (result of kh_put/get)
  int absent;

  // transcript accounting stuff
  khash_t(transcript2expr) *t2e = kh_init(transcript2expr);

  int aligned = 0;

  if(paf_file != NULL) {
    fprintf(stderr, "Reading from alignment file '%s'\n", paf_file);
    paf_file_t *p = paf_open(paf_file);
    if(!p) {
      fprintf(stderr, "Cannot open '%s'\n", paf_file);
      return 1;
    }
    paf_rec_t r;
    int ret = paf_read(p, &r);
    int best_matches = -1;
    char *q = NULL;
    kvec_t(khint_t) best_targets;
    kv_init(best_targets);
    while(ret == 0) {
      if(!q || strcmp(q, r.qn) != 0) {
        if(q) {
          for (int j = 0; j < kv_size(best_targets); ++j) {
            kh_val(t2e, kv_A(best_targets, j)).ct = kh_val(t2e, kv_A(best_targets, j)).ct + (1.0/kv_size(best_targets));
          }
          free(q);
        }
        q = malloc(strlen(r.qn)+1);
        q[strlen(r.qn)] = '\0';
        strcpy(q, r.qn);
        best_matches = -1;
        best_targets.n = 0; // soft reset/init
        aligned++;
      }
      if(best_matches == -1) best_matches = r.ml;
      if(r.ml == best_matches) {
        bin = kh_put(transcript2expr, t2e, r.tn, &absent);
        if(absent) {
          char *t = malloc(strlen(r.tn)+1);
          t[strlen(r.tn)] = '\0';
          strcpy(t, r.tn);
          kh_key(t2e, bin) = t;
          kh_val(t2e, bin).ct = 0;
        }
        kv_push(khint_t, best_targets, bin);
      }
      ret = paf_read(p, &r);
    }
    // last alignment
    if(q) {
      for (int j = 0; j < kv_size(best_targets); ++j) {
        kh_val(t2e, kv_A(best_targets, j)).ct = kh_val(t2e, kv_A(best_targets, j)).ct + (1.0/kv_size(best_targets));
      }
      free(q);
    }
    // clean up PAF reading
    kv_destroy(best_targets);
    paf_close(p);
    fprintf(stderr, "Closing file '%s'\n", paf_file);
  } else {
    // open index reader
    fprintf(stderr, "Building mm2 index...\n");
    mm_idx_reader_t *r = mm_idx_reader_open(ref_fasta, &iopt, 0);
    if(!r) {
      fprintf(stderr, "Cannot open '%s'\n", ref_fasta);
      return 1;
    }
    mm_idx_t *mi;
    while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
      // open (or re-open) the query file -- needs to be re-read through for each part of the index
      f = gzopen(read_fasta, "r");
      if(!f) {
        fprintf(stderr, "Cannot open '%s'\n", read_fasta);
        return 1;
      }
      ks = kseq_init(f); 

      fprintf(stderr, "Processing mm2 index...\n");
      mm_mapopt_update(&mopt, mi);
      mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
      int n = 0;

      while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
        mm_reg1_t *reg;
        int j, i, n_reg;
        if(verbose) {
          fprintf(stderr, "Processing read %d (%s, %u bp)\n", n, ks->name.s, ks->seq.l);
        }
        reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
        if(verbose) {
          fprintf(stderr, "  %d raw alignments\n", n_reg);
        }
        if(n_reg > 0) {
          aligned++;
        }

        int best_matches = -1;
        int n_best_matches = 0;
        for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
          mm_reg1_t *r = &reg[j];
          if(verbose) {
            fprintf(stderr, "    %s: %d - %d\n", mi->seq[r->rid].name, r->rs, r->re);
          }
          if(best_matches == -1) best_matches = r->mlen;
          if(r->mlen < best_matches) { // stop when we see any sub-best alignments
            free(r->p);
            break;
          }
          n_best_matches++;
          free(r->p);
        }

        // add partial ct to each based on # best matches
        for (j = 0; j < n_best_matches; ++j) {
          mm_reg1_t *r = &reg[j];
          bin = kh_put(transcript2expr, t2e, mi->seq[r->rid].name, &absent);
          if(absent) {
            // we have to copy the target name because it will be overwritten by the next alignment
            char *t = malloc(strlen(mi->seq[r->rid].name)+1);
            t[strlen(mi->seq[r->rid].name)] = '\0';
            strcpy(t, mi->seq[r->rid].name);
            kh_key(t2e, bin) = t;
            kh_val(t2e, bin).ct = 0;
          }
          kh_val(t2e, bin).ct = kh_val(t2e, bin).ct + (1.0/n_best_matches);
          free(r->p);
        }

        n++;
        free(reg);
      }
      mm_tbuf_destroy(tbuf);
      mm_idx_destroy(mi);
      kseq_destroy(ks); // close the query file
      gzclose(f);
    }
    mm_idx_reader_close(r); // close the index reader
  }

  const char* t;
  fprintf(stdout, "transcript\treads\tTPM\n");
  // print read counts per transcript
  double tot = 0.0;
  for (bin = 0; bin < kh_end(t2e); ++bin) {
    if (kh_exist(t2e, bin)) {
      tot = tot + kh_val(t2e, bin).ct;
    }
  }
  fprintf(stderr, "%d reads aligned\n", aligned);
  fprintf(stderr, "%lf total read ct assigned to transcripts\n", tot);
  for (bin = 0; bin < kh_end(t2e); ++bin) {
    if (kh_exist(t2e, bin)) {
      t = kh_key(t2e, bin);
      double ct = kh_val(t2e, bin).ct;
      fprintf(stdout, "%s\t%lf\t%lf\n", t, ct, ct*1000000/tot);
    }
  }

  // --- clean up memory ---

  // clean up t2e (transcript2expr)
  for (bin = 0; bin < kh_end(t2e); ++bin) {
    if (kh_exist(t2e, bin)) {
      free((char*)kh_key(t2e, bin));
    }
  }
  kh_destroy(transcript2expr, t2e);

  return 0;
}

// <o><
