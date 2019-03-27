minnow
======

A lightweight but complete tool to generate transcript abundances from ONT Direct RNA sequences


Installation
------------

minnow uses klib and minimap2 headers and requires zlib development files, and make

You can build minnow from scratch as follows on most \*nix systems,
but minimap2 - and thus minnow - are optimized for x86-64 CPUs: 

    git clone http://github.com/txje/minnow
    cd minnow
    mkdir incl
    cd incl
    git clone http://github.com/attractivechaos/klib
    git clone http://github.com/lh3/minimap2
    cd minimap2
    make
    cd ../..
    make

Usage
-----

minnow requires a reference transcriptome (currently NOT a whole genome) and long reads.
Both may be FASTA or FASTQ and optionally gzipped.

For humans, you probably want [the RefSeq transcriptome](http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz).

    ./minnow -r transcripts.fa -q reads.fq > tpm.tsv

    Options:
      -q: FASTA/Q[.gz] file with reads (required)
      -r: Reference FASTA/Q[.gz] (required)
      -t: Threads (default: 1)
      -c: 'careful': more accurate but slower
      -v, --verbose: verbose
      -h, --help: show this
      --version: show version information

Output
------

Since most of your reads should represent full-length or nearly full-length transcripts, we consider 1 read == 1 transcript
Reads are assigned to the best-fit transcript. If multiple transcripts tie (this happens less often if you use --careful),
they are each assigned a partial read. TPM is just normalized as though there were exactly 1m reads.

    transcript      reads     TPM
    NM_001191005.2  6.533333  13.586062
    XM_011530291.3  3.000000  6.238498
    XR_001749299.1  0.166667  0.346583
    ...

Features
--------

Right now, that's all you get, but in the future I might add spliced
alignment so that a reference genome and ideally some annotation can
be provided instead of a transcriptome.
