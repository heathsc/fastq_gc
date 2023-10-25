# fastq_gc
Utility to estimate GC distribution from FASTQ file.

# Changes
0.2.3 - Add output file option.  Require >1 non-overlapping hash hit for control sequences
0.2.2 - Move to 27-mers for bisulfite analysis (2-mers for regular sequences)
0.2.1 - Check that read end is known if stranded bisulfite type specified
0.2.0 - Change output to use serde/serde_json
0.1.1 - Allow specification of bisulfite type
0.1.0 - Initial commit