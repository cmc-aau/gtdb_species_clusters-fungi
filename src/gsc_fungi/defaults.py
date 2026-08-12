# criteria for QC'ing genomes
QC_MAX_PERC_DUP_MARKERS = 10.0
QC_MIN_PERC_MARKERS = 30.0

# criteria for defining ANI-based species clusters
ANI_SP = 95.0
AF_SP = 50.0

# minimum ANI estimate from skani sketch estimate for genome
# pair to be fully processsed
SKANI_PREFILTER_THRESHOLD = 85.0

# skani presets selecting for different compression factors (k-mer sampling rates)
# None = default which sets c to 125
SKANI_PRESET = None
