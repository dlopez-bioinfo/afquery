from config import *  # noqa: F401, F403

# Override scales for smoke test (tiny)
SYNTH_SCALES = [100]  # 1 scale only
BUILD_THREAD_COUNTS = [1, 4]  # 2 thread counts
ANNOTATE_THREAD_COUNTS = [1, 4]
ANNOTATE_VARIANT_COUNTS = [100]
QUERY_COLD_REPS = 1
QUERY_WARM_REPS = 3
QUERY_REGION_REPS = 2
QUERY_BATCH_REPS = 2
BCFTOOLS_REPS = 1
