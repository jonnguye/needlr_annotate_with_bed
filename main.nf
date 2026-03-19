nextflow.enable.dsl=2

// --- Parameters ---
// You can override these via command line: --input_tsvs 'path/*.txt'
params.input_tsvs = "/path/to/your/input_tsvs/*.tsv"
params.output_dir = "final_results"

// ID-based: [ file, [indices], [labels], prefix ]
params.id_annotations = [
  [ file("/oak/stanford/groups/smontgom/jonnguye/cardioid_sv/pipeline/linked_peaks.bed.gz"),
    [0,1,2,3,6],
    ["Chr","Start","End","Gene","Cluster"],
    "CARDIOID_ATAC"
  ]
]

// Coord-based: [ file, label_prefix ]
params.coord_annotations = [
    [file("/oak/stanford/groups/smontgom/jonnguye/cardioid_sv/pipeline/002_1-fire-v0.1.1-peaks.bed.gz"), "002_1_FIRE"],
    [file("/oak/stanford/groups/smontgom/jonnguye/cardioid_sv/pipeline/002_3-fire-v0.1.1-peaks.bed.gz"), "002_3_FIRE"]
]

process ANNOTATE_WITH_IDS {
    tag "ID: ${sample_id}-${prefix}"

    input:
    tuple val(sample_id), path(tsv), path(bed), val(col_indices), val(col_labels), val(prefix)

    output:
    tuple val(sample_id), path("tmp_${sample_id}_${prefix}.tsv")

    script:
    // bedtools uses 1-based column indices; params are 0-based
    def indices = col_indices.collect { it + 1 }.join(',')
    def ops     = col_indices.collect { 'collapse' }.join(',')

    """
    ml py-pybedtools/0.9.0_py39 py-pandas/2.0.1_py39
    
    cat > annotate_ids.py << 'SCRIPT_END'
import pandas as pd
import pybedtools
import sys

tsv_file = '${tsv}'
bed_file = '${bed}'
indices = "${indices}"
ops = "${ops}"
sample_id = "${sample_id}"
prefix = "${prefix}"

# Load TSV
df = pd.read_csv(tsv_file, sep='\\t', header=0)
df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')
df.iloc[:, 2] = pd.to_numeric(df.iloc[:, 2], errors='coerce')
df = df.sort_values(by=[df.columns[0], df.columns[1], df.columns[2]])

# Use BED3 only for -a to avoid any weirdness with extra columns
a_df = df.iloc[:, :3].copy()
a_df.columns = ["chrom","start","end"]
a = pybedtools.BedTool.from_dataframe(a_df)

b = pybedtools.BedTool(bed_file).sort()

# Collapse BED fields we need (chr,start,end,gene,cluster)
res_bt = a.map(b, c=indices, o=ops)
res_df = res_bt.to_dataframe(disable_auto_names=True, header=None)

# res_df columns: A_chrom,A_start,A_end, collapsed_fields...
chr_col     = 3 + 0
start_col   = 3 + 1
end_col     = 3 + 2
gene_col    = 3 + 3
cluster_col = 3 + 4

def split_list(x):
    if pd.isna(x) or x == '.':
        return []
    return str(x).split(',')

hits = []
counts = []

for i in range(len(res_df)):
    chrs = split_list(res_df.iat[i, chr_col])
    sts  = split_list(res_df.iat[i, start_col])
    ens  = split_list(res_df.iat[i, end_col])
    ges  = split_list(res_df.iat[i, gene_col])
    cls  = split_list(res_df.iat[i, cluster_col])

    n = min(len(chrs), len(sts), len(ens), len(ges), len(cls))
    if n == 0:
        hits.append('None')
        counts.append(0)
        continue

    entries = [f"{chrs[j]}:{sts[j]}-{ens[j]},{ges[j]},{cls[j]}" for j in range(n)]
    hits.append(';'.join(entries))
    counts.append(len(entries))

final_df = df.copy()
final_df[f"{prefix}_Hits"] = hits
final_df[f"{prefix}_Count"] = counts

final_df.to_csv(f"tmp_{sample_id}_{prefix}.tsv", sep='\\t', index=False)
SCRIPT_END

    python3 annotate_ids.py
    """
}

process ANNOTATE_BY_COORDS {
    tag "Coord: ${sample_id}-${prefix}"
    
    input:
    tuple val(sample_id), path(tsv), path(bed), val(prefix)

    output:
    tuple val(sample_id), path("tmp_${sample_id}_${prefix}.tsv")

    script:
    """
    ml py-pybedtools/0.9.0_py39 py-pandas/2.0.1_py39
    
    cat > annotate_coords.py << 'SCRIPT_END'
import pandas as pd
import pybedtools

tsv_file = '${tsv}'
bed_file = '${bed}'
sample_id = "${sample_id}"
prefix = "${prefix}"

# Load and Sort TSV
df = pd.read_csv(tsv_file, sep='\\t')
df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')
df.iloc[:, 2] = pd.to_numeric(df.iloc[:, 2], errors='coerce')
df = df.sort_values(by=[df.columns[0], df.columns[1], df.columns[2]])

a = pybedtools.BedTool.from_dataframe(df)

# Load messy BED, create coord-ID, and sort
b_df = pd.read_csv(bed_file, sep='\\t', header=None, usecols=[0,1,2], names=['c','s','e'])
b_df['id'] = b_df['c'].astype(str) + ':' + b_df['s'].astype(str) + '-' + b_df['e'].astype(str)
b = pybedtools.BedTool.from_dataframe(b_df).sort()

# Map using our generated ID (column 4)
res_bt = a.map(b, c=4, o=['distinct', 'count'])
res_df = res_bt.to_dataframe(disable_auto_names=True, header=None)

final_df = res_df.iloc[:, :len(df.columns) + 2]
final_df.columns = list(df.columns) + [f"{prefix}_Coords", f"{prefix}_Count"]

final_df[f"{prefix}_Coords"] = final_df[f"{prefix}_Coords"].replace({'.': 'None'})
final_df[f"{prefix}_Count"] = pd.to_numeric(final_df[f"{prefix}_Count"]).fillna(0).astype(int)

final_df.to_csv(f"tmp_{sample_id}_{prefix}.tsv", sep='\\t', index=False)
SCRIPT_END

    python3 annotate_coords.py
    """
}

process JOIN_ANNOTATIONS {
    tag "Join: ${sample_id}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(files)

    output:
    path "${sample_id}_final_annotated.tsv"

    script:
    """
    ml py-pandas/2.0.1_py39
    
    cat > join_annotations.py << 'SCRIPT_END'
import pandas as pd

file_list = "${files}".split()
sample_id = "${sample_id}"

# Start with the first file
main_df = pd.read_csv(file_list[0], sep='\\t')

for f in file_list[1:]:
    curr_df = pd.read_csv(f, sep='\\t')
    # Join only the new annotation/count columns
    new_cols = [c for c in curr_df.columns if c not in main_df.columns]
    main_df = pd.concat([main_df, curr_df[new_cols]], axis=1)

main_df.to_csv(f"{sample_id}_final_annotated.tsv", sep='\\t', index=False)
SCRIPT_END

    python3 join_annotations.py
    """
}

// --- Workflow ---

workflow {
    // 1. Channel of [baseName, path]
    tsv_ch = Channel.fromPath(params.input_tsvs)
        .map { it -> [it.baseName, it] }

    // 2. Parallel Annotations
    // Using .combine() creates a flat tuple for the process input
    id_results = ANNOTATE_WITH_IDS(
        tsv_ch.combine(Channel.fromList(params.id_annotations))
    )
    
    coord_results = ANNOTATE_BY_COORDS(
        tsv_ch.combine(Channel.fromList(params.coord_annotations))
    )

    // 3. Merge streams and group by Sample ID
    id_results.mix(coord_results)
        .groupTuple()
        .set { grouped_ch }

    JOIN_ANNOTATIONS(grouped_ch)
}
