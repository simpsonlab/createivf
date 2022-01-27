import sys

#sample_metadata = "/.mounts/labs/ont/createivf/analysis/metadata/sample_cytogenetics.tsv"
sample_metadata = "/.mounts/labs/simpsonlab/sw/createivf/bin/metadata/sample_cytogenetics.tsv"

target_sample = sys.argv[1]
bp_id = int(sys.argv[2])

# Get breakpoints from metadata file
target = ""
with open(sample_metadata) as f:
    for line in f:
        (sampleid, region1, region2) = line.rstrip().split()

        # Handle targets like PGD15_mother
        if target_sample.find(sampleid) != -1:
            if bp_id == 1:
                target = region1
            else:
                target = region2
            break

# add chr prefix
target = "chr" + target

bp_chr = list()
bp_coords = list()

#cytobands = "/.mounts/labs/ont/createivf/analysis/metadata/hg38.cytoBand.composite.txt"
cytobands = "/.mounts/labs/simpsonlab/sw/createivf/bin/metadata/hg38.cytoBand.composite.txt"
with open(cytobands) as f:
    for line in f:
        
        fields = line.rstrip().split()
        if len(fields) != 6:
            continue

        (chrom, start, end, _, _, key) = fields

        start = int(start)
        end = int(end)

        # handle ambiguous regions
        if key.find(target) != -1:
            bp_chr.append(chrom)
            bp_coords.append(start)
            bp_coords.append(end)
        
bp_coords = sorted(bp_coords)

chrom = bp_chr[0]
start = bp_coords[0]
end = bp_coords[-1]

print "%s:%d-%d" % (chrom, start, end)

