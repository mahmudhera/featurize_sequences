# Script paths
Cmd_bedtools = '/home/jl2791/.conda/envs/work/bin/bedtools'#'singularity exec --cleanenv --bind /home/$USER,/scratch/$USER,/projectsp/f_ak1833_1,/projects,/tmp /projects/community/singularity.images/NGSgadgets/bedtools/bedtools_v2.30.0.sif bedtools'#'/home/jl2791/.conda/envs/jlenv/bin/bedtools'
Cmd_deepbind = '/scratch/jl2791/programs/deepbind/deepbind'
Cmd_bed_overlap = 'python3 /home/jl2791/scripts/Overlaps2/Overlaps2.py'
Cmd_fimo = '/home/jl2791/.conda/envs/work/bin/fimo'
Cmd_big_wig_average = '/home/jl2791/.conda/envs/work/bin/bigWigAverageOverBed'
#### Jiayi Added: featurizer_paths.py
Cmd_sei_api = 'python3 /home/jl2791/scripts/featurizer/DeepSea_sei_API.py'
Cmd_Rscript = "/home/jl2791/scripts/featurizer/bin/rscript_singularity.sh"
DNAshapeR_cli = "/home/jl2791/scripts/featurizer/DNAshapeR_cli.R"



# Large file paths
#Genome = '/projectsp/f_ak1833_1/jliu/data/genome/hg19.fa'
Genome = '/scratch/jl2791/data/genome/staging/hg19.fa'

# Project subpaths for convenience
Resources = '/home/jl2791/scripts/featurizer/resources/'

Phylo_annotations = Resources + '/hg19.100way.phyloP100way.bw'

# Files in project subpaths
Chrom_sizes = Resources + 'hg19.chrom.sizes.txt'
Gene_symbols = Resources + 'ENSG_GENEsymbol.txt'
Deepbind_ids = Resources + 'HomoSapiens_SELEX_ChIP.IDS'
Motif_files = {
    'encode' : Resources + 'encode_motifs.meme',
    'hg19' : Resources + 'hg19_motifs.meme'
}
Epi_annotations = Resources + 'all_annotations.tab'
Gene_expressions = Resources + 'gene_expressions/'
Genes = Resources + 'hgTables_genes.tab'
Exons = Resources + 'hgTables_exons.bed'
Introns = Resources + 'hgTables_introns.bed'
Promoters = Resources + 'hgTables_tx_-2500_500.bed'