import pandas as pd, numpy as np, scipy as sp
import subprocess, os, requests, re, time, zipfile, io, shutil, string, random, itertools
from featurizer_util import *
from featurizer_paths import *
from featurizer_sei import run_sei

# Note: any undefined path is probably in paths.py ; any undefined variable is probably in util.py

### Data featurization pipeline

def generate_features(features, working_dir, **args):
    args.setdefault('max_num_threads', 1)
    assert(not bool(set(args) - set(['max_num_threads', 'cell_type', 'vcf_sequence_length'])))
    vcf = working_dir + 'locs.vcf'
    bed = working_dir + 'locs.bed'
    fasta = working_dir + 'locs.fasta'

    has_vcf, has_bed, has_fasta = map(os.path.exists, [vcf, bed, fasta])
    
    assert has_vcf or has_bed or has_fasta, 'At least one of %s is required' % [vcf, bed, fasta]

    # generate bed from vcf / vcf_sequence_length if needed
    locs_df = pd.DataFrame()
    if has_vcf:
        locs_df[Vcf_header] = read_vcf(vcf)
        if not has_bed:
            assert 'vcf_sequence_length' in args, 'vcf_sequence_length must be specified as an argument if using vcf and no bed file is provided'
            L = args['vcf_sequence_length']
            # bed_pos = pos - 1; start = bed_pos - L // 2 + 1 = pos - L // 2; end = pos + L // 2
            locs_df['start'] = locs_df['pos'] - L // 2
            locs_df['end'] = locs_df['pos'] + L // 2
            write_bed(locs_df, bed, columns=Bed_header_name)
    if has_bed:
        bed_df = read_bed(bed)
        locs_df[Bed_header_name] = bed_df[Bed_header_name]
    if has_fasta:
        fasta_df = read_fasta(fasta)
        fasta_df['sequence'] = fasta_df['sequence'].apply(str.upper)
        locs_df[Fasta_header] = fasta_df[Fasta_header]
    elif has_bed or has_vcf:
        fasta_df = bed_df_to_fasta_df(locs_df)
        locs_df['sequence'] = fasta_df['sequence'].apply(str.upper)
        if has_vcf:
            locs_df['sequence'] = locs_df.apply(lambda row: substitute(row['sequence'], (row['pos'] - 1) - row['start'], row['ref'], row['alt']), axis=1)
        write_fasta(locs_df, fasta)
    
    # Columns are name, pos, start, end, sequence, ref, and alt
    # Since vcf is one-indexed and bed is zero-indexed:
    #   pos is one-indexed
    #   start and end are zero-indexed
    args['locs'] = locs_df

    # call featurization functions (run_*)

    for feature in features:
        if not os.path.exists(working_dir + feature + '.csv'):
            print('running', feature)
            eval('run_' + feature)(working_dir, args)
            print('finished', feature)
        else:
            print('already finished', feature)

from concurrent.futures import ProcessPoolExecutor

def split_fasta_(locs, output_base, get_num_splits):
    N = len(locs)
    num_splits = get_num_splits(N)
    locs_splits = np.array_split(locs, num_splits)
    for i, locs_split in enumerate(locs_splits):
        write_fasta(locs_split, output_base % i)
    return num_splits

def submit_deepsea_(input_file, output_dir):
    princeton = 'http://deepsea.princeton.edu'
    print('submitting %s to %s' % (input_file, princeton))
    with requests.session() as client:
        r = client.get(princeton + '/job/analysis/create/#')
        files = {'file' : open(input_file, 'rb')}
        data = {
            'csrfmiddlewaretoken' : client.cookies['csrftoken'],
            'submit' : 'Submit Job',
            'file_type': input_file.split('.')[-1]
        }
        r = client.post(princeton + '/job/analysis/create/#', data=data, files=files)
        time.sleep(5)
        iterator = re.compile('(/job/analysis/results/[\w-]+/)').finditer(r.text)
        try:
            output_url = princeton + next(iterator).group()
        except StopIteration:
            err = 'error running DeepSea; try manually uploading file to %s/job/analysis/create/' % (princeton)
            print(err)
            raise RuntimeError(err)
        time.sleep(5)
        r = client.get(output_url)
        while 'Job in progress' in r.text:
            print('waiting 10 seconds...')
            time.sleep(10)
            r = client.get(output_url)
        iterator = re.compile('(/media/output/[\w-]+/jobs.zip)').finditer(r.text)
        download_url = princeton + next(iterator).group()
        print('downloading and extracting to %s...' % output_dir)
        with io.BytesIO(client.get(download_url).content) as f:
            with zipfile.ZipFile(f, 'r') as z:
                z.extractall(output_dir)
        print('finished extracting')

def run_deepsea(working_dir, args):
    deepsea_dir = make_dir(working_dir + 'deepsea/')

    vcf = working_dir + 'locs.vcf'
    deepsea_bed = working_dir + 'locs_deepsea.bed'
    deepsea_fasta = working_dir + 'locs_deepsea.fasta'
    has_vcf, has_bed, has_fasta = map(os.path.exists, [vcf, deepsea_bed, deepsea_fasta])

    if has_fasta: # submit fasta
        fasta_df = read_fasta(deepsea_fasta)
        output_base = deepsea_dir + 'output_%s'
        def get_output_df_fasta(path):
            output_df = pd.read_csv(path, index_col=0)
            output_df['name'] = output_df['name'].apply(lambda x: x[1:])
            return output_df.set_index('name')
        num_splits = split_fasta_(fasta_df, output_base + '.fasta', lambda N: (N - 1) // 2000 + 1)
        input_fastas = [output_base % i + '.fasta' for i in range(num_splits)]
        output_dirs = [output_base % i + '/' for i in range(num_splits)]
        with ProcessPoolExecutor(max_workers=args['max_num_threads']) as executor:
            executor.map(submit_deepsea_, input_fastas, output_dirs)
        output_dfs = [get_output_df_fasta(output_dir + 'infile.fasta.out') for output_dir in output_dirs]
        output_df = pd.concat(output_dfs, axis=0)
    elif has_bed or not has_vcf: # submit bed
        if has_bed:
            in_file = deepsea_bed
        else:
            in_file = working_dir + 'locs.bed'
        deepsea_output = deepsea_dir + 'infile.bed.out'
        if not os.path.exists(deepsea_output):
            submit_deepsea_(in_file, deepsea_dir)
        else:
            print('Using previously deepsea server output %s' % deepsea_output)
        output_df = pd.read_csv(deepsea_output, index_col=0)
        # deepsea does not return name if bed file is submitted, we have to determine the names ourselves by sorting all the average positions
        output_df['center'] = (output_df['start'] + output_df['end']) / 2
        locs_df = args['locs'].copy()
        locs_df['center'] = (locs_df['start'] + locs_df['end']) / 2
        output_df.sort_values(['chr', 'center'], inplace=True)
        locs_sorted = locs_df.sort_values(['chr', 'center']).reset_index(drop=True)
        output_df['name'] = locs_sorted['name']
        delete_columns(output_df, Bed_header + ['center']).set_index('name', inplace=True)
    else: # submit vcf
        deepsea_output = deepsea_dir + 'infile.vcf.out.alt'
        if not os.path.exists(deepsea_output):
            submit_deepsea_(vcf, deepsea_dir)
        else:
            print('Using previously deepsea server output %s' % deepsea_output)
        output_df = pd.read_csv(deepsea_output, index_col=0)
        output_df = args['locs'][Vcf_header].merge(delete_columns(output_df, ['name']), on=['chr', 'pos', 'ref', 'alt'], how='inner')
        delete_columns(output_df, ['chr', 'pos', 'ref', 'alt']).set_index('name', inplace=True)
    output_df.reindex(index=args['locs']['name']).to_csv(working_dir + 'deepsea.csv', float_format='%.5g')

def run_5mer(working_dir, args):
    k = 5
    all_kmers = [''.join(reversed(tup)) for tup in itertools.product('ACGT', repeat=k)]
    kmer_to_index = { kmer : index for index, kmer in enumerate(all_kmers) }
    rows = []
    def count_kmers(sequence):
        counts = [0] * len(all_kmers)
        for i in range(len(sequence) - k + 1):
            counts[kmer_to_index[sequence[i : i + k]]] += 1
        return counts
    kmer_counts = []
    for i, row in args['locs'].iterrows():
        sequence = row['sequence'].upper()
        if 'N' not in sequence:
            kmer_counts.append(count_kmers(sequence))
    df = pd.DataFrame(kmer_counts, columns=all_kmers, index=args['locs']['name'])
    df.to_csv(working_dir + '5mer.csv')

def run_deepbind(working_dir, args):
    deepbind_dir = make_dir(working_dir + 'deepbind/')
    locs = args['locs']

    name_base = deepbind_dir + 'locs_%s'
    num_splits = split_fasta_(locs, name_base + '.fasta', lambda N: max(1, min(N // 20, args['max_num_threads'])))
    name_bases = [name_base % i for i in range(num_splits)]
    print('parallelizing %s jobs for %s' % (num_splits, ', '.join(name_bases)))
    def run_batch(i):
        name_base = name_bases[i]
        cmd = 'cat %s | %s %s > %s' % (name_base + '.fasta', Cmd_deepbind, Deepbind_ids, name_base + '.tab')
        return shell(cmd, wait=False)
    parallel_execution(run_batch, num_splits, ignore_outputs=True, error_msg='error in parallelizing deepbind')
    df = pd.concat([pd.read_csv(name_base + '.tab', sep='\t') for name_base in name_bases], axis=0)
    df.index = locs['name']
    df.to_csv(working_dir + 'deepbind.csv', float_format='%.5g')

def run_epigenetic(working_dir, args):
    epi_dir = make_dir(working_dir + 'epigenetic/')
    locs_unique = args['locs'].drop_duplicates(subset=Bed_header)
    unique_bed = epi_dir + 'locs_unique.bed'
    write_bed(locs_unique, unique_bed, columns=Bed_header)
    if not os.path.exists(epi_dir + 'Overlap_Table.tsv'):
        cmd = 'cd %s; %s --annotation_list %s --input_bed %s --out %s --threads %s --additional_text_labels ENCODEMarks --additional_text hg19_data' % (epi_dir, Cmd_bed_overlap, Epi_annotations, unique_bed, epi_dir, args['max_num_threads'])
        shell(cmd)
    epi_df = pd.read_csv(epi_dir + 'Overlap_Table.tsv', sep='\t', index_col=None)
    epi_df = epi_df.drop_duplicates(subset=Bed_header).set_index(Bed_header)
    epi_df = epi_df[[col for col in epi_df.columns if col.endswith('_OverlapBP')]]
    locs = args['locs'].set_index(Bed_header)[['name']]
    locs.merge(epi_df, left_index=True, right_index=True) \
        .set_index('name').loc[locs['name']] \
        .to_csv(working_dir + 'epigenetic.csv', float_format='%.5g')

def run_fimo_(working_dir, locs, motif, max_num_threads):
    fimo_dir = make_dir(working_dir + 'fimo_%s/' % motif)
    fimo_path = fimo_dir + 'fimo.txt'
    if os.path.exists(fimo_path):
        print('already ran fimo')
        return pd.read_csv(fimo_path, sep='\t', index_col=False)
    print('running fimo')
    motif_file = Motif_files[motif]

    name_base = fimo_dir + 'locs_%s'
    num_splits = split_fasta_(locs, name_base + '.fasta', lambda N: max(1, min(N // 20, max_num_threads)))
    name_bases = [name_base % i for i in range(num_splits)]
    print('parallelizing %s jobs for %s' % (num_splits, ', '.join(name_bases)))
    def run_batch(i):
        name_base = name_bases[i]
        cmd = '%s --thresh 0.0001 -verbosity 1 --oc %s %s %s' % (Cmd_fimo, name_base + '/', motif_file, name_base + '.fasta')
        return shell(cmd, wait=False)
    parallel_execution(run_batch, num_splits, ignore_outputs=True, error_msg='error in parallelizing fimo')
    fimo_df = pd.concat([pd.read_csv(name_base + '/fimo.tsv', sep='\t', index_col=False) for name_base in name_bases], axis=0) 
    #HERE: Jiayi: fimo.txt changed into fimo.tsv; read_csv to 
    fimo_df.to_csv(fimo_path, sep='\t', index=False)
    print('finished fimo')
    return fimo_df

def run_encode_matrix(working_dir, args):
    locs = args['locs']
    with open(Motif_files['encode'], 'r') as f:
        motif_list = [line.rstrip().split(' ')[1] for line in f if line.startswith('MOTIF')]
    motif_df = locs['name'].to_frame().set_index('name')
    for motif in motif_list:
        motif_df[motif] = 0
    fimo_df = run_fimo_(working_dir, locs, 'encode', args['max_num_threads'])
    for i, row in fimo_df.iterrows():
        motif_df.loc[row['sequence_name'], row['motif_id']] += row['stop'] - row['start'] + 1 #HERE
    motif_df.to_csv(working_dir + 'encode_matrix.csv', float_format='%.5g')

def find_max_overlap(intervals, padding):
    starts, ends = zip(*intervals)
    starts = sorted(starts)
    ends = np.array(sorted(ends)) + padding
    i_start = i_end = 0
    max_overlap = 0
    curr = 0
    while i_start < len(starts) and i_end < len(ends):
        if ends[i_end] <= starts[i_start]:
            curr -= 1
            i_end += 1
        else:
            curr += 1
            i_start += 1
            max_overlap = max(max_overlap, curr)
    return max_overlap

def find_max_window(starts):
    max_window = 0
    i_start = i_stop = 0
    N = len(starts)
    while i_stop < N:
        while i_stop < N and starts[i_stop] < starts[i_start] + 20:
            i_stop += 1
        max_window = max(max_window, i_stop - i_start)
        i_start += 1
    return max_window

def run_fimo_summary(working_dir, args):
    cutoff = 1e-4
    locs = args['locs']
    motifs = ['encode', 'hg19']
    features = ['motifs', 'max_window']
    summary_df = pd.DataFrame(data=0, index=locs['name'], columns=['_'.join((m, f)) for m in motifs for f in features])
    for motif in 'encode', 'hg19':
        motif_columns = ['_'.join((motif, f)) for f in features]
        fimo_df = run_fimo_(working_dir, locs, motif, args['max_num_threads'])
        for name, sequence_group in fimo_df[fimo_df['p-value'] < cutoff].groupby('sequence_name'):
            motif_locs = [motif_group.loc[motif_group['p-value'].idxmin(), 'start'] for _, motif_group in sequence_group.groupby('motif_id')]#HERE
            max_window = find_max_window(motif_locs)
            summary_df.loc[name, motif_columns] = (len(motif_locs), max_window)
    summary_df.to_csv(working_dir + 'fimo_summary.csv', float_format='%.5g')

def run_polyA_polyT_GC(working_dir, args):
    locs = args['locs']
    def count_max_consecutive(sequence, base, max_consecutive=50):
        return min(max_consecutive, max([0] + [sum(1 for x in group) for c, group in itertools.groupby(sequence) if c == base]))
    df = locs['sequence'].apply(str.upper).apply(lambda seq: pd.Series({
        'polyA' : count_max_consecutive(seq, 'A'),
        'polyT' : count_max_consecutive(seq, 'T'),
        'GC' : sum(x in ['G', 'C'] for x in seq)
    }))
    df.index = locs['name']
    df.to_csv(working_dir + 'polyA_polyT_GC.csv', float_format='%.5g')

def run_dna_shape(working_dir, args):
    dna_shape_dir = make_dir(working_dir + 'dna_shape/')
    rohs_lab = 'https://rohslab.cmb.usc.edu/DNAshape/'
    types = ['HelT', 'MGW', 'ProT', 'Roll']
    if not all(map(lambda ext: list_dir(dna_shape_dir, ext), types)):
        with requests.session() as client:
            files = {'seqfile' : open(working_dir + 'locs.fasta', 'rb')}
            data = {
                'delimiter' : '1',
                'entriesPerLine' : '20',
                'submit_button' : 'Submit',
            }
            print('querying %s server' % rohs_lab)
            r = client.post(rohs_lab + 'serverBackend.php', data=data, files=files, verify=False)
            p = re.compile(r'(Download.php\?filename=/tmp/\w+.zip)')
            iterator = p.finditer(r.text)
            try:
                output_url = rohs_lab + next(iterator).group()
            except StopIteration:
                raise RuntimeError('error running DNA shape prediction, try manually uploading file to %s to check correct format' % (rohs_lab))
            print('extracting from %s' % output_url)
            time.sleep(1)
            with io.BytesIO(client.get(output_url).content) as f:
                with zipfile.ZipFile(f, 'r') as z:
                    z.extractall(dna_shape_dir)
            print('finished extracting file')
    types_means = []
    for t in types:
        file_path = list_dir(dna_shape_dir, t)[0]
        names = []
        t_means = []
        with open(file_path, 'r') as f:
            line = f.readline()
            while line:
                names.append(line.rstrip()[1:])
                values = []
                line = f.readline()
                while line and line[0] != '>':
                    values.extend([float(val) for val in line.rstrip().split(',') if val != 'NA'])
                    line = f.readline()
                t_means.append(np.mean(values))
        types_means.append(t_means)
    df = pd.DataFrame(data=np.array(types_means).T, index=names, columns=types).loc[args['locs']['name']]
    df.to_csv(working_dir + 'dna_shape.csv', float_format='%.5g')

def run_conservation(working_dir, args):
    conservation_dir = make_dir(working_dir + 'conservation/')
    output_file = conservation_dir + 'output.tab'
    shell('%s %s %s %s' % (Cmd_big_wig_average, Phylo_annotations, working_dir + 'locs.bed', output_file), ignore_error=True)
    df = pd.read_csv(output_file, sep='\t', index_col=0, usecols=[0, 5], names=['name', 'conservation'])
    df.to_csv(working_dir + 'conservation.csv', float_format='%.5g')

def get_genes_tpm(cell_type):
    cell_paths = [path for name, path in list_dir(Gene_expressions, 'txt', return_name=True) if name.lower().startswith(cell_type.lower())]
    assert len(cell_paths) > 0, 'Need to download gene quantification file for cell type %s from https://www.encodeproject.org/experiments/' % cell_type
    genes_tpm = pd.concat([pd.read_csv(path, sep='\t', usecols=['gene_id', 'TPM'], index_col='gene_id') for path in cell_paths], axis=1)
    genes_tpm['TPM_mean'] = genes_tpm.apply(np.mean, axis=1)

    # change index from gene_id to gene_symbol
    ensembl_id_to_symbol = pd.read_csv(Gene_symbols, sep='\t', index_col='Ensembl_ID')
    genes_tpm['symbol'] = genes_tpm.index.map(lambda gene_id: ensembl_id_to_symbol['Gene'].get(gene_id.split('.')[0], np.NaN))
    genes_tpm = genes_tpm[~pd.isnull(genes_tpm['symbol'])].set_index('symbol')
    genes_tpm = genes_tpm[~genes_tpm.index.duplicated(keep='first')] # IRF9 is mapped twice, keep first
    return genes_tpm

def get_deepbind_id_symbol():
    deepbind_genes = pd.read_csv(Deepbind_ids, sep=' ', names=['gene_id', 'gene_symbol'], index_col='gene_id')
    deepbind_genes['gene_symbol'] = deepbind_genes['gene_symbol'].apply(lambda s: s.lstrip('#'))
    return deepbind_genes

def get_deepbind_annotation():
    db_annotation = pd.read_csv(Resources + 'deepbind_tf_annotations.txt', sep='\t', index_col=0)
    db_annotation['gene_symbol'] = db_annotation['Protein']
    import re
    def search(details):
        res = re.search('CloneType=([\w-]+)', details)
        if res is None:
            res = re.search('CellLine=([\w-]+)', details)
        return res.group(1)
    db_annotation['cell_type'] = db_annotation['Experiment Details'].apply(search)
    return db_annotation[['gene_symbol', 'cell_type']]

def get_epigenetic_tf_annotation():
    annotation = pd.DataFrame(columns=['cell_type', 'tf'])
    with open(Resources + 'encode_tf_annotations.txt', 'r') as f:
        for line in f:
            line = line.rstrip()
            arr = line.split('\t')
            name = arr[0].split('.')[0]
            details = dict(x.split('=') for x in arr[1].split('; '))
            annotation.loc[name, ['cell_type', 'tf']] = details['cell'], details['antibody']
    return annotation

def run_anonymize_tf(working_dir, args):
    genes_tpm = get_genes_tpm(args['cell_type'])

    deepbind_genes = get_deepbind_id_symbol()
    deepbind_genes['TPM'] = deepbind_genes['gene_symbol'].apply(lambda s: genes_tpm['TPM_mean'].get(s, np.NaN))

    # tf anonymization
    q20, q40, q60, q80 = deepbind_genes['TPM'].quantile([0.2, 0.4, 0.6, 0.8])
    top_middle_bottom = np.array(deepbind_genes['TPM'].apply(
        lambda t: pd.Series([t <= q20, q40 < t <= q60, q80 < t])
    ))
    top_middle_bottom_counts = np.sum(top_middle_bottom, axis=0)
    deepbind_path = working_dir + 'deepbind.csv'
    print('running deepbind to quantify binding for sequences')
    if os.path.exists(deepbind_path):
        print('already finished')
    else:
        run_deepbind(working_dir, args)
        print('finished deepbind')
    deepbind_df = pd.read_csv(deepbind_path, index_col=0)
    df_quantiles = deepbind_df.quantile(0.9)
    # number of [low, medium, high] gene expressions with active binding for the sequence
    tf_bindings = deepbind_df.apply(lambda r: pd.Series( # r is (515,); top_middle_bottom is (515, 3)
        data=np.sum(np.expand_dims(r > df_quantiles, axis=1) & top_middle_bottom, axis=0),
        index=['num_low_tf', 'num_med_tf', 'num_high_tf']
    ), axis=1)
    tf_bindings.index = deepbind_df.index
    tf_bindings.to_csv(working_dir + 'anonymize_tf.csv', float_format='%.5g')

def run_closest_gene(working_dir, args):
    closest_gene_dir = make_dir(working_dir + 'closest_gene/')
    great_root = 'http://great.stanford.edu/public/'
    great_output_path = closest_gene_dir + 'great.txt'
    with requests.session() as client:
        files = { 'fgFile' : open(working_dir + 'locs.bed', 'rb') }
        data = {
            'species' : 'hg19',
            'rule' : 'oneClosest',
            'fgChoice' : 'file',
            'bgChoice' : 'wholeGenome',
        }
        cgi_bin = great_root + 'cgi-bin/'
        res = client.post(cgi_bin + 'greatWeb.php', data=data, files=files)
        p = re.compile(r'(showAllDetails.php\?.*foreName=locs.bed)')
        iterator = p.finditer(res.text)
        download_url = cgi_bin + next(iterator).group().replace('showAllDetails', 'downloadAssociations') + '&table=region'
        res = client.get(download_url)
        with open(great_output_path, 'wb') as f:
            f.write(res.content)
    genes_tpm = get_genes_tpm(args['cell_type'])
    closest_gene = pd.read_csv(great_output_path, sep='\t', skiprows=1, names=['name', 'gene'], index_col='name')
    closest_gene['TPM_closest_gene'] = closest_gene['gene'].apply(lambda g: genes_tpm['TPM_mean'].get(g.split(' ')[0], 0))
    delete_columns(closest_gene, ['gene']).to_csv(working_dir + 'closest_gene.csv', float_format='%.5g')

def run_intron_exon_promoter(working_dir, args):
    eip_dir = make_dir(working_dir + 'intron_exon_promoter/')
    overlap_columns = ('intron', 'exon', 'promoter')
    eip_df = pd.DataFrame(0, index=args['locs']['name'], columns=overlap_columns)
    for label, ref_file in zip(overlap_columns, (Introns, Exons, Promoters)):
        overlap_file = eip_dir + '%s_overlap.bed' % label
        cmd = '%s intersect -a %s -b %s -u > %s' % (Cmd_bedtools, working_dir + 'locs.bed', ref_file, overlap_file)
        shell(cmd)
        overlap_bed = read_bed(overlap_file)
        eip_df.loc[overlap_bed['name'], label] = 1
    eip_df['distal'] = ((eip_df['intron'] | eip_df['exon'] | eip_df['promoter']) == 0).astype(int)
    eip_df.to_csv(working_dir + 'intron_exon_promoter.csv', float_format='%.5g')