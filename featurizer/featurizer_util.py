import numpy as np, pandas as pd
import subprocess, os, re, time, zipfile, gzip, io, shutil, string, random, itertools, pickle, json, codecs

from featurizer_paths import *

# Constants
Bed_header = ['chr', 'start', 'end']
Bed_header_name = Bed_header + ['name']
Fasta_header = ['name', 'sequence']
Vcf_header = ['chr', 'pos', 'name', 'ref', 'alt']
Chroms = sorted(['chr' + x for x in ([str(i) for i in range(1, 23)] + ['X', 'Y'])])
Bases = ['A', 'C', 'G', 'T', 'N']

def get_name(path, ext=True):
    name = os.path.basename(path)
    if ext: return name
    else: return os.path.splitext(name)[0]

def list_dir(dir_, ext, return_name=False):
    if ext == '/':
        criteria = lambda x: os.path.isdir(os.path.join(dir_, x))
        strip_ext = lambda x: x
        post_process = lambda x: x + '/'
    else:
        if ext[0] != '.': ext = '.' + ext.lower()
        criteria = lambda x: x[-len(ext):].lower() == ext
        strip_ext = lambda x: get_name(x, ext=False)
        post_process = lambda x: x
    files = (f for f in os.listdir(dir_) if criteria(f))
    files = sorted((strip_ext(file), post_process(os.path.join(dir_, file))) for file in files)
    if return_name: return files
    else: return [path for file, path in files]

def load_json(path):
    with open(path, 'rb') as f:
        return json.load(f)

def save_json(dict_, path):
    with open(path, 'wb') as f:
        json.dump(dict_, codecs.getwriter('utf-8')(f), indent=4, sort_keys=True)

def load_text(path):
    with open(path, 'r') as f:
        return f.read()

def save_text(string, path):
    with open(path, 'w+') as f:
        f.write(string)

def load_pickle(path):
    with open(path, 'rb') as f:
        return pickle.load(f)

def save_pickle(obj, path):
    with open(path, 'wb') as f:
        pickle.dump(obj, f)

def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def remove(path):
    if not os.path.exists(path):
        return
    elif os.path.isfile(path):
        os.remove(path)
    else:
        shutil.rmtree(path)

def wget(link, output_directory):
    cmd = 'wget %s -P %s' % (path, output_directory)
    shell(cmd)
    output_path = os.path.join(os.path.basename(link))
    if not os.path.exists(output_path): raise RuntimeError('Failed to run %s' % cmd)
    return output_path

def extract(input_path, output_path=None):
    if input_path[-3:] == '.gz':
        if not output_path:
            output_path = input_path[:-3]
        with gzip.open(input_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                f_out.write(f_in.read())
    else:
        raise RuntimeError('Don\'t know file extension for ' + input_path)

def get_temp_file(ext, N=20):
    if not ext.startswith('.'):
        ext = '.' + ext
    return ''.join(np.random.choice(list(string.digits), N)) + ext

def get_temp_dir(N=10):
    return get_temp_file('/', N)

def shell(cmd, wait=True, ignore_error=True):
    if type(cmd) != str:
        cmd = ' '.join(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if not wait:
        return process
    out, err = process.communicate()
    if err:
        if ignore_error == 2:
            pass
        elif ignore_error:
            print(err.decode())
        else:
            print(err.decode())
            raise RuntimeError('Error in command line call')
    return out.decode()

def parallel_execution(process_fn, n_jobs, ignore_outputs=False, error_msg=''):
    '''
    process_fn: fn that takes in the thread index (0 to n_job - 1) returns a process
    '''
    processes = [process_fn(i) for i in range(n_jobs)]
    outs, errs = zip(*[process.communicate() for process in processes])
    error = False
    for err in errs:
        if err:
            print(err.decode('UTF-8'))
            error = True
    if error: raise RuntimeError(error_msg)
    if ignore_outputs: return
    for out in outs:
        print(out.decode('UTF-8'))

def attributes(obj):
    import inspect, pprint
    pprint.pprint(inspect.getmembers(obj, lambda a: not inspect.isroutine(a)))


def substitute(seq, pos, ref, alt):
    pos_end = pos + len(ref)
    assert seq[pos : pos_end] == ref, 'sequence %s and ref %s do not match.' % (seq[pos : pos_end], ref)
    return seq[:pos] + alt + seq[pos_end:]

def resize_length(df, new_length, start='start', end='end', new_start='start', new_end='end'):
    '''
    resize the sequence around the center. Center is center for odd length, right of center for even length.
    '''
    convert_columns_to_numeric(df, [start, end])
    middle = (df[start] + df[end]) // 2
    front_length = new_length // 2
    back_length = new_length - front_length
    df[new_start] = middle - front_length
    df[new_end] = middle + back_length
    return df

def append_bed_columns(df, column):
    df[Bed_header] = df[column].str.extract('(\w+):(\d+)-(\d+)', expand=True)
    for column in 'start', 'end':
        df[column] = df[column].astype(int)
    return df

def append_bed_summary_column(df, column):
    df[column] = df.apply(lambda row: '%s:%s-%s' % (row['chr'], row['start'], row['end']), axis=1)
    return df

def merge_sequences(ref_seq, alt_seq, alt_offset):
    return ''.join([
        ref_seq[:alt_offset],
        alt_seq,
        ref_seq[alt_offset + len(alt_seq):]
    ])

def read_fasta(input_fasta, has_name=False):
    data = []
    with open(input_fasta, 'r') as f:
        line = f.readline()
        while line:
            name = line.rstrip()[1:]
            if has_name:
                name = name.split('::')[0]
            data.append((name, f.readline().rstrip()))
            line = f.readline()
    return pd.DataFrame(data, columns=Fasta_header)

def read_bed(input_bed):
    if len(read_header(input_bed, sep='\t')) == 3:
        names = Bed_header
    else:
        names = Bed_header_name
    return pd.read_csv(input_bed, sep='\t', names=names)

def read_vcf(input_vcf):
    return pd.read_csv(input_vcf, sep='\t', names=Vcf_header, keep_default_na=False)

def read_header(path, sep=None):
    with open(path, 'r') as f:
        return f.readline().rstrip().split(sep)

def write_bed(df, output_bed, columns=Bed_header):
    df[columns].to_csv(output_bed, sep='\t', index=False, header=False)

def write_fasta(df, output_fasta, columns=Fasta_header):
    with open(output_fasta, 'w+') as f:
        for i, row in df.iterrows():
            f.write('>%s\n%s\n' % tuple(row[columns]))

def write_vcf(df, output_vcf, columns=Vcf_header):
    df[Vcf_header].to_csv(output_vcf, sep='\t', index=False, header=False)

def bed_df_to_fasta_df(bed_df, columns=Bed_header_name):
    temp_bed = get_temp_file('bed')
    temp_fasta = get_temp_file('fasta')
    write_bed(bed_df, temp_bed, columns=columns)
    bed_to_fasta(temp_bed, temp_fasta)
    fasta_df = read_fasta(temp_fasta)
    os.remove(temp_bed)
    os.remove(temp_fasta)
    return fasta_df

def bed_to_fasta(bed_path, fasta_path):
    name_flag = ''
    if len(read_header(bed_path, '\t')) > 3:
        name_flag = ' -name'
    shell('%s getfasta -fi %s -bed %s -fo %s' % (Cmd_bedtools, Genome, bed_path, fasta_path) + name_flag)
    if name_flag:
        df = read_fasta(fasta_path)
        df['name'] = df['name'].apply(lambda x: x.split('::')[0])
        write_fasta(df, fasta_path)

def convert_columns_to_numeric(df, columns):
    df[columns] = df[columns].apply(pd.to_numeric)
    return df

def delete_columns(df, columns):
    for column in columns:
        del df[column]
    return df

def append_bed_columns(df, column):
    df[Bed_header] = df[column].str.extract('(\w+):(\d+)-(\d+)', expand=True)
    for column in 'start', 'end':
        df[column] = df[column].astype(int)
    return df
