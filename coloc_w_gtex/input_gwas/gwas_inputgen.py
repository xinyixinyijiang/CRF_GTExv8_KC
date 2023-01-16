import argparse
import pandas as pd

def gwas_format(args):
    # TODO: format to fastenloc input
    if args.sep == None:
        df = pd.read_csv(args.gwas_fm_raw, sep='\t')
    else:
        df = pd.read_csv(args.gwas_fm_raw, sep=args.sep)
    n_uniqsnps = len(df[args.snpid_col].drop_duplicates(keep='first'))
    n_lines = len(df)
    print(f'reading gwas file - finished - {n_lines}lines,{n_uniqsnps}snps')
    if args.gtexid_ref != None:
        ref = pd.read_csv(args.gtexid_ref, sep='\t')
        ref.columns = [args.snpid_col, 'gtexid']
        ref = ref.drop_duplicates(subset = [args.snpid_col], keep='first')
        df = df.merge(ref, how = 'left')
        print('gwas id switched to gtex id')
    df['signalid'] = 'locus'+df[args.locusid_col].apply(str)+'_cs'+df[args.csid_col].apply(str)
    if args.gtexid_ref == None:
        df = df[[args.snpid_col, 'signalid', args.pip_col]]
    else:
        df = df[['gtexid', 'signalid', args.pip_col]]
    df.to_csv(args.output, compression = 'gzip', sep='\t', index=False, header=None)
    print(f'fastenloc gwasinput has been output to {args.output}')
    return None

def main(args):
    gwas_format(args)
    return 'gwas_format - finished'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-gwas_fm_raw', required=True, type=str, help='GWAS raw file')
    parser.add_argument('-snpid_col', required=True, type=str, help='column name for the SNP ids')
    parser.add_argument('-locusid_col', required=True, type=str, help='column name for the locus ids')
    parser.add_argument('-csid_col', required=True, type=str, help='column name for the cs ids')
    parser.add_argument('-pip_col', required=True, type=str, help='column name for the pip')
    parser.add_argument(
        '-gtexid_ref', type=str, 
        help='File that contain two columns, the first columns is gwassnpid, the second columns is gtexsnpid, tab deliminated, should with header')
    parser.add_argument('-output', required=True, type=str, help='output gzip (thus suffix should be .gz) file name/path')
    parser.add_argument('-sep', type=str, help='if file is not tab deliminated, please specify')
    args = parser.parse_args()
    main(args)