import argparse
import os

def finemap_dapg(args):
    # Read files
    z = open(args.dapg_z, 'r')
    zf = z.readlines()
    zf = [x.rstrip('\n') for x in zf]
    ld = open(args.dapg_ld, 'r')
    ldf = ld.readlines()
    ldf = [x.rstrip('\n') for x in ldf]
    if args.loci_name != None:
        ln = open(args.loci_name, 'r')
        lnf = ln.readlines()
        lnf = [x.rstrip('\n') for x in lnf]
    if args.dapg_ens != None:
        if args.dapg_ens not in [str(x) for x in range(1,11)]:
            ens = open(args.dapg_ens, 'r')
            ensf = ens.readlines()
            ensf = [x.rstrip('\n') for x in ensf]
    # run
    for l in range(args.loci[0], args.loci[1]):
        if 'lnf' in locals().keys():
            fn = f'{lnf[l-1]}'
        else:
            fn = f'locus{l}'
        locusid = l
        opdir = f'{args.output}/{fn}'
        if not os.path.exists(opdir):
            os.makedirs(opdir)

        dapg_co = f'{args.dapg} -d_z {zf[l-1]} -d_ld {ldf[l-1]}'
        if args.dapg_log:
            dapg_co = dapg_co +f' -l {opdir}/{fn}.log'
        if args.dapg_msize != None:
            dapg_co = dapg_co +f' -msize {args.dapg_msize}'
        if args.dapg_t != None:
            dapg_co = dapg_co +f' -t {args.dapg_t}'
        if args.dapg_ldcont != None:
            dapg_co = dapg_co +f' -ld_control {args.dapg_ldcont}'
        if args.dapg_ens != None:
            if args.dapg_ens in [str(x) for x in range(1,11)]:
                dapg_co = dapg_co +f' -ens {args.dapg_ens}'
            else:
                dapg_co = dapg_co +f' -ens {ensf[l-1]}'

        dapg_co = dapg_co + f' -o {opdir}/{fn}.output'
        print('running dapg using command below')
        print(dapg_co)
        os.system(dapg_co)

        # summarise dapg_results
        r1 = os.system(f'grep \"\[\" {opdir}/{fn}.output > {opdir}/{fn}.model.summary')
        r2 = os.system(f'grep \"((\" {opdir}/{fn}.output > {opdir}/{fn}.snp.summary')
        r3 = os.system('grep \"\{\" %s/%s.output > %s/%s.signal.cluster.summary'%(opdir, fn, opdir, fn))
        for i in [r1, r2, r3]:
            if i != 0:
                print(i)
        # get credible sets information
        if args.cs_prob!=None and args.getcs_script != None:
            r4 = os.system(f'perl {args.getcs_script} -d {opdir}/{fn}.output -p {args.cs_prob} > {opdir}/{fn}.{args.cs_prob}cs')
            if r4 != 0:
                print (r4)
        if r2 or r3 or r4:
            print (f'Error! Locus {fn}')
        else:
            print(f'Finemap for {fn} done!')
    return None

def main(args):
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    finemap_dapg(args)
    return 'Finemapping using dapg done'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-dapg', required=True, type=str, help='path for excutable dapg')
    parser.add_argument(
        '-loci', required=True, type=int, nargs=2,
        help="""loci for running.
        Should have two loci ids for input, the first one is for start loci(included), the second one is for end loci (not included).
        For example, if input 1 116, the script will run dapg for locus 1 to 115. 
        Id of the loci are the row number (start from 1) of the file provided through -dapg_z/-dapg_ld""")
    parser.add_argument(
        '-loci_name', type=str, 
        help="""
        Optional: for providing the name for loci, otherwise the name for output files will be the loci id, 
        i.e. the row number (start from 1) of the file provided through -dapg_z/-dapg_ld.""")
    parser.add_argument(
        '-dapg_z', required=True, type=str, 
        help="""
        a file contains zfile paths provided to dap-g; one path per row, 
        order should in consistant with the file provided for -dapg_ld""")
    parser.add_argument(
        '-dapg_ld', required=True, type=str, 
        help="""
        a file contains ldfile paths provided to dap-g; one path per row, 
        order should in consistant with the file provided for -dapg_z""")
    parser.add_argument(
        '-dapg_log', type = bool,
        help="Optional: whether use dapg -l for logging"
)
    parser.add_argument(
        '-dapg_msize', type = int,
        help="Optional: dapg msize parameter"
    )
    parser.add_argument(
        '-dapg_t', type = int,
        help="Optional: dapg -t parameter, i.e. number of cores"
    )
    parser.add_argument(
        '-dapg_ldcont', type = float,
        help="Optional: dapg -ld_control parameter"
    )
    parser.add_argument(
        '-dapg_ens', type = str,
        help="Optional: dapg -ens parameter"
    )
    parser.add_argument(
        '-output', type = str,
        help="Optional: output directory"
    )
    parser.add_argument(
        '-getcs_script', type = str,
        help="Optional: path to get_credible_set.pl (provided in the github dap/utility folder)"
    )
    parser.add_argument(
        '-cs_prob', type = float,
        help="Optional: percentage of cs to be summarized, e.g. 0.95"
    )
    args = parser.parse_args()
    main(args)