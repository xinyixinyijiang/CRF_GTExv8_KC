import os
import pandas as pd

def summary_dapg_res(res, csprob):
    fl = os.listdir(res)
    output = []
    no_cs = []
    for d in fl:
        try:
            df = pd.read_csv(f'{res}/{d}/{d}.{csprob}cs', sep='\s+', header=None)
            df.columns = ['cs','snp', 'pip']
            df['locus'] = d
            cslist = list(set(df['cs']))
            cslist.sort()
            df['cs'] = df['cs'].replace(cslist, range(1,len(set(df['cs']))+1))
            output.append(df)
        except pd.errors.EmptyDataError:
            print(d, f'no {csprob} cs')
            no_cs.append(d)
    output = pd.concat(output)
    output['chr'] = output['locus'].apply(lambda x: int(x.split('-')[1].strip('chr')))
    output['start'] = output['locus'].apply(lambda x: int(x.split('-')[2]))
    output['end'] = output['locus'].apply(lambda x: int(x.split('-')[3]))
    output['locus'] = output['locus'].apply(lambda x: int(x.split('-')[0].strip('locus')))
    output.to_csv(f'./{res}_pip_summary.tsv', sep='\t', index=False)
    if len(no_cs) > 0:
        f = open(f'{res}_without_{csprob}cs.txt', 'w')
        for i in no_cs:
            f.write(i+'\n')
        f.close()
    print(f'summary_dapg_res in {res} with csprob {csprob} - finished')
    print(f'file wrote to ./{res}_pip_summary.tsv')
    if len(no_cs) > 0:
        print(f'file wrote to .{res}_without_{csprob}cs.txt')
    return  
