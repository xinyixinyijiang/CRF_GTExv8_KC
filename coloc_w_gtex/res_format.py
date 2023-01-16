import pandas as pd
import os
import sys

def enrichment_data():
    # enrichment level of colocalization in different qtl celltypes
    dfs = [pd.read_csv(f'{RESFOLDER}/{QTLTYPE}.{FMTYPE}.{tsu}.enloc.enrich.out', sep='\s+', comment='#') for tsu in TISSUELIST]
    enrichlevel = [df.iloc[[1],[1]].values[0][0] for df in dfs]
    stderr = [df.iloc[[1],[2]].values[0][0] for df in dfs]
    data = pd.DataFrame({
        'Tissue': TISSUELIST, 
        'Enrichment Estimate (log-odds ratio, with shrinkage)': enrichlevel,
        'stderr': stderr
        })
    data.to_csv('%s/alltsu.%s.%s.enloc.enrich.out.tsv.gz'%(OUTPUTFOLDER, QTLTYPE,FMTYPE), sep='\t' , index=False, compression='gzip')

def findinfo(t, info):
    loc = t.find(info)
    t=t[loc:]
    t=t.split(';')[0]
    loc = t.find(' ')
    t=t[loc+1:]
    t=t.rstrip('"').lstrip('"')
    return t

def add_genesymbol():
    rcp_cs_dfs = [pd.read_csv(f'{RESFOLDER}/{QTLTYPE}.{FMTYPE}.{tsu}.enloc.sig.out', sep='\s+') for tsu in TISSUELIST]
    if QTLTYPE == 'sqtl':
        INTRONREF = pd.read_csv(INTRONREFPATH,sep='\t')
    for i in range(len(rcp_cs_dfs)):
        df = rcp_cs_dfs[i].copy()
        df['tsu'] = TISSUELIST[i]
        if QTLTYPE == 'eqtl':
            df['Geneid'] = df['Signal'].apply(lambda x: '.'.join(x.split(':')[0].split('_')))
            df = df.merge(GENEREF, how='left')
        elif QTLTYPE == 'sqtl':
            # res intron format: intron_10_103585420_103588685
            # ref intron format: chr22:16601407:16602064
            df['intron'] = df['Signal'].apply(lambda x: 'chr'+':'.join(x.split(':')[0].split('_')[1:]))
            df = df.merge(INTRONREF, how = 'left')
            df = df.groupby(['Signal', 'Num_SNP', 'CPIP_qtl', 'CPIP_gwas_marginal','CPIP_gwas_qtl_prior', 'RCP', 'tsu', 'intron','clu']).agg(list)
            df['Geneid'] = df['Geneid'].apply(lambda x: ','.join(x))
            df['Symbol'] = df['Symbol'].apply(lambda x: ','.join(x))
            df['gene_type'] = df['gene_type'].apply(lambda x: ','.join(x))
            df = pd.DataFrame(df.to_records())
        rcp_cs_dfs[i] = df.copy()
    rcp_cs = pd.concat(rcp_cs_dfs)
    # if QTLTYPE=='sqtl':
    #     rcp_cs = pd.DataFrame(rcp_cs.to_records())
    #     print(df['intron'])
    rcp_cs.to_csv('%s/alltsu.%s.%s.enloc.sig.out.tsv.gz'%(OUTPUTFOLDER, QTLTYPE,FMTYPE), sep='\t', index=False, compression='gzip')

    rcp_snp_dfs = [pd.read_csv(f'{RESFOLDER}/{QTLTYPE}.{FMTYPE}.{tsu}.enloc.snp.out', sep='\s+') for tsu in TISSUELIST]
    for i in range(len(rcp_snp_dfs)):
        df = rcp_snp_dfs[i].copy()
        df['tsu'] = TISSUELIST[i]
        if QTLTYPE == 'sqtl':
            df = df.merge(rcp_cs[['RCP','Signal','tsu','intron','clu','Geneid','Symbol','gene_type']])
        elif QTLTYPE == 'eqtl':
            df = df.merge(rcp_cs[['RCP','Signal','tsu','Symbol']])
        rcp_snp_dfs[i] = df.copy()
    rcp_snp = pd.concat(rcp_snp_dfs)
    rcp_snp.to_csv('%s/alltsu.%s.%s.enloc.snp.out.tsv.gz'%(OUTPUTFOLDER,QTLTYPE,FMTYPE), sep='\t', index=False, compression='gzip')
    return[rcp_cs, rcp_snp]

if __name__ == "__main__" :
    QTLTYPE=sys.argv[1]
    FMTYPE=sys.argv[2]
    RESFOLDER = sys.argv[3]
    OUTPUTFOLDER = sys.argv[4]
    GENEREFPATH = sys.argv[5]
    TISSUELISTPATH = sys.argv[6]
    INTRONREFPATH = sys.argv[7]
    # read generef
    GENEREF = pd.read_csv(GENEREFPATH,sep='\t',comment='#',header=None)
    GENEREF=pd.DataFrame(GENEREF[8])
    GENEREF['Geneid']=GENEREF[8].apply(lambda x: findinfo(x, 'gene_id'))
    GENEREF['Symbol']=GENEREF[8].apply(lambda x: findinfo(x, 'gene_name'))
    del GENEREF[8]
    GENEREF=GENEREF.drop_duplicates(keep='first')
    # read tissueslist
    TISSUELIST = list(pd.read_csv(TISSUELISTPATH, header = None)[0])
    # make output folder
    os.system(f'mkdir -p {OUTPUTFOLDER}')
    # perfrom analysis
    rcp_cs, rcp_snp = add_genesymbol()
    enrichment_data()