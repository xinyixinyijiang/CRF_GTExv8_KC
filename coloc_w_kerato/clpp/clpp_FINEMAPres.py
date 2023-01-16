import pandas as pd 
crf =  pd.read_csv('../../data/crf/FINEMAP_credible_sets.tsv',sep='\t')
k = pd.read_csv('./FINEMAP/summary/keratoconus_finemaping_res.tsv',sep='\t')
crf = crf[['GenomicLocus','chr','start','end','in_credible_set_of','SNP','posterior_probability']]
crf.columns = ['crfloci','chr','start','end','in_credible_set_of','SNP','prob']
k.columns = ['SNP','prob','cs','crfloci','chr','start','end']
k = k[['crfloci','chr','start','end','cs','SNP','prob']]

loci = pd.read_csv('../crf_kerato_overlaploci.tsv', sep='\t', header = None)
crfloci = list(loci[3].apply(lambda x: int(x.split('crf_locus')[1])))
loci['loci'] = 'kerato'+loci[7].apply(lambda x: x.split('kerato_locus')[1])+'_'+'crf'+loci[3].apply(lambda x: x.split('crf_locus')[1])
loci['crfloci'] = crfloci

crf = crf[crf.crfloci.isin(loci.crfloci)]
k = k[k.crfloci.isin(loci.crfloci)]
loci.rename(columns = {0: 'chr',1: 'start', 2: 'end'}, inplace=True)
crf = crf.merge(loci[['chr', 'start', 'end', 'loci']])
crf.rename(columns = {'loci': 'overlaping_loci'}, inplace=True)
k = k.merge(loci[['chr', 'start', 'end', 'loci']])
k.rename(columns = {'loci': 'overlaping_loci'}, inplace=True)

del k['crfloci']
del crf['crfloci']
k.columns = ['chr', 'start', 'end', 'kerato_cs', 'rsid', 'kerato_pip', 'overlapping_loci']
crf.columns = ['chr', 'start', 'end', 'crf_cs', 'rsid', 'crf_pip', 'overlapping_loci']
for i in loci.loci:
    leadsnps = list(set(crf[crf.overlapping_loci == i]['crf_cs']))
    for j in range(len(leadsnps)):
        crf.crf_cs.replace({leadsnps[j]:j+1}, inplace=True)

res = crf.merge(k,how='outer')

#### NULL reason #### 
#### check reason why keratoconus is na (whether "not selected by FINEMAP" or "not in keratoconus summary statistics")
res['crf_cs'].fillna('not selected by FINEMAP', inplace=True)
res['crf_pip'].fillna('not selected by FINEMAP', inplace=True)

ids = res[res.kerato_cs.isna()].index
keratosumstats = pd.read_csv('../../data/keratoconus/42003_2021_1784_MOESM4_ESM.txt', sep='\t')
keratosumstats = keratosumstats['variant_id']
instats = res.loc[ids][res.loc[ids,'rsid'].isin(keratosumstats)].index
res.loc[instats,'kerato_cs'] = 'not selected by FINEMAP'
res.loc[instats,'kerato_pip'] = 'not selected by FINEMAP'
res['kerato_cs'].fillna('not in keratoconus GWAS summary stats', inplace=True)
res['kerato_pip'].fillna('not in keratoconus GWAS summary stats', inplace=True)
res.to_csv('FINEMAP_res.tsv', sep='\t',index=False)


#### calculate clpp #### 
crf['CLPP'] = 'NA'
inid = res[
    (~res.crf_cs.isin(['not selected by FINEMAP']))&
    (~res.kerato_cs.isin(['not selected by FINEMAP','not in keratoconus GWAS summary stats']))
    ].index
for i in loci.loci:
    df = res.loc[inid]
    df = df[df['overlapping_loci']==i]
    css = list(set(df['crf_cs']))
    for j in css:
        df1 = df[df['crf_cs']==j]
        df1['1-kpipxcrfpip'] = 1 - df1['crf_pip'] * df1['kerato_pip']
        clpp = 1
        for k in df1.index:
            clpp = clpp*df1['1-kpipxcrfpip'][k]
        clpp = 1-clpp
        res.loc[res[(res['overlapping_loci']==i)&(res['crf_cs']==j)].index,'CLPP'] = clpp

res.CLPP.fillna('no overlap between crf cs and kerato cs', inplace=True)
res.to_csv('clpp.tsv', sep='\t', index=False)