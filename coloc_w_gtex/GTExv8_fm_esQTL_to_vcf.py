import os
import pandas as pd
fo = open("tissuelist", "r")
for line in fo.readlines():
  line = line.strip()
  # unzip the data downloaded from https://zenodo.org/record/3517189#.Y-0UwcfP2Ul
  os.system('gunzip -c %s.variants_pip.txt.gz > ./hakyim_original_data/%s.variants_pip.txt'%(line,line))
  os.system('gunzip -c %s.clusters.txt.gz > ./hakyim_original_data/%s.clusters.txt'%(line,line))
  os.system('mkdir ./fenloc_gtexv8/%s'%line)
  # read and reform the data
  pip = pd.read_csv('./hakyim_original_data/%s.variants_pip.txt'%line,sep='\t')
  cluster = pd.read_csv('./hakyim_original_data/%s.clusters.txt'%line,sep='\t')
  pip['rank']=pip['rank'].apply(lambda x: '(('+str(x)+'))')
  cluster['cluster']=cluster['cluster'].apply(lambda x: '{'+str(x)+'}')
  pip.set_index(['gene'],inplace=True)
  cluster.set_index(['gene'],inplace=True)
  genes = list(cluster.index.drop_duplicates())
  for gene in genes:
    output=pip.loc[gene]
    clusterid = cluster.loc[gene,'cluster']
    
    if type(output) == pd.Series:
      output=pd.DataFrame(output).transpose()
    
    if type(clusterid) == str:
        clusterid=int(clusterid.rstrip('}').lstrip('{'))
        output=output[output['cluster_id']==clusterid]
    else:
        clusterid=clusterid.apply(lambda x: x.rstrip('}').lstrip('{')).astype(int)
        output=output[output['cluster_id'].isin(clusterid)]
    gene1=gene.replace('.','_')
    if type(output) == pd.DataFrame:
        output.to_csv('%s.pipcluster'%gene1,index=False,sep='\t')
    elif type(output) == pd.Series:
        pd.DataFrame(output).transpose().to_csv('%s.pipcluster'%gene1,index=False,sep='\t')
    else:
        print('error')
    output2=cluster.loc[gene]
    if type(output2) == pd.DataFrame:
        output2.to_csv('%s.cluster'%gene1,index=False,sep='\t')
    elif type(output2) == pd.Series:
        pd.DataFrame(output2).transpose().to_csv('%s.cluster'%gene1,index=False,sep='\t')
    else:
        print('error')
    os.system('cat %s.cluster >> %s.pipcluster'%(gene1,gene1))
    os.system('rm %s.cluster'%gene1)
    os.system('mv %s.pipcluster ./fenloc_gtexv8/%s/'%(gene1,line))
  
    
  #make vcf file
  pip['#CHROM']=pip['variant_id'].apply(lambda x:x.split('_')[0])
  pip['POS']=pip['variant_id'].apply(lambda x:x.split('_')[1])
  pip['REF']=pip['variant_id'].apply(lambda x:x.split('_')[2])
  pip['ALT']=pip['variant_id'].apply(lambda x:x.split('_')[3])
  pip=pip[['#CHROM','POS','variant_id','REF','ALT']]
  pip.to_csv('./fenloc_gtexv8/%s/snps.vcf'%(line),sep='\t',index=False)
  os.system('gzip ./fenloc_gtexv8/%s/snps.vcf'%(line))
  os.system('./summarize_dap2enloc.pl -dir ./fenloc_gtexv8/%s/ -vcf ./fenloc_gtexv8/%s/snps.vcf.gz -tissue %s | gzip - > ./fenloc_gtexv8/%s/fastenloc.eqtl.annotation.%s.vcf.gz'%(line,line,line,line,line))
  
  # drop duplicates 
  df = pd.read_csv('./fenloc_gtexv8/%s/fastenloc.eqtl.annotation.%s.vcf.gz'%(line,line),sep='\t',header=None)
  df = df.drop_duplicates(keep='first')
  df.to_csv('./fenloc_gtexv8/fastenloc.eqtl.annotation.%s.vcf.gz'%(line),sep='\t',index=False,header=None,compression='gzip')
  print('%s done'%(line))