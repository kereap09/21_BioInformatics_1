#!/usr/bin/env python
# coding: utf-8

# ## Transcript count table 마련하기
# 
# 제공된 alignment와 annotation 파일들을 이용해서 transcript별 read count를 구해봅니다. 실제 연구에서는 multi-mapping 등을 고려해야 하지만, 여기서는 단순화해서 모두 무시합니다. 데이터가 있는 곳으로 이동해서 작업합시다.

# In[1]:


import pandas as pd
from matplotlib import pyplot as plt
import numpy as np


# In[2]:


cnts = pd.read_csv('read-counts.txt', sep='\t', comment='#', index_col=0)
cnts.head()


# In[3]:


cnts['clip_enrichment'] = cnts['CLIP-35L33G.bam'] / cnts['RNA-control.bam']
cnts['rden_change'] = (cnts['RPF-siLin28a.bam'] / cnts['RNA-siLin28a.bam']) / (cnts['RPF-siLuc.bam'] / cnts['RNA-siLuc.bam'])
cnts.head()


# In[5]:


fig, ax = plt.subplots(1, 1, figsize=(5, 5))
plt.grid()
ax.scatter(np.log2(cnts['clip_enrichment']),
           np.log2(cnts['rden_change']),
           s=0.2,
           c = "black")
plt.xlabel("CLIP_enrichment(log2)")
plt.ylabel("rden_change(log2)")


# 이 그림에는 문제가 많이 있습니다. 논문의 그림처럼 한 번 만들어 봅시다~

# ## Protein localization 반영하기
# 
# 논문 그림에서는 protein localization을 반영해서 색으로 나타냈습니다. 이것을 한 번 구현해 봅시다. 각 gene product의 localization을 따로 구하려면 좀 귀찮기 때문에 미리 모아 두었습니다. UniProt에서 매핑된 데이터를 토대로 아주 명확한 것만 일부 정리한 것입니다. 이것도 마찬가지로 처음엔 에러가 나기도 하니까, 다시 실행해보면 거의 잘 됩니다~

# In[4]:


import ssl
ssl._create_default_https_context = ssl._create_unverified_context
mouselocal = pd.read_csv('https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt', sep='\t')
mouselocal.head()


# 이제 이 localization 데이터와 위에서 만든 scatter를 결합해서 논문 그림과 비슷하게 만들어 봅시다.

# In[ ]:





# ## 이제 시작

# In[ ]:





# In[115]:


integ_df = cnts.loc[:,['clip_enrichment','rden_change']]
integ_df["gene_id"] = [x.split(".")[0] for x in integ_df.index]


# In[30]:


integ_df


# In[53]:


# integ_df = integ_df.merge(mouselocal, on="gene_id", how="right")
# integ_df


# In[116]:


integ_df = integ_df.merge(mouselocal, on="gene_id", how="right")
integ_df


# In[98]:


from scipy.stats import pearsonr


# In[117]:


integ_df.replace([np.inf, -np.inf], np.nan, inplace=True)


# In[ ]:





# In[119]:


integ_df_f = integ_df.dropna()


# In[120]:


integ_df_f


# In[121]:


corr, _ = pearsonr(integ_df_f["clip_enrichment"], integ_df_f["rden_change"])
corr


# In[135]:


fig, ax = plt.subplots(1, 1, figsize=(5, 5))
plt.grid()



ax.scatter(np.log2(integ_df_f.loc[integ_df["type"]=="nucleus",'clip_enrichment']),
           np.log2(integ_df_f.loc[integ_df["type"]=="nucleus",'rden_change']),
           s=0.8,
           label = "Nucleus",
           c = "blue")
ax.scatter(np.log2(integ_df_f.loc[integ_df["type"]=="integral membrane",'clip_enrichment']),
           np.log2(integ_df_f.loc[integ_df["type"]=="integral membrane",'rden_change']),
           s=0.8,
           label = "integral membrane",
           c = "red")
ax.scatter(np.log2(integ_df_f.loc[integ_df["type"]=="cytoplasm",'clip_enrichment']),
           np.log2(integ_df_f.loc[integ_df["type"]=="cytoplasm",'rden_change']),
           s=0.8,
           label = "cytoplasm",
           c = "green")
# for color in ['nucleus:blue', 'integral membrane:red', 'cytoplasm:green']:
#     ax.scatter(np.log2(integ_df['clip_enrichment']),
#                np.log2(integ_df['rden_change']),
#                c=color, s=0.2, edgecolors='none')

# ax.scatter(np.log2(integ_df['clip_enrichment']),
#            np.log2(integ_df['rden_change']),
#            s=0.2,
#            c = integ_df['type'])

plt.legend()
plt.xlabel("LIN28A CLIP enrichment(log2)")
plt.ylabel("Ribosome density change \nupon LIN28a knockdown(log2)")
plt.text(6, -5,'r = ' + str(round(corr,4)),
     horizontalalignment='center',
     verticalalignment='center')


# In[ ]:




