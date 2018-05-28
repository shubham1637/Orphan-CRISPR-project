
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd


# In[33]:


data_path = 'Desktop/Ira/lefties.csv'
data = pd.read_csv(data_path, sep =',', names = ['ID', 'V1', 'GenRef','start','end','length','repeat','repeatlength','V8','N_spacers','Known','Bacteria','V12','V13','V14'])
#data = data.iloc[1:,]
data = data.loc[:, "V1":]
data


# In[67]:


contig_id = []
start =[]
end = []
for i in range(data.shape[0]-1):
    contig_id.append(data.GenRef[i+1])
    start.append(int(data.start[i+1]))
    end.append(int(data.end[i+1]))
start[0:9]


# #CodeStyle

# In[65]:


import subprocess
for i in n:
    
#subprocess.call(
        #"blastdbcmd -db " + Database + " -entry_batch " + ClusterGIListFileName + " > " + ClusterFASTA, shell = True)


