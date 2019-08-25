import pandas as pd

fileN= 'Data/web-Google.txt'

DF = pd.read_csv(fileN,sep='\t')

DF = DF.sort_values('FromNodeId')

print(max( DF['FromNodeId'].max(),DF['ToNodeId'].max()))

DF.to_csv('new',sep='\t',index=False)

