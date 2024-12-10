# import pandas as pd
# import numpy as np
# df=pd.read_pickle('F:/incomplete modle/11_CDIMC-net/IJCAI_final_code_handwritten_release1/IJCAI_final_code_handwritten_tochen/data/handwritten-5view_2_0.3_aelr_0.01_aeproches_500_pretrained_model.pkl')
# print(df)

import pickle
fr=open('F:/incomplete modle/11_CDIMC-net/IJCAI_final_code_handwritten_release/IJCAI_final_code_handwritten_tochen/data/handwritten-5view_2_0.3_aelr_0.01_aeproches_500_pretrained_model.pkl','rb')
inf = pickle.load(fr)
doc = open('1.txt', 'a')
print(inf, file=doc)

# import torch
# torch.load('handwritten-5view_2_0.3_aelr_0.01_aeproches_500_pretrained_model.pkl')