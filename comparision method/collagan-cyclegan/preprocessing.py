import matplotlib.pyplot as plt
import nibabel as nib
from PIL import Image
from tqdm import tqdm
import pandas as pd
import os
import numpy as np

# os.listdir('./Dataset/MICCAI_BraTS_2019_Data_Training/LGG/')

dataflie='dataverse_files'

types = ['breast_metastasis/','glioma_HG/','lung_metastasis/']
paths = []
path = f'./Dataset/{dataflie}/'
for i in types:
  p = os.listdir(path+i)
  paths+=[path+i+x+'/'+x for x in p]

# paths.sort(key=lambda x: int(x.split('-')[-1].split('.')[0]))

im_types = ['-t1.nii.gz','-t1ce.nii.gz','-t2.nii.gz','-flair.nii.gz']
file_map = {'-t1.nii.gz':'T1/','-t1ce.nii.gz':'T2/','-t2.nii.gz':'T1CE/','-flair.nii.gz':'FLAIR/'}
# slices = [20,55,100,110]

first_case_path = paths[0] + im_types[0]  # 例如�?'./Dataset/MICCAI_BraTS_2019_Data_Training/LGG/case01/case01_t1.nii'
im = nib.load(first_case_path)
print(im.get_fdata().shape)


from random import shuffle
ids = list(range(len(paths)))
shuffle(ids)
valid_ids = list(range(len(paths)))
train_ids = ids[:5]
# valid_ids = ids[:len(ids)//5]
# train_ids = ids[len(ids)//5:]
save_path = f'./Dataset/Preprocessed/{dataflie}/'

paths = np.array(paths)

for im_type in im_types:
  for path in tqdm(paths[train_ids]):
    img = nib.load(path+im_type)
    slices=img.get_fdata().shape[2]
    for s in range(slices):
      Image.fromarray(img.get_fdata()[:,:,s]).save(save_path+'Train/'+file_map[im_type]+path[path.rfind('/')+1:]+'_'+str(s)+'.tiff')

for im_type in im_types:
  aaa = []
  for path in tqdm(paths[valid_ids]):
    img = nib.load(path+im_type)
    slices = img.get_fdata().shape[2]
    aaa.append(slices)
    for s in range(slices):
      Image.fromarray(img.get_fdata()[:,:,s]).convert('L').save(save_path+'/Valid/'+file_map[im_type]+path[path.rfind('/')+1:]+'_'+str(s)+'.png')
  df = pd.DataFrame(aaa)
  # 将DataFrame写入CSV文件
  df.to_csv(f'./Dataset/Preprocessed/{dataflie}_slice_data.csv', header=False, index=False)
len(os.listdir(save_path+'Valid/T1')),len(os.listdir(save_path+'Valid/T1CE')),len(os.listdir(save_path+'Valid/T2')),len(os.listdir(save_path+'Valid/FLAIR'))

len(os.listdir(save_path+'Train/T1')),len(os.listdir(save_path+'Train/T1CE')),len(os.listdir(save_path+'Train/T2')),len(os.listdir(save_path+'Train/FLAIR'))


# for path in tqdm(paths[valid_ids]):
#   img = nib.load(path + im_types[3])
#   print(path)
#   print(img.shape)