import matplotlib.pyplot as plt
# import nibabel as nib
from PIL import Image
# from tqdm import tqdm_notebook
from tqdm import notebook
import os
import numpy as np
import pytorch_ssim
from cyclegan_model import GeneratorG,GeneratorF,Discriminator
import nibabel as nib
import torch
from torch import nn, optim
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm_notebook as tqdm
from time import time
import torch.nn.parallel
import torch.optim as optim
import torch.utils.data
import torchvision
# from torch.autograd import Variable
import random
from torch.nn.utils import spectral_norm
# from scipy.stats import truncnorm
import torch as th
from torchvision import transforms
import re

# torch.backends.cudnn.enabled = False

def LSLoss(y, yhat):
    return torch.mean((y - yhat) ** 2).to(device)


class BRATSDataSet(Dataset):

    def __init__(self, path, transform=transforms.ToTensor(), types=['T1/', 'T2/', 'T1CE/', 'FLAIR/']):
        self.path = path
        self.type = types
        self.names = os.listdir(path + 'T1/')
        self.transforms = transform
        if len(self.names) == 0:
            raise RuntimeError("Found 0 files in {}".format(path))

    def __getitem__(self, idx):

        if torch.is_tensor(idx):
            idx = idx.tolist()

        t1 = self.transforms(Image.open(self.path + self.type[0] + self.names[idx]))
        t2 = self.transforms(Image.open(self.path + self.type[1] + self.names[idx]))
        t1ce = self.transforms(Image.open(self.path + self.type[2] + self.names[idx]))
        flair = self.transforms(Image.open(self.path + self.type[3] + self.names[idx]))

        return {'T1': t1, 'T2': t2, 'T1CE': t1ce, 'FLAIR': flair}

    def __len__(self):
        return len(self.names)


dataflie = 'MICCAI_BraTS_2019_Data_Training'
# train_path = f'/home/wcx/new deep learn/Brain-Image-Augmentation-using-GAN-in-pytorch-master/Dataset/Preprocessed/{dataflie}/Train/'
valid_path = f'/home/wcx/new deep learn/Brain-Image-Augmentation-using-GAN-in-pytorch-master/Dataset/Preprocessed/{dataflie}/Valid/'
save_path = f'/home/wcx/new deep learn/Brain-Image-Augmentation-using-GAN-in-pytorch-master/CycleGAN/3D/{dataflie}/'


bs = 30


# train_dataset = BRATSDataSet(train_path)
valid_dataset = BRATSDataSet(valid_path)
file_names = valid_dataset.names
data_loader = torch.utils.data.DataLoader(valid_dataset, shuffle=True, batch_size=bs)


device = "cuda"

retrain = True
dataflie1 = 'MICCAI_BraTS_2019_Data_Training'
if retrain:
    G=torch.load(f'./model_save/FLAIR_{dataflie}_G_final.pt')
    # G=torch.load(f'./model_save/T1_{dataflie1}_G_final.pt')
    # G=torch.load(f'./model_save/T2_{dataflie}_G_final.pt')
    # G=torch.load(f'./model_save/T1CE_{dataflie}_G_final.pt')
else:
    Gen = GeneratorG(5, True).to(device)
    # Dis = Discriminator(True).to(device)





for i, imgs in notebook.tqdm(enumerate(data_loader), total=len(data_loader)):
    print('aaa')

    start_idx = i * bs


    # 获取当前批次的文件名
    # batch_file_names = file_names[start_idx].replace("_55.png", "")
    batch_file_names = re.sub(r'_\d+\.png$', '', file_names[start_idx])

    t1 = imgs['T1'].to(device)
    t2 = imgs['T2'].to(device)
    t1ce = imgs['T1CE'].to(device)
    # flair = imgs['FLAIR'].to(device)

    X = torch.cat([t1, t2, t1ce], dim=1).to(device)
    yhat = G(X)

    # t1_recon_np = yhat.cpu().detach().numpy().squeeze()  # 去除多余的维度
    # t2_recon_np = yhat.cpu().detach().numpy().squeeze()
    # t1ce_recon_np = yhat.cpu().detach().numpy().squeeze()
    flair_recon_np = yhat.cpu().detach().numpy().squeeze()

    # 调整形状为 {240, 240, 35}
    # t1_recon_np = t1_recon_np.transpose((1, 2, 0))  # 从 {35, 240, 240} 到 {240, 240, 35}
    # t2_recon_np = t2_recon_np.transpose((1, 2, 0))
    # t1ce_recon_np = t1ce_recon_np.transpose((1, 2, 0))
    flair_recon_np = flair_recon_np.transpose((1, 2, 0))

    # 保存为 NIfTI 格式
      # 指定保存路径



    folder_path = os.path.join(save_path, batch_file_names)

    # 判断文件夹是否存在，如果不存在则创建
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    else:
        print(f"'{batch_file_names}' 已存在于路径 {folder_path}")

    # nib.save(nib.Nifti1Image(t1_recon_np, np.eye(4)), os.path.join(folder_path, f'{batch_file_names}_t1.nii'))
    # nib.save(nib.Nifti1Image(t2_recon_np, np.eye(4)), os.path.join(folder_path, f'{batch_file_names}_t2.nii'))
    # nib.save(nib.Nifti1Image(t1ce_recon_np, np.eye(4)), os.path.join(folder_path, f'{batch_file_names}_t1ce.nii'))
    nib.save(nib.Nifti1Image(flair_recon_np, np.eye(4)), os.path.join(folder_path, f'{batch_file_names}_flair.nii'))




