import matplotlib.pyplot as plt
# import nibabel as nib
from PIL import Image
# from tqdm import tqdm_notebook
from tqdm import notebook
import os
import numpy as np
import pytorch_ssim
from collagan_model import Generator, Discriminator
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
        # self.names = sorted(self.names, key=natural_sort_key)
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
def natural_sort_key(file_name):
    # 从文件名提取数字部分并转换为整数
    return int(''.join(filter(str.isdigit, file_name)))

dataflie = 'MICCAI_BraTS_2019_Data_Training'
train_path = f'/home/wcx/new deep learn/Brain-Image-Augmentation-using-GAN-in-pytorch-master/Dataset/Preprocessed/{dataflie}/Train/'
valid_path = f'/home/wcx/new deep learn/Brain-Image-Augmentation-using-GAN-in-pytorch-master/Dataset/Preprocessed/{dataflie}/Valid/'
save_path = f'/home/wcx/new deep learn/Brain-Image-Augmentation-using-GAN-in-pytorch-master/CollaGAN/3D/{dataflie}/'


bs = 8
gen_lr = 1e-5
dis_lr = 1e-3
lambda_gen_ls = 0.5
lambda_l1_cyc = 10
lambda_l1 = 1
lambda_l2_cyc = 10
lambda_l2 = 0
lambda_ce_gen = 15
lambda_ssim = 1

# train_dataset = BRATSDataSet(train_path)
valid_dataset = BRATSDataSet(valid_path)

file_names = valid_dataset.names
data_loader = torch.utils.data.DataLoader(valid_dataset, shuffle=True, batch_size=bs)

a = next(iter(data_loader))
for i in a:
    a[i] = a[i][0].numpy().transpose(1, 2, 0).reshape(240, 240)

fig, axs = plt.subplots(1, 4, figsize=(10, 10))
for i, j in enumerate(a):
    axs[i].imshow(a[j], cmap='Greys_r')
    axs[i].set_title(j)

device = "cuda"

retrain = True
dataflie1='MICCAI_BraTS_2019_Data_Training'
if retrain:
    Gen = torch.load(f'./model_save/{dataflie1}_generator.pth').to(device)
    Dis = torch.load(f'./model_save/{dataflie1}_discriminator.pth').to(device)
else:
    Gen = Generator(5, True).to(device)
    Dis = Discriminator(True).to(device)

optimizer_g = optim.Adam(Gen.parameters(), lr=gen_lr, betas=(0.5, 0.999))
optimizer_d = optim.Adam(Dis.parameters(), lr=dis_lr, betas=(0.5, 0.999))

l2 = nn.MSELoss()
l1 = nn.L1Loss()
ssim_loss = pytorch_ssim.SSIM()
ce = nn.CrossEntropyLoss()

ssim_lambda = lambda x, y, r: -torch.log((1.0 + x) / 2.0) if y != r else torch.zeros_like(x, device=device)
l1_lambda = lambda x, y, r: x if y != r else torch.zeros_like(x, device=device)
l2_lambda = lambda x, y, r: x if y != r else torch.zeros_like(x, device=device)

real = torch.ones((32, 1), device=device)
fake = torch.zeros((32, 1), device=device)


for i, imgs in notebook.tqdm(enumerate(data_loader), total=len(data_loader)):
    print('aaa')

    start_idx = i * bs
    print(file_names[start_idx:start_idx+bs])

    # 获取当前批次的文件名
    # batch_file_names = file_names[start_idx].replace("_55.png", "")
    batch_file_names = re.sub(r'_\d+\.png$', '', file_names[start_idx])
    t1 = imgs['T1'].to(device)
    t2 = imgs['T2'].to(device)
    t1ce = imgs['T1CE'].to(device)
    flair = imgs['FLAIR'].to(device)
    # if a.shape[0]<bs:
    #   break
    # print(a.shape[0], b.shape, c.shape, d.shape, mask.shape)

    mask_0 = torch.zeros((bs, 4, 240, 240), device=device)
    mask_0[:, 0:1, :, :] = 1
    mask_1 = torch.zeros((bs, 4, 240, 240), device=device)
    mask_1[:, 1:2, :, :] = 1
    mask_2 = torch.zeros((bs, 4, 240, 240), device=device)
    mask_2[:, 2:3, :, :] = 1
    mask_3 = torch.zeros((bs, 4, 240, 240), device=device)
    mask_3[:, 3:4, :, :] = 1

    dummy = torch.zeros((bs, 1, 240, 240), device=device)

    t1_recon = Gen(torch.cat([dummy, t2, t1ce, flair, mask_0], dim=1))
    t2_recon = Gen(torch.cat([t1, dummy, t1ce, flair, mask_1], dim=1))
    t1ce_recon = Gen(torch.cat([t1, t2, dummy, flair, mask_2], dim=1))
    flair_recon = Gen(torch.cat([t1, t2, t1ce, dummy, mask_3], dim=1))

    t1_recon_np = t1_recon.cpu().detach().numpy().squeeze()  # 去除多余的维�?
    t2_recon_np = t2_recon.cpu().detach().numpy().squeeze()
    t1ce_recon_np = t1ce_recon.cpu().detach().numpy().squeeze()
    flair_recon_np = flair_recon.cpu().detach().numpy().squeeze()

    # 调整形状�? {240, 240, 35}
    t1_recon_np = t1_recon_np.transpose((1, 2, 0))  # �? {35, 240, 240} �? {240, 240, 35}
    t2_recon_np = t2_recon_np.transpose((1, 2, 0))
    t1ce_recon_np = t1ce_recon_np.transpose((1, 2, 0))
    flair_recon_np = flair_recon_np.transpose((1, 2, 0))

    # 保存�? NIfTI 格式
      # 指定保存路径



    folder_path = os.path.join(save_path, batch_file_names)

    # 判断文件夹是否存在，如果不存在则创建
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    else:
        print(f"{batch_file_names}' 已存在于路径 {folder_path}")

    nib.save(nib.Nifti1Image(t1_recon_np, np.eye(4)), os.path.join(folder_path, f'{batch_file_names}_t1.nii'))
    nib.save(nib.Nifti1Image(t2_recon_np, np.eye(4)), os.path.join(folder_path, f'{batch_file_names}_t2.nii'))
    nib.save(nib.Nifti1Image(t1ce_recon_np, np.eye(4)), os.path.join(folder_path, f'{batch_file_names}_t1ce.nii'))
    nib.save(nib.Nifti1Image(flair_recon_np, np.eye(4)), os.path.join(folder_path, f'{batch_file_names}_flair.nii'))




