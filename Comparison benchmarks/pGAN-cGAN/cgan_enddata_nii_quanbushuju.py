import os
import numpy as np
import nibabel as nib
from PIL import Image
import pandas as pd

def natural_sort_key(file_name):
    # 从文件名提取数字部分并转换为整数
    return int(''.join(filter(str.isdigit, file_name)))

def process_images(path, method, image_type,panduan):
    # 文件夹路径
    input_folder = f'/home/wcx/new deep learn/pGAN-cGAN-master/result/{path}/{method}/test_latest/t2A_flairB_images'
    output_folder = f'/home/wcx/new deep learn/pGAN-cGAN-master/3D/{path}/{method}/{image_type}'
    csv_path = f'/home/wcx/new deep learn/pGAN-cGAN-master/datasets/{path}_name_data.csv'
    # 确保输出文件夹存在
    os.makedirs(output_folder, exist_ok=True)

    # 读取 CSV 文件
    brats_names = pd.read_csv(csv_path, header=None)[0].tolist()  # 读取第一列并转换为列表


    # 获取所有图片文件名，并按顺序排序

    image_files = sorted(os.listdir(input_folder), key=natural_sort_key)
    image_files = [name for name in image_files if panduan in name]

    # 读取 batch_sizes 从 CSV 文件
    batch_sizes_file = pd.read_csv(f'/home/wcx/new deep learn/pGAN-cGAN-master/datasets/{path}_slice_data.csv', header=None)
    batch_sizes = batch_sizes_file[0].astype(int).tolist()

    # 读取图片并按批次处理
    num_images = len(image_files)
    current_index = 0
    name=0
    for batch_size in batch_sizes:
        # 确保当前批次的图片数量不超过总图片数量
        current_batch_size = min(batch_size, num_images - current_index)

        # 获取当前批次的文件名
        batch_files = image_files[current_index:current_index + current_batch_size]
        current_index += current_batch_size

        # 如果当前批次没有图片，跳过处理
        if not batch_files:
            continue

        # 初始化一个列表用于存储当前批次的图像数据
        batch_images = []

        for file in batch_files:
            file_path = os.path.join(input_folder, file)

            # 使用 PIL 读取图像
            image = Image.open(file_path).convert('L')  # 转换为灰度图像

            # 将图像转换为 NumPy 数组
            image_array = np.array(image)

            # 确保图像维度为 [H, W]，增加通道维度
            image_array = image_array[:, :, np.newaxis]  # 增加通道维度

            batch_images.append(image_array)

        # 将当前批次的图像堆叠成一个 3D 数组 [H, W, D]，D 为批次大小
        batch_array = np.concatenate(batch_images, axis=2)  # 沿着深度方向堆叠

        # 创建 NIfTI 图像对象
        nii_image = nib.Nifti1Image(batch_array, affine=np.eye(4))

        # 生成保存文件名
        output_filename = f'{brats_names[name]}.nii'  # 使用 CSV 第一列的名称
        output_path = os.path.join(output_folder, output_filename)
        name+=1
        # 保存 NIfTI 文件
        nib.save(nii_image, output_path)

        print(f'Saved {output_path}')
    print('1')

# 调用函数，处理不同类型的图像
path = 'dataverse'
method = 'cGAN_run'
panduan="fake_B"
# process_images(path, method, 't1ce',panduan)
# process_images(path, method, 't2',panduan)  # 处理 t1ce 类型的图像
# process_images(path, method, 't1',panduan)
process_images(path, method, 'flair',panduan)
