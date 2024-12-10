import os
import numpy as np
import nibabel as nib
from PIL import Image
import pandas as pd

def natural_sort_key(file_name):

    return int(''.join(filter(str.isdigit, file_name)))

def process_images(path, method, image_type):

    input_folder = f'/home/wcx/new deep learn/pGAN-cGAN-master/result/{path}/{method}/test_latest/{image_type}_images'
    output_folder = f'/home/wcx/new deep learn/pGAN-cGAN-master/3D/{path}/{method}/{image_type}'
    csv_path = f'/home/wcx/new deep learn/pGAN-cGAN-master/datasets/{path}_name_data.csv'

    os.makedirs(output_folder, exist_ok=True)


    brats_names = pd.read_csv(csv_path, header=None)[0].tolist()


    image_files = sorted(os.listdir(input_folder), key=natural_sort_key)


    batch_sizes_file = pd.read_csv(f'/home/wcx/new deep learn/pGAN-cGAN-master/datasets/{path}_slice_data.csv', header=None)
    batch_sizes = batch_sizes_file[0].astype(int).tolist()


    num_images = len(image_files)
    current_index = 0
    name=0
    for batch_size in batch_sizes:

        current_batch_size = min(batch_size, num_images - current_index)


        batch_files = image_files[current_index:current_index + current_batch_size]
        current_index += current_batch_size

        if not batch_files:
            continue


        batch_images = []

        for file in batch_files:
            file_path = os.path.join(input_folder, file)


            image = Image.open(file_path).convert('L')


            image_array = np.array(image)


            image_array = image_array[:, :, np.newaxis]

            batch_images.append(image_array)


        batch_array = np.concatenate(batch_images, axis=2)


        nii_image = nib.Nifti1Image(batch_array, affine=np.eye(4))

        output_filename = f'{brats_names[name]}.nii'
        output_path = os.path.join(output_folder, output_filename)
        name+=1

        nib.save(nii_image, output_path)

        print(f'Saved {output_path}')
    print('1')


path = 'dataverse'
method = 'pGAN_run'
process_images(path, method, 't1ce')
process_images(path, method, 't2')
process_images(path, method, 't1')
process_images(path, method, 'flair')
