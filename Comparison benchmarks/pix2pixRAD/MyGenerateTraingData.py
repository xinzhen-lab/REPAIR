from __future__ import absolute_import, division, print_function, unicode_literals
from __future__ import print_function
import tensorflow as tf

import numpy as np
import os
import time

from matplotlib import pyplot as plt
from IPython import display

from PIL import Image




import nibabel as nib
from skimage.transform import resize as rsz
import imageio
from sklearn.metrics import mean_squared_error
from skimage.metrics import structural_similarity as ssim
import pandas as pd

#### https://github.com/giemmecci/pix2pixRAD/tree/master

################################ added by Xin #############################################

# The following code is to generate the PNG images for training

train_data_path = r'C:\Users\THINKPAD\Desktop\pix2pixRAD\mydataset\myBRATS\train'

image_rows = 256
image_cols = 256

patients = os.listdir(train_data_path)
print(patients)

for patient in patients:
    print(f'WORKING ON PATIENTS: {patient}')

    patient_path = os.path.join(train_data_path, patient)
    scans = os.listdir(patient_path)

    ## Create the folders to store pngs and nifti
    PNGs_dir = 't2-to-t1ce-pngs'

    png_dir = os.path.join(train_data_path, PNGs_dir)

    # os.mkdir(png_dir)

    png_path = os.path.join(train_data_path, png_dir)


    # Create the numpy arrays to generate pngs
    for scan in scans:
        if 't2.nii' in scan:
            t2path = os.path.join(patient_path, scan)
            print(t2path)

            t2_img = nib.load(t2path)
            t2_imgdata = t2_img.get_fdata()
            t2_imgdata = t2_imgdata.astype(float)
            t2_affine = t2_img.affine

            resized_t2_imgdata = rsz(t2_imgdata, (256, 256), order=1, preserve_range=True)

            resized_t2_imgdata -= resized_t2_imgdata.min()
            resized_t2_imgdata *= 1. / np.max(resized_t2_imgdata)
            resized_t2_imgdata *= 255
            resized_t2_imgdata = np.flipud(resized_t2_imgdata)

        if 't1ce.nii' in scan:
            flairpath = os.path.join(patient_path, scan)
            print(flairpath)

            flair_img = nib.load(flairpath)
            flair_imgdata = flair_img.get_fdata()
            flair_imgdata = flair_imgdata.astype(float)
            # the FLAIR affine will be used to build the nifti from gan's output
            flair_affine = flair_img.affine

            resized_flair_imgdata = rsz(flair_imgdata, (256, 256), order=1, preserve_range=True)

            resized_flair_imgdata -= resized_flair_imgdata.min()
            resized_flair_imgdata *= 1. / np.max(resized_flair_imgdata)
            resized_flair_imgdata *= 255
            resized_flair_imgdata = np.flipud(resized_flair_imgdata)

    ## use the numpy arrays to generate pngs
    for i in range(0, resized_t2_imgdata.shape[2]):
        slice_t2 = resized_t2_imgdata[:, :, i]
        slice_flair = resized_flair_imgdata[:, :, i]

        slice_t2 = np.rot90(slice_t2, 3)
        slice_flair = np.rot90(slice_flair, 3)

        ## generates the input for pix2pix stack: INPUT -> TARGET
        fusion = np.hstack((slice_t2, slice_flair)).astype('uint8')
        fn = str(patient) + '_' + str(i).zfill(3) + '.png'
        outdir = png_path
        outname = os.path.join(outdir, fn)
        imageio.imwrite(outname, fusion, compress_level=0)


    # # Predictions input
    # test_prediction = tf.data.Dataset.list_files(outdir + '/*.png', shuffle=False)
    # test_prediction = test_prediction.map(load_image_test)
    # test_prediction = test_prediction.batch(BATCH_SIZE)


############################################################################################

print('Generation of training data done!')