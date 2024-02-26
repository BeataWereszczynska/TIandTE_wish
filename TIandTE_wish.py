# -*- coding: utf-8 -*-
"""
Calculating theoretical MRI images with both TI (T1-weighting) and TE (T2-weighting)
of choise, from separate T1-weighted and T2-weighted sets of images.
For Agilent SEMS (with IR) and MEMS .fid data acquired for the same slices and matrice size.

@author: Beata Wereszczyńska
"""

# .......... User defined parameters ........................................................

path_semsIR = 'sems_20190407_07.fid'  # SEMS-IR .fid folder location [str]
path_mems = 'mems_20190406_01.fid'    # MEMS .fid folder location [str]
slices_semsIR = [0,1]                 # pick slices from IR experiment (first = 0) [list]
slices_mems = [5,2]                   # pick corresponding slices from MEMS experiment [list]
TI_wish = [100, 400, 1100, 2000]      # list of TI values (ms) for theoretical MRI images
TE_wish = [1, 10, 30, 60, 100]        # list of TE values (ms) for theoretical MRI images


glob_var = 0           # save the new images in a python global variable? [int]
                       # 0 - run without saving anything as a global variable
                       # 1 - run with saving the new images in a global variable
                       # 2 - run with saving the new images and the maps in global variables

# .......... End of user defined parameters .................................................

import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import shutil
import warnings
import joblib
from joblib.externals.loky import get_reusable_executor


def reorder_forT1IR(images, nI, number_of_images):
    shape = images.shape
    images = np.reshape(images, (nI, int(number_of_images/nI), shape[1], shape[2]))
    images = np.transpose(images, (1, 0, 2, 3))
    images = np.reshape(images, shape)
    return images


def get_images(number_of_images, echoes, ni_ne, slices): #interleave
    
    # dividing the data into k-spaces
    kspaces = np.empty((number_of_images, int(echoes.shape[0]/number_of_images), 
                        echoes.shape[1]), dtype=np.complex_)
    
    for kspace_no in range(0,number_of_images):
        
        kspaces[kspace_no] = echoes[kspace_no : echoes.shape[0] : number_of_images, :]
    
    # reordering T1-w data
    if slices is slices_semsIR:
        kspaces = reorder_forT1IR(kspaces, ni_ne, number_of_images)
        
    # reshaping to 4D do distinguish slices direction from and TI or TE direction
    kspaces = kspaces.reshape(int(kspaces.shape[0]/ni_ne), ni_ne, kspaces.shape[1], kspaces.shape[2])
    
    
    # picking slices
    kspaces = list(kspaces[i] for i in slices)
    kspaces = np.array(kspaces)
        
    # back to 3D data shape
    kspaces = kspaces.reshape(kspaces.shape[1] * kspaces.shape[0], kspaces.shape[2], kspaces.shape[3])
    
    # calculating images from k-spaces
    images = np.empty(kspaces.shape)
    
    for kspace_no in range(0,kspaces.shape[0]):
        ft = np.fft.fft2(kspaces[kspace_no])
        ft = np.fft.fftshift(ft)              # fixing problem with corner being center of the image
        ft = np.transpose(np.flip(ft, (1,0))) # matching geometry with VnmrJ-calculated image (still a bit shifted)
        images[kspace_no] = abs(ft)
        
    return images


def T1_function(x, T1, Mo, C, a):               # y = SI, x = TI, a = approx. 2
    x = np.array(x)
    y = abs(Mo * (1 - a* np.exp(-x/T1)) + C)
    return y


def T2_function(x, T2, Mo, C):                        # y = SI, x = TE
    x = np.array(x)
    y = Mo * np.exp(-x/T2) + C
    return y

def T1_T2_function(TI, TE, T1, T2, Mo, C1, C2):           # y = SI
    y = abs(Mo * (1 - 2* np.exp(-TI/T1)) + C1) * np.exp(-TE/T2) + C2
    return y


def calculate_maps(images, T_train, function):
    
    T_maps = []
    Momaps = []
    Cmaps = []
    nE_nI = len(T_train)
    

    for i in range(int(images.shape[0]/nE_nI)):
    
        a = nE_nI*i
        b = nE_nI*(i+1)
        
        slice1 = images[a:b]
        
        T_list = []
        Molist = []
        Clist = []
        
        for k in range(0,slice1.shape[1]):
                        
            def task(j):
                points = slice1[:, k, j]
                
                if function == T1_function:
                    bounds = ([0.001, points.max()*0.9, -(points.max()/100+1), 1.85],
                              [7000, 2*points.max()+1, points.max()/100+1, 2.05])
                else:
                    bounds = ([0.001, points.max()*0.9, -(points.max()/100+1)], 
                              [4000, 2*points.max()+1, points.max()/100+1])
        
                try:
                    parameters = curve_fit(function, T_train, points, 
                                           bounds=bounds, maxfev = 1000)[0]
                except RuntimeError:
                    if function == T1_function:
                        parameters = [0.000001,points.max(),0,2]
                    else:
                        parameters = [0.000001,points.max(),0]
                    
                return parameters
            
            with joblib.parallel_backend(backend="loky"):
                result = joblib.Parallel(n_jobs=-1)(joblib.delayed(task)(j) for j in range(0,slice1.shape[2]))
                
            result = np.array(result)
            T_list.append(result[:, 0])
            Molist.append(result[:, 1])
            Clist.append(result[:, 2])
        
        T_maps.append(np.reshape(np.array(T_list), (slice1.shape[1],slice1.shape[2])))
        Momaps.append(np.reshape(np.array(Molist), (slice1.shape[1],slice1.shape[2])))
        Cmaps.append(np.reshape(np.array(Clist), (slice1.shape[1],slice1.shape[2])))
                
    get_reusable_executor().shutdown() # close joblib processes
    
    return T_maps, Momaps, Cmaps


def theoret_MRI(function, TI_wish, TE_wish, T1_maps, T2_maps, Momaps, Cmaps_T1, Cmaps_T2):
    
    TheoretImgs = []
    out_folder = 'Theoretical_MRI'
    shutil.rmtree(out_folder, ignore_errors=True)  # removing residual output folder with content
    os.makedirs(out_folder)                        # creating new output folder
    
    
    for i in range(len(TI_wish)):
        
        for k in range(len(TE_wish)):
        
            for j in range(len(Momaps)):
                
                SI_image = function(TI_wish[i], TE_wish[k], T1_maps[j], T2_maps[j],
                                    Momaps[j], Cmaps_T1[j], Cmaps_T2[j])
                    
                TheoretImgs.append(SI_image)
                
                SI_image *= 255.0/SI_image.max()
                plt.imsave(fname=f'Theoretical_MRI/slice_{j}_TI_{TI_wish[i]}ms_TE_{TE_wish[k]}ms.png', 
                           arr=SI_image, cmap='gray')
    
    return TheoretImgs



def main(path_semsIR, path_mems, TI_wish, TE_wish, slices_semsIR, slices_mems, glob_var):
    
    # deleting global variables - already have them as local
    del globals()['path_semsIR']
    del globals()['path_mems']
    del globals()['TI_wish']
    del globals()['TE_wish']
    del globals()['slices_mems']
    del globals()['glob_var']
    
        
    # k-space data import (T1-w) with supressed nmrglue warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        paramsT1, echoesT1 = ng.agilent.read(dir=path_semsIR)
    
    # k-space data import (T2-w) with supressed nmrglue warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        paramsT2, echoesT2 = ng.agilent.read(dir=path_mems)
    
    del path_semsIR, path_mems

    # checking SEMS-IR data
    condition1 = paramsT1['procpar']['layout']['values'][0] == 'sems' and paramsT1['procpar']['ir']['values'][0] == 'y'
    
    if condition1:        
        print('(1/4) Recognised T1-weighted (SEMS-IR) images.')
        
    else:
        print('Error: no SEMS-IR data.')
        
    # checking MEMS data
    condition2 = paramsT2['procpar']['layout']['values'][0] == 'mems' and paramsT2['procpar']['ir']['values'][0] == 'n'
    
    if condition2:
        print('(2/4) Recognised T2-weighted (MEMS) images.')
        
    else:
        print('Error: no MEMS data.') 
        
        
    # Calculations   
    
    if condition1 and condition2:
        print('(3/4) Calculations in progress - please have patience...')
        
        # Calculations 1: T1, Mo(T1) and C(T1) maps
        
        # parameters of use    
        TI_train = np.array([eval(i) for i in (paramsT1['procpar']['ti']['values'])])*1000
        number_of_images = paramsT1['ntraces'] * len(TI_train)
        del paramsT1
        
        # calculating images from the data
        images = get_images(number_of_images, echoesT1, len(TI_train), slices_semsIR)
        del echoesT1
        
        # calculating parametric maps
        T1maps, Momaps_T1, Cmaps_T1 = calculate_maps(images, TI_train, T1_function)
        del TI_train
        
        
        # Calculations 2: T2, Mo(T2) and C(T2) maps
                
        # parameters of use  
        TE_train = np.array([eval(i) for i in (paramsT2['procpar']['TE']['values'])])
        number_of_images = paramsT2['ntraces']
        del paramsT2
               
        # calculating images from the data
        images = get_images(number_of_images, echoesT2, len(TE_train), slices_mems)
        del echoesT2, number_of_images
        del globals()['slices_semsIR']
        
        # calculating parametric maps
        T2maps, Momaps_T2, Cmaps_T2 = calculate_maps(images, TE_train, T2_function)
        del TE_train, images
        
        
        # Calculations 3: mean Mo maps
        Momaps_mean = np.mean(np.array([Momaps_T1, Momaps_T2]), axis=0)
        del Momaps_T1, Momaps_T2
        
        
        # Calculations 4: calculating and saving theoretical images
        TheoretImgs = theoret_MRI(T1_T2_function, TI_wish, TE_wish, T1maps, T2maps, 
                                  Momaps_mean, Cmaps_T1, Cmaps_T2)
        print('(4/4) New image(s) saved.')
        
                
        # Saving the data in global variables (optional)
        if glob_var == 1:
            return TheoretImgs
        
        elif glob_var == 2:
            return TheoretImgs, T1maps, T2maps, Momaps_mean
            
        
    # message after data error
    else:
        print('Calculations not passible.')


if __name__ == "__main__":
    
    if glob_var == 1:
        TheoretImgs = main(path_semsIR, path_mems, TI_wish, TE_wish, 
                           slices_semsIR, slices_mems, glob_var)
    
    elif glob_var == 2:
        TheoretImgs, T1_maps, T2_maps, Momaps = main(path_semsIR, path_mems, TI_wish, TE_wish, 
                                                     slices_semsIR, slices_mems, glob_var)
        
    else:
        main(path_semsIR, path_mems, TI_wish, TE_wish, slices_semsIR, slices_mems, glob_var)
