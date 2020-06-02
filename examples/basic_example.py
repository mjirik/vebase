#! /usr/bin/python
# -*- coding: utf-8 -*-

# import logging
# logger = logging.getLogger(__name__)
# from loguru import logger
# import pytest
# import os.path
import vebase.livseg
import io3d.datasets
import numpy as np

p1 = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/"
p2 = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/portalvein/"
p1_i = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/image_1"
p2_i = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/image_2"

p1   = io3d.datasets.join_path("medical/orig/3Dircadb1.1/MASKS_DICOM/liver/", get_root=True)
p2   = io3d.datasets.join_path("medical/orig/3Dircadb1.1/MASKS_DICOM/portalvein/", get_root=True)
p1_i = io3d.datasets.join_path("medical/orig/3Dircadb1.1/MASKS_DICOM/liver/image_1", get_root=True)
p2_i = io3d.datasets.join_path("medical/orig/3Dircadb1.1/MASKS_DICOM/liver/image_2", get_root=True)
# io3d.datasets.join_path("medical/orig/3Dircadb1.1/")

# print(f"toto je promenna p1={p1}")
# print(f"{=p1}")

liver, porta, vox_sz, volumecekr = vebase.livseg.load_vdata(p2,p1,p1_i,p2_i) # --- load ncsry data


datap_liver = io3d.read(p1, dataplus_format=True)
datap_porta = io3d.read(p2, dataplus_format=True)

liver = ((datap_liver["data3d"] > 0 )).astype(np.int)
porta = ((datap_porta["data3d"] > 0 )).astype(np.int)
vox_sz = 1
volumecekr = np.sum(liver)



voda_ = vebase.livseg.voda_sk(liver, porta, vox_sz, volumecekr) # --- build volumes
tree_red = vebase.livseg.tree_reduction(voda_[0],voda_[1], voda_[2],volumecekr,1) # --- choose important nodes of portal vein



print(tree_red)



"""
priprava na utery funkcni kod.. jen pridat vebase.livseg.etc..

p1_7 = "D:/liver_datasets_outputs/wokrin/3/MASKS_DICOM/liver/"
p2_7 = "D:/liver_datasets_outputs/wokrin/3/MASKS_DICOM/portalvein/"
p1_i_7 = "D:/liver_datasets_outputs/wokrin/3/MASKS_DICOM/liver/image_1"
p2_i_7 = "D:/liver_datasets_outputs/wokrin/3/MASKS_DICOM/liver/image_2"

l_data_7 = load_vdata(p2_7,p1_7,p1_i_7,p2_i_7) # --- load ncsry data

voda__7 = voda_sk(*l_data_7) # --- build volumes
#couinaud + graphviz
#tree_red_7 = tree_reduction(voda__7[0],voda__7[1], voda__7[2],l_data_7[3],1) #0/1 $$$ graphvis array of nodes!if want to specific ones $$$ filtr consts. lowerone, higher one set to 47% 76% 
#manual targeting
#tree_red_7 = tree_reduction(voda__7[0],voda__7[1], voda__7[2],l_data_7[3],1,res_div_nodes = [1],seg_params = [0.35,0.8]) 
"""
