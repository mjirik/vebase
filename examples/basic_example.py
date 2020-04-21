#! /usr/bin/python
# -*- coding: utf-8 -*-

# import logging
# logger = logging.getLogger(__name__)
# from loguru import logger
# import pytest
# import os.path
import vebase.livseg
import io3d.datasets

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

l_data = vebase.livseg.load_vdata(p2,p1,p1_i,p2_i) # --- load ncsry data
voda_ = vebase.livseg.vebavoda_sk(l_data) # --- build volumes
tree_red = vebase.livseg.tree_reduction(voda_[0],voda_[1], voda_[2],l_data[3]) # --- choose important nodes of portal vein

print(tree_red)
