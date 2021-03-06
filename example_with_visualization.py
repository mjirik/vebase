from loguru import logger

import vebase
import vebase.livseg
import io3d.datasets

from vebase.livseg import load_vdata, voda_sk, tree_reduction #, plot_tree_reduction_3d
#steps to execute livseg.py functions

p1 = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/"
p2 = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/portalvein/"
p1_i = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/image_1"
p2_i = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/image_2"
p1   = io3d.datasets.join_path("medical/orig/3Dircadb1.1/MASKS_DICOM/liver/", get_root=True) + "/"
p2   = io3d.datasets.join_path("medical/orig/3Dircadb1.1/MASKS_DICOM/portalvein/", get_root=True) + "/"
p1_i = io3d.datasets.join_path("medical/orig/3Dircadb1.1/MASKS_DICOM/liver/", get_root=True) + "/image_1"
p2_i = io3d.datasets.join_path("medical/orig/3Dircadb1.1/MASKS_DICOM/liver/", get_root=True) + "/image_2"

logger.debug(f"{p1}, {p2}")
logger.debug(f"{p1_i}, {p2_i})")
# --- load ncsry data
l_data = load_vdata(p2,p1,p1_i,p2_i)

#--- build volumes
voda_ = voda_sk(*l_data)

#--- choose important nodes of portal vein
# tree_red = tree_reduction(voda_[0],voda_[1], voda_[2],l_data[3], 1) #last argument 1 - vizualize 2d png of tree 0 go without graphviz
# kwargs = {"other_param": 4, "A":4}
tree_red = tree_reduction(*voda_, 1) #last argument 1 - vizualize 2d png of tree 0 go without graphviz
# tree_red = tree_reduction(*voda_, 1) #last argument 1 - vizualize 2d png of tree 0 go without graphviz

#figures 2d 3d choice, 2d choice needs specific slice - if u want to plot, tree_reduction(a,b,c,d) must be executed!
tree_red = tree_reduction(voda_[0],voda_[1], voda_[2],l_data[3],0) 
# plot_tree_reduction_3d(tree_red)

#tree_red = tree_reduction(voda_[0],voda_[1], voda_[2],l_data[3],0) 
#plot_tree_reduction_2d(tree_red,66) #choose from z axis interval
