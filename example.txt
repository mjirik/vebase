#steps to execute livseg.py functions

"""
p1 = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/"
p2 = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/portalvein/"
p1_i = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/image_1"
p2_i = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/image_2"

# --- load ncsry data
l_data = load_vdata(p2,p1,p1_i,p2_i)

#--- build volumes
voda_ = voda_sk(l_data) 

#--- choose important nodes of portal vein
tree_red = tree_reduction(voda_[0],voda_[1], voda_[2],l_data[3], 1) #last argument 1 - vizualize 2d png of tree 0 go without graphviz 

#figures 2d 3d choice, 2d choice needs specific slice - if u want to plot, tree_reduction(a,b,c,d) must be executed!
tree_red = tree_reduction(voda_[0],voda_[1], voda_[2],l_data[3],0) 
plot_tree_reduction_3d(tree_red)

#tree_red = tree_reduction(voda_[0],voda_[1], voda_[2],l_data[3],0) 
#plot_tree_reduction_2d(tree_red,66) #choose from z axis interval
"""
