import matplotlib.pyplot as plt
import matplotlib.image as mpimage
from mpl_toolkits import mplot3d
from loguru import logger

import pydicom as pyd

import numpy as np 

from skimage import exposure
import skimage.morphology as morp
from skimage.filters import rank

import os, glob
from pathlib import Path

from datetime import datetime

import natsort 

import skelet3d
import imma.labeled as imlb

from graphviz import Digraph

from collections import deque
#possible to remove some if needend


def load_vdata(p_mask_path, l_mask_path, l_mask_path_img1, l_mask_path_img2, accepted_suffixes=(".dcm")):
    """
    #load 3d masks liver and porta
    #use path arguments in format 
    #**win "D:liver_datasets_outputs/MASKS_DICOM/liver/" similar for porta
    #and two consecutive images like image5, image6
    #exmpl: "D:liver_datasets_outputs/MASKS_DICOM/liver/img5"
    #example: load_vdata(path1_porta, path2_liver, liver_image103, liverimage104)
    
    FUNCTION returns array: [volumecticDataPorta, volumetricDataLiver, spacingvalue]
    """    
    ds1 = pyd.read_file(l_mask_path_img1)
    ds2 = pyd.read_file(l_mask_path_img2)
    vox_sz = ds2.ImagePositionPatient[2]-ds1.ImagePositionPatient[2]
    logger.debug("masks readed")
    #liver_mask
    PathDicom = l_mask_path
    lstFilesDCM = [] 
    lstFilesDCM = os.listdir(PathDicom) 
    sortedlist = natsort.natsorted(lstFilesDCM, reverse=False)
    fdcml = []
    for x in range(0,len(sortedlist)):
        if Path(sortedlist[x]).suffix in accepted_suffixes:
            fdcml.append(l_mask_path + sortedlist[x])
    # loop through all the DICOM files
    dcmmatrx = []
    for filenameDCM in fdcml:
        # read the file
        ds = pyd.read_file(filenameDCM)
        # store the raw image data
        dcmmatrx.append(ds.pixel_array)
        logger.debug(f"fn={filenameDCM}")
    test = dcmmatrx[1]
    tempsz = (np.shape(test))
    ylab = tempsz[1]
    xlab = tempsz[0]
    volume_data_liver = np.zeros([len(dcmmatrx), xlab ,ylab], dtype=np.int)
    counter = 0
    volumecekr = 0
    for pic in range(0, len(dcmmatrx)):
        temppi = dcmmatrx[pic]
        tempbox = [np.zeros([xlab, ylab])]
        logger.debug(f"liver frame={pic}")
        for x in range(0, xlab):
            for y in range(0, ylab):
                if(temppi[x, y]) == 255:
                    volumecekr = volumecekr + 1
                    volume_data_liver [int(counter), int(x), int(y)] = int(1)
        counter = counter + 1
    #porta_mask
    logger.debug("porta mask")
    PathDicom = p_mask_path
    lstFilesDCM = [] 
    lstFilesDCM = os.listdir(PathDicom) 
    sortedlist = natsort.natsorted(lstFilesDCM,reverse=False)
    fdcml = []
    for x in range(0,len(sortedlist)):
        if Path(sortedlist[x]).suffix in accepted_suffixes:
            fdcml.append(p_mask_path + sortedlist[x])
    # loop through all the DICOM files
    dcmmatrx = []
    for filenameDCM in fdcml:
        # read the file
        ds = pyd.read_file(filenameDCM)
        # store the raw image data
        dcmmatrx.append(ds.pixel_array)
    test = dcmmatrx[1]
    tempsz = (np.shape(test))
    ylab = tempsz[1]
    xlab = tempsz[0]
    volume_data_porta = np.zeros([len(dcmmatrx), xlab ,ylab], dtype=np.int)
    counter = 0
    for pic in range(0, len(dcmmatrx)):
        temppi = dcmmatrx[pic]
        tempbox = [np.zeros([xlab, ylab])]
        for x in range(0, xlab):
            for y in range(0, ylab):
                if(temppi[x, y]) == 255:
                    volume_data_porta [int(counter), int(x), int(y)] = int(1)
        counter = counter + 1
    vox_1 = ds1["PixelSpacing"][0]
    vox_2 = ds1["PixelSpacing"][1]
    print("Approximate volume of chosen liver is:", vox_sz*vox_1*vox_2*volumecekr," mm^3")
    load_array = []
    load_array.append(volume_data_porta)
    load_array.append(volume_data_liver)
    vox_sz = 1 
    load_array.append(vox_sz)
    load_array.append(volumecekr)
    return load_array



# def voda_sk(l_d):
def voda_sk(organ_seg, liver_seg, voxelsize, cr):
    """
    build skeleton of porta
    voda_sk(xyz) build distance 3D map and compute volumes of liverpart for each skeleton element of chosen dataset
    voda_sk(xyz) takes argument from load_vdata - it need volumetric data of porta, liver and vox size (not ncsry..)
    """
    l_d = [organ_seg, liver_seg, voxelsize, cr]

    tempsz = (np.shape(l_d[0][1]))
    ylab = tempsz[1]
    xlab = tempsz[0]
    
    volume_data = l_d[0]
    #skeletonization - params are saved in var stats
    skelet = skelet3d.skelet3d(volume_data)
    skan = skelet3d.SkeletonAnalyser(skelet, volume_data = volume_data, voxelsize_mm = [1, 1, 1])
    stats = skan.skeleton_analysis()
    edge_number = 1

    #build of 3D distance map
    linesarr_x = []
    linesarr_y = []
    linesarr_z = []
    temparr = skan.sklabel
    counter = 0
    tempslx_l = 0
    for i in range(0,len(l_d[0])):
        for j in range(0, xlab):
            for k in range(0, ylab):
                if (temparr[i][j][k] != 0):
                    linesarr_x.append(k)
                    linesarr_y.append(j)
                    counter = counter + 1
        for z in range(0, counter):
            tempslx_l = i * l_d[2] #third value in l_data is space between slices of current datasete
            linesarr_z.append(tempslx_l)
            counter = 0
    distance_map = imlb.distance_segmentation(temparr)
    nodeArr = []
    for i in range(1, (len(stats)) + 1):
        #print(i)
        try:
            nodeArr.append((stats[i]['nodeA_ZYX']))
        except:
            None
        try:
            nodeArr.append((stats[i]['nodeB_ZYX']))
        except:
            None
    zdata = []
    xdata = []
    ydata = []
    temp_xp = 2
    temp_yp = 1
    temp_zp = 0
    for i in range(len(nodeArr)):
        zdata.append(nodeArr[i][temp_zp])
        ydata.append(nodeArr[i][temp_yp])
        xdata.append(nodeArr[i][temp_xp])
        zdata[i] = zdata[i]*l_d[2]
    #build list of all elements in skeleton
    list_= []
    for slx in range(len(l_d[1])):
        for x in range(0, xlab):
            for y in range(0, ylab):
                if(distance_map[slx, x, y]) < 1e5:
                    list_.append(distance_map[slx, x, y])
    set(list_)
    list_of_areas = set(list_)
    list_of_areas_arr = list(list_of_areas)
    list_of_areas_arr_edges = []
    for x in range(0, len(list_of_areas_arr)):
        if(list_of_areas_arr[x] > 0):
            list_of_areas_arr_edges.append(list_of_areas_arr[x])

    #plot data if needed          
    slicearr_p = []
    zte_p = []
    xte_p = []
    yte_p = []
    counter = 0

    for slx in range(len(l_d[0])):
        for x in range(0, xlab):
            for y in range(0, ylab):
                if(l_d[0][slx, x, y]) == 1:
                    xte_p.append(y)
                    yte_p.append(x)
                    counter = counter + 1
        for z in range(0, counter):
            tempslx_p = slx * l_d[2]
            zte_p.append(tempslx_p)
        counter = 0
    slicearr_p.append(xte_p)
    slicearr_p.append(yte_p)
    slicearr_p.append(zte_p)

    slicearr_l = []
    zte_l = []
    xte_l = []
    yte_l = []
    counter = 0

    for slx in range(len(l_d[1])):
        for x in range(0, xlab):
            for y in range(0, ylab):
                if(l_d[1][slx, x, y]) == 1:
                    xte_l.append(y)
                    yte_l.append(x)
                    counter = counter + 1
        for z in range(0, counter):
            tempslx_l = slx * l_d[2]
            zte_l.append(tempslx_l)
        counter = 0
    slicearr_l.append(xte_l)
    slicearr_l.append(yte_l)
    slicearr_l.append(zte_l)

    dist_map_x = []
    dist_map_y = []
    dist_map_z = []
    counter_sl = 0

    #compute of ALL VOLUMES - takes time - control print for each skelet element
    dist_map_final_liver_vol_temp = []
    dist_map_final_liver_vol = []

    loar_l = len(list_of_areas_arr)
    if(loar_l % 2 == 0):
        loar_l = loar_l/2
    else:
        loar_ = loar_l/2 + 1


    for area in range(0, loar_l):
        temp_area = list_of_areas_arr[area]
        for slc in range(0, len(l_d[1])):
            for x in range(0, xlab):
                for y in range(0, ylab):
                    if ((l_d[1][slc, x, y]) == 1):
                        if(((distance_map[slc, x, y]) == temp_area)):
                            dist_map_x.append(y)
                            dist_map_y.append(x)
                            counter_sl = counter_sl + 1
            for slic in range(0, counter_sl):
                dist_map_z.append(slc)
            counter_sl = 0
        dist_map_final_liver_vol_temp.append(dist_map_z)
        dist_map_final_liver_vol_temp.append(dist_map_y)
        dist_map_final_liver_vol_temp.append(dist_map_x)
        dist_map_final_liver_vol_temp.append(temp_area)
        dist_map_final_liver_vol.append(dist_map_final_liver_vol_temp)
        dist_map_final_liver_vol_temp = []
        dist_map_x = []
        dist_map_y = []
        dist_map_z = []
        print("computing volume of area: ", temp_area, "Num. of areas in this dataset: ", len(list_of_areas_arr))
        print("$$$$$-------------------------------------------------------------------$$$$$")

    # ret_array = []
    # ret_array.append(stats)
    # ret_array.append(list_of_areas_arr_edges)
    # ret_array.append(dist_map_final_liver_vol)
    return stats, list_of_areas_arr_edges, dist_map_final_liver_vol
    # return ret_array

def create_labeling(input_image_shape) -> np.ndarray:
    """
    Create ndarray with segment labeling
    """


    labeled = np.zeros_like(input_image_shape)

    return labeled
    
def tree_reduction(stats, list_of_areas_arr_edges, dist_map_final_liver_vol, volumecekr, tree_graph_check, res_div_nodes = [1],seg_params = [0.47,0.76]):
    
    if tree_graph_check == 1:
        from graphviz import Digraph
    """
    tree_reduction merge volumetric data to tree and remove pot. unnecessary nodes - return array with pot. IMPORTANT nodes
    represent connections in var mattest
    """
    #testlistn contains every conncted pair of nodes
    #mattest prepresent structure of stats 2D matrix
    id_ = []
    for i in range(0, (len(stats))):
        try:
            id_.append(stats[i+1]["nodeIdA"])
            id_.append(stats[i+1]["nodeIdB"])
        except:
            None
    id_res = list(set(id_))
    id_res.sort(reverse = True)
    #matrix tree representation
    edgcounter = []
    for i in range(0, len(stats)):
        i = i + 1
        for j in range(0, len((stats[i]["connectedEdgesA"]))):
            try:
                if(stats[i]["connectedEdgesA"][j]) != None:
                    edgcounter.append(stats[i]["connectedEdgesA"][j])
                else:
                    None
            except:
                None
                try:
                    if(stats[i]["connectedEdgesB"]) != None:
                        edgcounter.append(stats[i]["connectedEdgesB"][j])
                    else:
                        None
                except:
                    None
    mattest = np.zeros([len(id_res) + 1, len(list_of_areas_arr_edges) + 1], dtype=np.int)
    err = 0
    for i in range(0, (len(stats))+1):
        try:
            mattest[i + 1][0] = id_res[i]
        except:
            err = err + 1
        try:
            mattest[i + 2][0] = id_res[i+1]
        except:
            err = err + 1
    err = 0
    for i in range(0, len(list_of_areas_arr_edges)):
        try:
            mattest[0][i + 1] = list_of_areas_arr_edges[i]
        except:
            err = err + 1
    tempsz_m = (np.shape(mattest))
    ylab_ = tempsz_m[1]
    xlab_ = tempsz_m[0]
    for i in range(0, len(stats)):
        i = i + 1
        try:
            x_a = stats[i]["nodeIdA"]
        except:
            None
        try:
            x_b = stats[i]["nodeIdB"]
        except:
            None
        temp_e = stats[i]["id"]
        try:
            mattest[-x_a][temp_e] = 1
        except:
            None
        try:
            mattest[-x_b][temp_e] = 1
        except:
            None
    for i in range(0, len(stats)):
        i = i + 1
        try:
            x_a = stats[i]["nodeIdA"]
        except:
            None
        try:
            x_b = stats[i]["nodeIdB"]
        except:
            None

        try:
            v_edg_a = stats[i]["connectedEdgesA"]
        except:
            None
        try:
            v_edg_b = stats[i]["connectedEdgesB"]
        except:
            None

        for j in range(len(v_edg_a)):
            mattest[-x_a][v_edg_a[j]] = 1
        for k in range(len(v_edg_b)):
            mattest[-x_b][v_edg_b[k]] = 1
    area_id_vol = []
    for i in range(0, len(id_res)):
        temp = [id_res[i]]
        area_id_vol.append(temp)
    # find ,,deepest,, element, and its index in skeleton etc. (works for human liver if CTed in a standard way)
    temp = 1000
    rootid = 0
    _id_ = 0
    neighbourNode = 0
    edg_v_size = 0
    #find root
    for i in range (len(stats)):
        i = i + 1
        try:
            if stats[i]["nodeA_ZYX"][0] < temp:
                temp = stats[i]["nodeA_ZYX"][0]
                rootid = stats[i]["nodeIdA"]
                _id_ = stats[i]["id"]
                neighbourNode = stats[i]["nodeIdB"]
        except:
            None
        try:
            if stats[i]["nodeB_ZYX"][0] < temp:
                temp = stats[i]["nodeB_ZYX"][0]
                rootid = stats[i]["nodeIdB"]
                _id_ = stats[i]["id"]
                neighbourNode = stats[i]["nodeIdA"]
        except:
            None
        try:
            if len(stats[i]["connectedEdgesA"]) > edg_v_size:
                edg_v_size = len(stats[i]["connectedEdgesA"])
        except:
            None
        try:    
            if len(stats[i]["connectedEdgesB"]) > edg_v_size:
                edg_v_size = len(stats[i]["connectedEdgesB"])
        except:
            None
    root_node_p = temp
    realroot = rootid

    """
    portal vein tree builder with graphics representation *graphviz
    return pot. important nodes of skeleton
    """
    arofedges_buildvol = []
    for i in range(0,len(dist_map_final_liver_vol)):
        check = (dist_map_final_liver_vol[i][3])
        if (check > 0):
            arofedges_buildvol.append(dist_map_final_liver_vol[i])
    arofedges_buildvolfinal = list(arofedges_buildvol)
    testlistn = []
    tedgs = []
    for asfd in range(0, len(id_res)):
        for i in range(1, len(mattest[1, :])):
            if mattest[id_res[asfd], i] != 0:
                tedgs.append(i)
    #connected edges
    tnods = []
    for x in range(0, len(tedgs)):
        tmpnar = []
        for i in range(1, len(mattest[:, 1])):
            if (mattest[i, tedgs[x]]) != 0:
                tmpnar.append(i)
        tnods.append(tmpnar)
        tmpnar = []
    for x in range(len(tnods)):
        if tnods[x] in testlistn:
            None
        else:
            testlistn.append(tnods[x])
    nodelisttoadd = []
    for i in range(0, len(testlistn)):
        for j in range(0, len(testlistn[i])):
            chck = testlistn[i]
            if chck[j] in nodelisttoadd:
                None
            else:
                nodelisttoadd.append(chck[j])
    objcts_atrbs = []
    for i in range(0, len(id_res)):
        temp = [id_res[i]]
        objcts_atrbs.append(temp)
    matt_sz = np.shape(mattest)

    ta = []
    for i in range(1,-min(id_res)+1):
        ta.append(-i)
    uzly = []
    ttt = -1
    for j in range(0,len(ta)):
        for i in range(0,len(dist_map_final_liver_vol)):
            if dist_map_final_liver_vol[i][3] == ttt:
                uzly.append(dist_map_final_liver_vol[i])
                ttt = ttt - 1
    for i in range(0, len(objcts_atrbs)):
        objcts_atrbs[i].append(len(uzly[i][0])) 
    for i in range(0,len(testlistn)):
        sumtemp = 0
        edgindex = -1
        ind1 = testlistn[i][0]
        try:
            ind2 = testlistn[i][1]
        except:
            None
        for j in range(0,len(mattest[0][:])):
            if mattest[ind1][j] == 1 and mattest[ind2][j] == 1:
                edgindex = j
                try:
                    objcts_atrbs[ind1-1].append(float(len(arofedges_buildvolfinal[edgindex-1][0]))/2)
                except:
                    None
                try:
                    objcts_atrbs[ind2-1].append(float(len(arofedges_buildvolfinal[edgindex-1][0]))/2)
                except:
                    None
    for i in range(0, len(objcts_atrbs)):
        tmp = sum(objcts_atrbs[i][1:len(objcts_atrbs[i])])
        while(len(objcts_atrbs[i]) > 1):
            objcts_atrbs[i].remove(objcts_atrbs[i][-1])
        objcts_atrbs[i].append(tmp)
        tmp = 0
    testos = []
    for i in range(0, len(objcts_atrbs)):
        testos.append(objcts_atrbs[i][1])

    class Node(object): 
        # Constructor to create a new node 
        def __init__(self, key, volume): 
            self.val = key 
            self.vol = volume
            self.left = None
            self.right = None

        def display(self):
            lines, _, _, _ = self._display_aux()
            for line in lines:
                print(line)
                
    def build_segments(root, arr, lastnode):
        if root==None:
            return
        if root == lastnode:
            None
        else:
            try:
                build_segments(root.left,arr,lastnode)
            except:
                None
            #print(root.val, end=" ")
            try:
                arr.append(root.val)
            except:
                None
            try:
                build_segments(root.right,arr,lastnode)
            except:
                None
        return arr
    def build_last(root, arr):
        if root==None:
            return
        else:
            build_last(root.left,arr)
            arr.append(root.val)
            build_last(root.right,arr)
        return arr
    def build_last_only(root, arr):
        if root==None:
            return
        else:
            arr.append(root.val)
        return arr
    def insert(target, node, l_r):
        log = 0
        if target is None: 
            target = node
        elif target.left == None: 
            if l_r == "l": 
                if target.left is None: 
                    target.left = node 
                else: 
                    insert(target.left, node, l_r)
            else: 
                if target.right is None: 
                    target.right = node 
                else: 
                    insert(target.right, node, l_r) 
        elif target.right == None: 
            if l_r == "l": 
                if target.left is None: 
                    target.left = node 
                else: 
                    insert(target.left, node, l_r)
            else: 
                if target.right is None: 
                    target.right = node 
                else: 
                    insert(target.right, node, l_r) 
        elif target.right != None and target.left != None:
            insert(target.left, node, l_r)
            log = 1
        if log == 1:
            print("WARNING: trifurcation founded!!! location:",target.val)

    def visualize_tree_vol(tree):
        def add_nodes_edges(tree, dot=None):
            if dot is None:
                dot = Digraph()
                dot.node(name=str(tree), label=str(tree.vol))
            if tree.left:
                dot.node(name = str(tree.left) ,label = str(tree.left.vol))
                dot.edge(str(tree), str(tree.left))
                dot = add_nodes_edges(tree.left, dot=dot)
            if tree.right:
                dot.node(name = str(tree.right) ,label = str(tree.right.vol))
                dot.edge(str(tree), str(tree.right))
                dot = add_nodes_edges(tree.right, dot=dot)
            return dot
        dot = add_nodes_edges(tree)
        # display(dot)
        return dot

    def visualize_tree_val(tree):
        def add_nodes_edges(tree, dot=None):
            if dot is None:
                dot = Digraph()
                dot.node(name=str(tree), label=str(tree.val))
            if tree.left:
                dot.node(name = str(tree.left) ,label = str(tree.left.val))
                dot.edge(str(tree), str(tree.left))
                dot = add_nodes_edges(tree.left, dot=dot)
            if tree.right:
                dot.node(name = str(tree.right) ,label = str(tree.right.val))
                dot.edge(str(tree), str(tree.right))
                dot = add_nodes_edges(tree.right, dot=dot)
            return dot
        dot = add_nodes_edges(tree)
        # display(dot)
        return dot

    class SumTree:
        def toSumTree(self, root, los): 
             # Base case 
            if root is None:
                return  0
            old_val = root.vol
            root.vol=self.toSumTree(root.left,los) + self.toSumTree(root.right,los) 
            return  root.vol + old_val

        def print_val(self, root):
            print("Printed value of node: ", root.val)

        def ret_val(self, root):
            return root.val

        def traverse_inorder_nodeid(self,root):
            if root==None:
                return
            self.traverse_inorder_nodeid(root.left)
            #print(root.val, end=" ")
            self.traverse_inorder_nodeid(root.right)

        def traverse_inorder_nodevol(self,root,los):
            if root==None:
                return
            self.traverse_inorder_nodevol(root.left,los)
            #print(root.vol,end=" ")
            t = []
            t.append(root.val)
            t.append(root.vol)
            los.append(t)
            t = []
            self.traverse_inorder_nodevol(root.right,los) 

    def find_bro(arr, val):

        if len(arr) == 1:
            return arr[0]
        if (arr[0] == val):
            return arr[1]
        if (arr[1] == val):
            return arr[0]

    def sort_fcn(e):
        return len(e)

    def main(res_div_nodes = res_div_nodes):
        #clasification important___bifurcation / notbifurcation
        def find_all_div_nodes(los, seg_params):
            div_nodes = []
            for i in range(0,len(los)):
                tempvol = potroots[los[i][0]-1].vol
                try:
                    if potroots[los[i][0]-1].left.vol > 150000 and (potroots[los[i][0]-1].left.vol > tempvol*seg_params[0] and potroots[los[i][0]-1].left.vol < tempvol*seg_params[1]):
                        div_nodes.append(potroots[los[i][0]-1].val)
                except:
                    None
                try:
                    if potroots[los[i][0]-1].right.vol > 150000 and (potroots[los[i][0]-1].right.vol > tempvol*seg_params[0] and potroots[los[i][0]-1].right.vol < tempvol*seg_params[1]):
                        div_nodes.append(potroots[los[i][0]-1].val)
                except:
                    None
            return div_nodes

        st = SumTree()
        addednodeslist = []
        #LISTOFNODES!!!
        potroots = []
        for i in range (0, len(objcts_atrbs)):
            potroots.append(Node(-objcts_atrbs[i][0], objcts_atrbs[i][1]))
        #start of tree
        root = potroots[-rootid-1]
        addednodeslist.append(-rootid)
        wta = []
        addtemparr = []   
        #potentional new nodes
        for i in range(0, len(testlistn)):
            if st.ret_val(potroots[-rootid-1]) in testlistn[i]:
                addtemparr.append(testlistn[i])
        #adding root mates
        for i in range(0, len(addtemparr)):
            nodetoadd = find_bro(addtemparr[i], st.ret_val(potroots[-rootid-1]))
            if nodetoadd in addednodeslist:
                None
            else:
                if potroots[-rootid-1].left == None:
                    insert(potroots[-rootid-1], potroots[nodetoadd-1], "l")
                    addednodeslist.append(nodetoadd)
                    nodelisttoadd.remove(nodetoadd)
                    wta.append(nodetoadd-1)

                else:
                    insert(potroots[-rootid-1], potroots[nodetoadd-1], "r")
                    addednodeslist.append(nodetoadd)
                    nodelisttoadd.remove(nodetoadd)
                    wta.append(nodetoadd-1)
        for b in range(0,len(objcts_atrbs)):
            addtemparr = []  
            tempwta = []
            for a in range(0, len(wta)):
                addtemparr = []   
                #potentional new nodes
                for i in range(0, len(testlistn)):
                    if st.ret_val(potroots[wta[a]]) in testlistn[i]:
                        addtemparr.append(testlistn[i])
                #adding root mates
                for i in range(0, len(addtemparr)):
                    nodetoadd = find_bro(addtemparr[i], st.ret_val(potroots[wta[a]]))
                    if nodetoadd in addednodeslist:
                        None
                    else:
                        if potroots[wta[a]].left == None:
                            insert(potroots[wta[a]], potroots[nodetoadd-1], "l")
                            addednodeslist.append(nodetoadd)
                            nodelisttoadd.remove(nodetoadd)
                            tempwta.append(nodetoadd-1)
                        else:
                            insert(potroots[wta[a]], potroots[nodetoadd-1], "r")
                            addednodeslist.append(nodetoadd)
                            nodelisttoadd.remove(nodetoadd)
                            tempwta.append(nodetoadd-1)
            wta = tempwta[:]
            tempwta = []
        los = []
        st.toSumTree(root,los)
        st.traverse_inorder_nodevol(root, los)
        los_c = []
        for i in range(0,len(los)):
            los_c.append(los[i][1])
        tempmax = max(los_c)
        potdiv = 0.95*tempmax #key to find first bifurcation
        indexarray = []
        valarray = []
        for i in range(0,len(los_c)):
            if los_c[i] > potdiv:
                indexarray.append(i)
                valarray.append(los_c[i])
        minval_ = min(valarray)
        for i in range(0,len(los)):
            if los_c[i] == minval_:
                p_newroot = los_c[i]
                break
        indexofp_root = -9999
        for i in range(0,len(potroots)):
            if p_newroot == (potroots[i].vol):
                indexofp_root = i
        root = potroots[indexofp_root]
        st.traverse_inorder_nodeid(root)
        st.traverse_inorder_nodevol(root,los)
        
        if tree_graph_check == 1:
            #usin graphviz for visual validation - can be turned off fcn will retun only pot. important nodes
            #first bifuracation
            #tree_pic = visualize_tree_val(root)  #visualize_tree_vol/val 

            orgroot = potroots[-rootid-1]
            #whole tree
            tree_pic = visualize_tree_val(orgroot)
            # Save as png file
            tree_pic.format = 'png'
            tree_pic.view(filename = 'Vesseltree_', directory='C:/Users/user/Desktop/testprj')
            tree_pic = visualize_tree_vol(orgroot)
            # Save as png file
            tree_pic.format = 'png'
            tree_pic.view(filename = 'Vesseltree_1', directory='C:/Users/user/Desktop/testprj')
        orgroot = potroots[-rootid-1]
        
        if len(res_div_nodes) == 1:
            res_div_nodes = find_all_div_nodes(los,seg_params)
        else:
            None
            
        res_div_nodes_n = list(set(res_div_nodes))
        final_array_segments = []
        #fcn to get specific segments for plot
        for i in range(0,len(res_div_nodes_n)):
            temp_arx = []
            try:
                final_array_segments.append(build_segments(potroots[res_div_nodes_n[i]-1].left,temp_arx,root))
            except:
                None
            try:
                final_array_segments.append(build_segments(potroots[res_div_nodes_n[i]-1].right,temp_arx,root))
            except:
                None
        temp_arx = []
        final_array_segments.append(build_last_only(root,temp_arx))
        print("Index of first bifurcation:", root.val)
        if len(res_div_nodes) == 1:
            final_array_segments.sort(key=sort_fcn)
            final_array_segments.reverse()            
        else:
            None
        final_array_segments.sort(key=sort_fcn)
        final_array_segments.reverse()            

            
        if None in final_array_segments:
            final_array_segments.remove(None)
        res_div_nodes_n.append(final_array_segments)
        print("Chosen important bifurcations and their adjacent areas are: ",res_div_nodes_n)
        return res_div_nodes_n
    
    los = main()
    #function to check the difference between volume of liver and summs
    #objcts_atrbs[1]
    chc = 0
    for i in range(0, len(objcts_atrbs)):
        chc = chc + objcts_atrbs[i][1]
    print("Pixel difference between final volumetric tree and input volume of liver: ", volumecekr - chc, "pix")
    arr_return = []
    arr_return.append(los)
    return arr_return

#fct. which gives u hit about segmentation parametres
def seg_hinter(req_seg_num):
    if req_seg_num == 2: #1. bifurcation
        params = [0.6,0.885]
        return params
    elif req_seg_num == 4: #2. bifurcation
        params = [0.5,0.77]
        return params
    elif req_seg_num > 8:
        params = [0.1,1]
        return params
    else:
        print("execute seg fcn. without specific input -> aprox. 8 segments ")

def vein_b_viz(porta,segs,stats,dist_map_final_liver_vol):
    id_ = []
    for i in range(0, (len(stats))):
        try:
            id_.append(stats[i+1]["nodeIdA"])
            id_.append(stats[i+1]["nodeIdB"])
        except:
            None
    id_res = list(set(id_))
    id_res.sort(reverse = True)
    arofedges_buildvol = []
    for i in range(0,len(dist_map_final_liver_vol)):
        check = (dist_map_final_liver_vol[i][3])
        if (check > 0):
            arofedges_buildvol.append(dist_map_final_liver_vol[i])
    arofedges_buildvolfinal = list(arofedges_buildvol)

    ta = []
    for i in range(1,-min(id_res)+1):
        ta.append(-i)
        
    uzly = []
    ttt = -1
    for j in range(0,len(ta)):
        for i in range(0,len(dist_map_final_liver_vol)):
            if dist_map_final_liver_vol[i][3] == ttt:
                uzly.append(dist_map_final_liver_vol[i])
                ttt = ttt - 1
    arr_objcts = []
    seg_frac = []
    for i in range(0,len(segs)):
        seg_frac.append(segs[i])
    for i in range(0,len(seg_frac)):
        ata = []
        for j in range(0,len(seg_frac[i])):
            test_val = -seg_frac[i][j]
            for k in range(1,len(stats)):
                try:
                    if test_val == stats[k]["nodeIdA"] or test_val == stats[k]["nodeIdB"]:
                        if arofedges_buildvol[k-1] in ata:
                            None
                        else:
                            ata.append(arofedges_buildvol[k-1])
                except:
                    if test_val == stats[k]["nodeIdA"]:
                        if arofedges_buildvol[k-1] in ata:
                            None
                        else:
                            ata.append(arofedges_buildvol[k-1])
                    
            ata.append(uzly[-test_val-1])
        arr_objcts.append(ata)
    for i in range(0,len(arr_objcts)):
        for j in range(0,len(arr_objcts[i])):
            for k in range(0,len(arr_objcts[i][j][0])):
                if (porta[arr_objcts[i][j][0][k]][arr_objcts[i][j][1][k]][arr_objcts[i][j][2][k]]) != 1:
                    arr_objcts[i][j][0][k] = 1
                    arr_objcts[i][j][1][k] = 256
                    arr_objcts[i][j][2][k] = 256

    
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection='3d')
    
    color = [(0.1, 0.1, 0.1), (0.1, 0.1, 0.3),(0.1, 0.1, 0.6),(0.1, 0.1, 0.9),(0.1, 0.3, 0.1),(0.1, 0.6, 0.1),(0.1, 0.9, 0.1),(0.3, 0.1, 0.1),(0.6, 0.1, 0.1),(0.9, 0.1, 0.1),(0.3, 0.1, 0.3),(0.6, 0.1, 0.6),(0.9, 0.1, 0.9),(0.5, 0.5, 0.5)]

    color_check = 0
    for i in range(0,len(arr_objcts)):
        for cnt in range(0,len(arr_objcts[i])):
            try:
                ax.scatter3D(arr_objcts[i][cnt][2],arr_objcts[i][cnt][1],arr_objcts[i][cnt][0],color = color[i])
            except:
                None
    #fig.savefig('myimage_p.png', format='png', dpi=600)
         

def vein_b_viz_l(porta,segs,stats,dist_map_final_liver_vol):
    id_ = []
    for i in range(0, (len(stats))):
        try:
            id_.append(stats[i+1]["nodeIdA"])
            id_.append(stats[i+1]["nodeIdB"])
        except:
            None
    id_res = list(set(id_))
    id_res.sort(reverse = True)
    arofedges_buildvol = []
    for i in range(0,len(dist_map_final_liver_vol)):
        check = (dist_map_final_liver_vol[i][3])
        if (check > 0):
            arofedges_buildvol.append(dist_map_final_liver_vol[i])
    arofedges_buildvolfinal = list(arofedges_buildvol)

    ta = []
    for i in range(1,-min(id_res)+1):
        ta.append(-i)
        
    uzly = []
    ttt = -1
    for j in range(0,len(ta)):
        for i in range(0,len(dist_map_final_liver_vol)):
            if dist_map_final_liver_vol[i][3] == ttt:
                uzly.append(dist_map_final_liver_vol[i])
                ttt = ttt - 1
    #print((arofedges_buildvol[5]))
    #print((uzly[5]))
    
    arr_objcts = []
    
    seg_frac = []
    for i in range(0,len(segs)):
        seg_frac.append(segs[i])
    for i in range(0,len(seg_frac)):
        ata = []
        for j in range(0,len(seg_frac[i])):
            test_val = -seg_frac[i][j]
            for k in range(1,len(stats)):
                try:
                    if test_val == stats[k]["nodeIdA"] or test_val == stats[k]["nodeIdB"]:
                        if arofedges_buildvol[k-1] in ata:
                            None
                        else:
                            ata.append(arofedges_buildvol[k-1])
                except:
                    if test_val == stats[k]["nodeIdA"]:
                        if arofedges_buildvol[k-1] in ata:
                            None
                        else:
                            ata.append(arofedges_buildvol[k-1])
                            print(arofedges_buildvol[k-1])

            ata.append(uzly[-test_val-1])
        arr_objcts.append(ata)
    #print(arr_objcts[1])
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection='3d')
    
    color = [(0.1, 0.1, 0.1), (0.1, 0.1, 0.3),(0.1, 0.1, 0.6),(0.1, 0.1, 0.9),(0.1, 0.3, 0.1),(0.1, 0.6, 0.1),(0.1, 0.9, 0.1),(0.3, 0.1, 0.1),(0.6, 0.1, 0.1),(0.9, 0.1, 0.1),(0.3, 0.1, 0.3),(0.6, 0.1, 0.6),(0.9, 0.1, 0.9),(0.5, 0.5, 0.5)]

    color_check = 0
    for i in range(0,len(arr_objcts)):
        for cnt in range(0,len(arr_objcts[i])):
            try:
                ax.scatter3D(arr_objcts[i][cnt][2],arr_objcts[i][cnt][1],arr_objcts[i][cnt][0],color = color[i])
            except:
                None
    #fig.savefig('algoritmus_jatra.png', format='png', dpi=600)

def seg_3dnp(liver,segs,stats,dist_map_final_liver_vol,slices_n):
    id_ = []
    for i in range(0, (len(stats))):
        try:
            id_.append(stats[i+1]["nodeIdA"])
            id_.append(stats[i+1]["nodeIdB"])
        except:
            None
    id_res = list(set(id_))
    id_res.sort(reverse = True)
    arofedges_buildvol = []
    for i in range(0,len(dist_map_final_liver_vol)):
        check = (dist_map_final_liver_vol[i][3])
        if (check > 0):
            arofedges_buildvol.append(dist_map_final_liver_vol[i])
    arofedges_buildvolfinal = list(arofedges_buildvol)

    ta = []
    for i in range(1,-min(id_res)+1):
        ta.append(-i)
        
    uzly = []
    ttt = -1
    for j in range(0,len(ta)):
        for i in range(0,len(dist_map_final_liver_vol)):
            if dist_map_final_liver_vol[i][3] == ttt:
                uzly.append(dist_map_final_liver_vol[i])
                ttt = ttt - 1
    #print((arofedges_buildvol[5]))
    #print((uzly[5]))
    
    arr_objcts = []
    
    seg_frac = []
    for i in range(0,len(segs)):
        seg_frac.append(segs[i])
    for i in range(0,len(seg_frac)):
        ata = []
        for j in range(0,len(seg_frac[i])):
            test_val = -seg_frac[i][j]
            for k in range(1,len(stats)):
                try:
                    if test_val == stats[k]["nodeIdA"] or test_val == stats[k]["nodeIdB"]:
                        if arofedges_buildvol[k-1] in ata:
                            None
                        else:
                            ata.append(arofedges_buildvol[k-1])
                except:
                    if test_val == stats[k]["nodeIdA"]:
                        if arofedges_buildvol[k-1] in ata:
                            None
                        else:
                            ata.append(arofedges_buildvol[k-1])
                try:
                    if test_val == stats[k]["nodeIdB"]:
                        if arofedges_buildvol[k-1] in ata:
                            None
                        else:
                            ata.append(arofedges_buildvol[k-1])
                except:
                    None
            ata.append(uzly[-test_val-1])
        arr_objcts.append(ata)
    seg3d = np.zeros([slices_n, 512 ,512], dtype=np.int)
    index_arr = []
    index_i = 1
    for i in range(0,999):
        index_arr.append(i)
    for i in range(0,len(arr_objcts)):
        if i != 0:
            index_i = index_i + 1
        for j in range(0,len(arr_objcts[i])):
            for k in range(0,len(arr_objcts[i][j][0])):
                if i == 0:
                    if (liver[arr_objcts[i][j][0][k]][arr_objcts[i][j][1][k]][arr_objcts[i][j][2][k]]) == 1:
                        seg3d[arr_objcts[i][j][0][k],arr_objcts[i][j][1][k],arr_objcts[i][j][2][k]] = index_arr[index_i]
                        
    
    #fig = plt.figure(figsize=(10,10))
    #ax = plt.axes(projection='3d')
    
    color = [(0.1, 0.1, 0.1), (0.1, 0.1, 0.3),(0.1, 0.1, 0.6),(0.1, 0.1, 0.9),(0.1, 0.3, 0.1),(0.1, 0.6, 0.1),(0.1, 0.9, 0.1),(0.3, 0.1, 0.1),(0.6, 0.1, 0.1),(0.9, 0.1, 0.1),(0.3, 0.1, 0.3),(0.6, 0.1, 0.6),(0.9, 0.1, 0.9),(0.5, 0.5, 0.5)]


    color_check = 0
    #for i in range(0,len(arr_objcts)):
    #    for cnt in range(0,len(arr_objcts[i])):
    #        try:
    #            ax.scatter3D(arr_objcts[i][cnt][2],arr_objcts[i][cnt][1],arr_objcts[i][cnt][0],color = color[i])
    #        except:
    #            None
    #fig.savefig('jacard.png', format='png', dpi=600)
    return seg3d
         
