  
[![Build Status](https://travis-ci.org/mjirik/vebase.svg?branch=master)](https://travis-ci.org/mjirik/vebase)
[![Coverage Status](https://coveralls.io/repos/github/mjirik/vebase/badge.svg?branch=master)](https://coveralls.io/github/mjirik/vebase?branch=master)
[![PyPI version](https://badge.fury.io/py/vebase.svg)](http://badge.fury.io/py/vebase)


vebase

Vessel Basin Segmentation


```python
import vebase.livseg

p1 = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/"
p2 = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/portalvein/"
p1_i = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/image_1"
p2_i = "D:/liver_datasets_outputs/wokrin/1/MASKS_DICOM/liver/image_2"

# load ncsry data
l_data = vebase.livseg.load_vdata(p2,p1,p1_i,p2_i)
# build volumes
voda_ = vebase.livseg.vebavoda_sk(l_data)
# choose important nodes of portal vein
tree_red = vebase.livseg.tree_reduction(voda_[0],voda_[1], voda_[2],l_data[3]) 

print(tree_red)
```