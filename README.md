# Individual-tree isolation (treeiso) from terrestrial laser scanning

### **`CloudCompare Plugin was developed (created - May 2023, upgraded - February 2025)!`**

### **`Python version was developed now (January 2025)!`**

The University of Lethbridge - Department of Geography & Environment - Artemis Lab

Author - Zhouxin Xi (zhouxin.xi@uleth.ca) and Prof. Chris Hopkinson (c.hopkinson@uleth.ca)

Please cite:

<cite>Xi, Z.; Hopkinson, C. 3D Graph-Based Individual-Tree Isolation (_Treeiso_) from Terrestrial Laser Scanning Point Clouds. _Remote Sens_. **2022**, _14_, 6116. https://doi.org/10.3390/rs14236116</cite>


This tool relies on the cut-pursuit algorithm, please also consider citing:

<cite> Landrieu, L.; Obozinski, G. Cut Pursuit: Fast Algorithms to Learn Piecewise Constant Functions on General Weighted Graphs. SIAM J. Imaging Sci. 2017, 10, 1724–1766. [hal link](https://hal.archives-ouvertes.fr/hal-01306779)</cite>

<cite> Raguet, H.; Landrieu, L. Cut-pursuit algorithm for regularizing nonsmooth functionals with graph total variation. In Proceedings of the International Conference on Machine Learning, Stockholm, Sweden, 10–15 July 2018; Volume 80, pp. 4247–4256.</cite>

## Folder structure

    ├── data                                    # All raw data
    │   ├── LPine1_demo.laz                     # Example TLS data (just a subset from a TLS plot due to the file size limit)
    ├── Matlab                                  # treeiso Matlab source code 
    │   ├── treeiso.m                           # Main Matlab program
    │   ├── cutPursuit.m                        # Wrapper of L0 cut-pursuit mex code
    │   ├── jitknnsearch.m                      # Matlab jit accelerated knnsearch - fast knnsearch from small amount of points
    │   ├── overlapping.m                       # Used to calculate overlap ratio between a pair of crown convex hulls
    │   ├── L0_cut_pursuit.mexw64               # compiled L0 cut pursuit program
    │   ├── L0_cut_pursuit_segmentation.mexw64  # compiled L0 cut pursuit program
    │   ├── las2mat.mexw64                      # compiled las/laz file reading program
    │   ├── mat2las.mexw64                      # compiled las/laz file writing program
    ├── Python                                  # treeiso Python source code 
    │   ├── treeiso.py                          # Main python program
    │   ├── cut_pursuit_L2.py                   # Cut-pursuit python version (simplified and L2 norm only)
    │   ├── cut_pursuit_L2_replica_cpp.py       # Cut-pursuit python version (closer to the original C++ version but slower, only for information)
    ├── LICENSE
    └── README.md

### 1. Preprocessing requirements
If you use your own TLS data, it is recommended to clean the point clouds and remove the noise, e.g. using CloudCompare(Sober filter,NN=10,std=1)

It is also suggested to decimate the point cloud to reasonable resolution (~2cm)

Ground points must be removed prior to tree isolation


**
The mex files were all compiled under Windows 10 (64bit) Intel using Visual Studio 2019. If you need to compile on your own, please refer to the authors' sites for detailed steps. Links below:

[cutpursuit](https://github.com/loicland/cut-pursuit), [mastools](https://github.com/plitkey/matlas_tools)

### 2. treeiso

#### Matlab

Run Matlab code under Matlab\treeiso.m to isolate TLS trees

Type in the folder path where your laz files are located

The treeiso.m will exhaustively search all laz files and run the treeiso for each laz file

You can type in the data folder for a quick test. The treeiso.m will isolate trees from LPine1_demo.laz

Result laz file will be saved to Matlab\output\

The laz file will include three additional fields:

1. init_segs: show the initial segmentation from 3D cut-pursuit
2. interim_segs: show the 2nd-stage segmentation from 2D cut-pursuit
3. segs: final tree segmentation

Please be patient: processing a laz file of ~25MB costs about 10 minutes using Intel Core i7-9700K and 16GB RAM.

#### Python (Uploaded on Jan 25, 2025)

Run Python code under Python\treeiso.py to isolate TLS trees. 

`python treeiso.py`

Like the Matlab version, type in the folder path where your laz files are located.

The result laz file will be saved to the same input folder you provided.


#### Example isolated trees: 
|Raw TLS Example1| After treeiso isolation| Top view|
|:---:|:---:|:---:|
|<img width="635" alt="demo1_crop" src="https://user-images.githubusercontent.com/8785889/182312969-7c81949f-67fa-409b-bb24-73b094917c52.png">|<img width="635" alt="demo1_treeiso_crop" src="https://user-images.githubusercontent.com/8785889/182313130-7d6ba091-b2cb-4482-ae21-3ee2ec5bb34e.png">|<img width="700" alt="demo1_treeiso2_crop" src="https://user-images.githubusercontent.com/8785889/182313150-a14e7a3e-79b0-4d40-a400-ce89520e5e4d.png">|


|Raw TLS Example2| After treeiso isolation| Top view|
|:---:|:---:|:---:|
|<img width="1354" alt="TAspen1_crop" src="https://user-images.githubusercontent.com/8785889/182315836-78bdd338-2678-4d29-bcce-1442777e6918.png">|<img width="1354" alt="TAspen1_treeiso_crop" src="https://user-images.githubusercontent.com/8785889/182315950-f44937e2-0ec8-4b43-a980-0da9b7925534.png">|<img width="724" alt="TAspen1_treeiso2_crop" src="https://user-images.githubusercontent.com/8785889/182315979-d8e9c5ea-7a00-4490-8bab-c43fd6af9498.png">|


|Raw TLS Example3| After treeiso isolation| Top view|
|:---:|:---:|:---:|
|<img width="1228" alt="NPoplar2_crop" src="https://user-images.githubusercontent.com/8785889/182318128-ddeac093-f48d-4e69-a9b2-560f5a25335e.png">|<img width="1228" alt="NPoplar2_treeiso_crop" src="https://user-images.githubusercontent.com/8785889/182318163-290e0663-4e0a-455c-85d5-5160f00a8e33.png">|<img width="1026" alt="NPoplar2_treeiso2_crop" src="https://user-images.githubusercontent.com/8785889/182318398-fc428e17-1fd7-415c-9594-e0b9e3f0c340.png">|





