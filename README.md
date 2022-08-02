# Individual-tree isolation (treeiso) from terrestrial laser scanning

The University of Lethbridge - Department of Geography & Environment - Artemis Lab

Author - Zhouxin Xi (zhouxin.xi@uleth.ca) and Prof. Chris Hopkinson (c.hopkinson@uleth.ca)

## Folder structure

    ├── data                                    # All raw data
    │   ├── LPine1_demo.laz                     # Example TLS data (just a subset from a TLS plot due to the file size limit)
    ├── Matlab                                  # treeiso Matlab Source code 
    │   ├── treeiso.m                           # Main program
    │   ├── cutPursuit.m                        # Wrapper of L0 cut-pursuit mex code
    │   ├── jitknnsearch.m                      # Matlab jit accelerated knnsearch - fast knnsearch from small amount of points
    │   ├── overlapping.m                       # Used to calculate overlap ratio between a pair of crown convex hulls
    │   ├── L0_cut_pursuit.mexw64               # compiled L0 cut pursuit program
    │   ├── L0_cut_pursuit_segmentation.mexw64  # compiled L0 cut pursuit program
    │   ├── las2mat.mexw64                      # compiled las/laz file reading program
    │   ├── mat2las.mexw64                      # compiled las/laz file writing program
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
Run matlab code under Matlab\treeiso.m to isolate TLS trees

Type in the folder path where your laz files are located

The treeiso.m will exhaustively search all laz files and run the treeiso for each laz file

You can type in the data folder for a quick test. The treeiso.m will isolate trees from LPine1_demo.laz

Result laz file will be saved to Matlab\output\

The laz file will include three additional fields:

1. init_segs: show the initial segmentation from 3D cut-pursuit
2. interim_segs: show the 2nd-stage segmentation from 2D cut-pursuit
3. segs: final tree segmentation

Please be patient: processing a laz file of ~25MB costs about 10 minutes using Intel Core i7-9700K and 16GB RAM.

Example isolated trees:  
|Raw TLS Example1| After treeiso isolation| Top view|
|:---:|:---:|:---:|
|<img width="635" alt="demo1_crop" src="https://user-images.githubusercontent.com/8785889/182312969-7c81949f-67fa-409b-bb24-73b094917c52.png">|<img width="635" alt="demo1_treeiso_crop" src="https://user-images.githubusercontent.com/8785889/182313130-7d6ba091-b2cb-4482-ae21-3ee2ec5bb34e.png">|<img width="700" alt="demo1_treeiso2_crop" src="https://user-images.githubusercontent.com/8785889/182313150-a14e7a3e-79b0-4d40-a400-ce89520e5e4d.png">|


|Raw TLS Example2| After treeiso isolation| Top view|
|:---:|:---:|:---:|
|<img width="1354" alt="TAspen1_crop" src="https://user-images.githubusercontent.com/8785889/182315836-78bdd338-2678-4d29-bcce-1442777e6918.png">|<img width="1354" alt="TAspen1_treeiso_crop" src="https://user-images.githubusercontent.com/8785889/182315950-f44937e2-0ec8-4b43-a980-0da9b7925534.png">|<img width="724" alt="TAspen1_treeiso2_crop" src="https://user-images.githubusercontent.com/8785889/182315979-d8e9c5ea-7a00-4490-8bab-c43fd6af9498.png">|


|Raw TLS Example3| After treeiso isolation| Top view|
|:---:|:---:|:---:|
|<img width="1228" alt="NPoplar2_crop" src="https://user-images.githubusercontent.com/8785889/182318128-ddeac093-f48d-4e69-a9b2-560f5a25335e.png">|<img width="1228" alt="NPoplar2_treeiso_crop" src="https://user-images.githubusercontent.com/8785889/182318163-290e0663-4e0a-455c-85d5-5160f00a8e33.png">|<img width="1026" alt="NPoplar2_treeiso2_crop" src="https://user-images.githubusercontent.com/8785889/182318398-fc428e17-1fd7-415c-9594-e0b9e3f0c340.png">|





