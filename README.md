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

Ground points must be removed before tree isolation

The mex files were all compiled under Windows 10 (64bit) Intel using Visual Studio 2019. If you need to compile on your own, please refer to the authors' sites for detailed steps.

Links below:

[cutpursuit](https://github.com/loicland/cut-pursuit)

[mastools](https://github.com/plitkey/matlas_tools)

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

Please be patient: processing a laz file of ~25MB costs about 10 minutes using Intel(R)i7-9700K@3.60GHz and 16GB RAM.

Example isolated trees from a plot:  
|Raw TLS Example1| After treeiso isolation|
|:---:|:---:|
|<img width="407" alt="tree_itc4" src="https://user-images.githubusercontent.com/8785889/153339986-63e9495b-4951-4252-a089-803e50dcd0b6.png">|<img width="513" alt="tree_itc5" src="https://user-images.githubusercontent.com/8785889/153341308-42afce5f-f8ea-4179-b0b1-3aef60c176f3.png">|

|Raw TLS Example2| After treeiso isolation|
|:---:|:---:|
|<img width="407" alt="tree_itc4" src="https://user-images.githubusercontent.com/8785889/153339986-63e9495b-4951-4252-a089-803e50dcd0b6.png">|<img width="513" alt="tree_itc5" src="https://user-images.githubusercontent.com/8785889/153341308-42afce5f-f8ea-4179-b0b1-3aef60c176f3.png">|

|Raw TLS Example3| After treeiso isolation|
|:---:|:---:|
|<img width="407" alt="tree_itc4" src="https://user-images.githubusercontent.com/8785889/153339986-63e9495b-4951-4252-a089-803e50dcd0b6.png">|<img width="513" alt="tree_itc5" src="https://user-images.githubusercontent.com/8785889/153341308-42afce5f-f8ea-4179-b0b1-3aef60c176f3.png">|


