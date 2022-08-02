# artemis_treeiso

Individual-tree isolation (treeiso) from terrestrial laser scanning point clouds

The University of Lethbridge - Department of Geography & Environment - Artemis Lab

Author - Zhouxin Xi (zhouxin.xi@uleth.ca) and Prof. Chris Hopkinson (c.hopkinson@uleth.ca)

## Folder structure

    ├── data                                    # All raw data
    │   ├── LPine_plot1_reference_demo.laz      # Example TLS data
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

### 2. treeiso
Run matlab code under Matlab\treeiso.m to isolate TLS trees

Result laz file will be saved to Matlab\output\ folder

The laz file will include three additional fields:

init_segs (show the initial segmentation from cut-pursuit, 

Example isolated trees from a plot:  
|Example 1| Example 2|
|:---:|:---:|
|<img width="407" alt="tree_itc4" src="https://user-images.githubusercontent.com/8785889/153339986-63e9495b-4951-4252-a089-803e50dcd0b6.png">|<img width="513" alt="tree_itc5" src="https://user-images.githubusercontent.com/8785889/153341308-42afce5f-f8ea-4179-b0b1-3aef60c176f3.png">|
