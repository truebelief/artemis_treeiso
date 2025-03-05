
# Individual-tree isolation (treeiso) from terrestrial laser scanning

### **`CloudCompare Plugin (Developed May 2023)`**

A new upgrade in February 2025 boosts the processing speed by 5-15x. The plugin is now integrated into CloudCompare 2.14.alpha (Windows only).

### **`Python Version (Developed January 2025)`**

A new Python version accelerated by C++ has been developed.



### **Authors:**

Zhouxin Xi (zhouxin.xi@nrcan-rncan.gc.ca) and Chris Hopkinson (c.hopkinson@uleth.ca)

The University of Lethbridge - Department of Geography & Environment - Artemis Lab

Please cite:

> Xi, Z.; Hopkinson, C. 3D Graph-Based Individual-Tree Isolation (_Treeiso_) from Terrestrial Laser Scanning Point Clouds. _Remote Sens_. **2022**, _14_, 6116. https://doi.org/10.3390/rs14236116


This tool relies on the cut-pursuit algorithm, please also consider citing:

> Landrieu, L.; Obozinski, G. Cut Pursuit: Fast Algorithms to Learn Piecewise Constant Functions on General Weighted Graphs. SIAM J. Imaging Sci. 2017, 10, 1724–1766. [hal link](https://hal.archives-ouvertes.fr/hal-01306779)

> Raguet, H.; Landrieu, L. Cut-pursuit algorithm for regularizing nonsmooth functionals with graph total variation. In Proceedings of the International Conference on Machine Learning, Stockholm, Sweden, 10–15 July 2018; Volume 80, pp. 4247–4256.

## Folder Structure

```
├── data                                    # Raw data
│   └── LPine1_demo.laz                     # Example TLS data (subset due to file size limits)
├── Matlab                                  # treeiso Matlab source code
│   ├── treeiso.m                           # Main Matlab program
│   ├── cutPursuit.m                        # Wrapper for L0 cut-pursuit mex code
│   ├── jitknnsearch.m                      # JIT-accelerated knnsearch for small point sets
│   ├── overlapping.m                       # Calculates overlap ratio between crown convex hulls
│   ├── L0_cut_pursuit.mexw64               # Compiled L0 cut-pursuit executable
│   ├── L0_cut_pursuit_segmentation.mexw64  # Compiled L0 cut-pursuit segmentation executable
│   ├── las2mat.mexw64                      # Compiled LAS/LAZ file reader
│   └── mat2las.mexw64                      # Compiled LAS/LAZ file writer
├── Python                                  # treeiso Python source code
│   ├── treeiso.py                          # Main Python program
│   ├── cut_pursuit_L0.py                   # Simplified Python version of cut-pursuit (L0 norm only)
│   └── cut_pursuit_L0_replica_cpp.py       # Python version closer to the original C++ (slower; for reference)
├── PythonCPP                               # Python source code accelerated by C++ (Recommended)
│   ├── treeiso.py                          # Main Python program
│   └── cut_pursuit_L0.py                   # Backup simplified cut-pursuit (L0 norm only)
├── R                                       # R script
│   ├── treeiso_example.R                   # An R example
├── CloudCompare                            # CloudCompare plugin (external link)
├── cut_pursuit_py                          # C++ cut-pursuit Python binder (external link)
├── .gitmodules                             # External repository declarations
├── LICENSE
└── README.md
```

## 1. Preprocessing requirements

- **Noise Removal:** Clean your TLS data to remove noise (e.g., using CloudCompare with a Sober filter, NN=10, std=1).
- **Decimation:** Decimate the point cloud to a reasonable resolution (approximately 2 cm).
- **Ground Points:** Remove ground points before running tree isolation.

## 2. Matlab Version

1. **Execution:**  
   Run the Matlab script located at `Matlab/treeiso.m`.  
   When prompted, enter the folder path containing your LAS/LAZ files.  
   The script will search for all LAS/LAZ files in the specified directory and process each one.

2. **Testing:**  
   For a quick test, point the script to the data folder containing `LPine1_demo.laz`.

3. **Output:**  
   Processed LAS/LAZ files (with additional fields) will be saved in `Matlab/output/`. Each output file includes three additional fields:  
   - `init_segs`: Initial segmentation from 3D cut-pursuit  
   - `interim_segs`: Second-stage segmentation from 2D cut-pursuit  
   - `segs`: Final tree segmentation

4. **Performance:**  
   Processing a ~25 MB LAS/LAZ file takes about 10 minutes on an Intel Core i7-9700K with 16 GB RAM.

> **Note:**  
> The mex files were compiled on Windows 10 (64-bit) using Visual Studio 2019. If you need to compile them yourself, please refer to the following resources:  
> [cutpursuit](https://github.com/loicland/cut-pursuit) and [mastools](https://github.com/plitkey/matlas_tools).

## 3. Python Version

### Prerequisites

- Python 3.6–3.12 (newer versions are untested)

### 3.1 C++-Based Solution (Recommended) *(Updated March 2025)*

#### Installation

1. Open your terminal or command-line interface.
2. Clone or download the repository:
   ```bash
   git clone https://github.com/artemislab/treeiso.git
   cd treeiso/PythonCPP
   ```
3. Install the dependencies:
   ```bash
   pip install -r requirements.txt
   ```

#### Usage

Run the `treeiso.py` script in the `PythonCPP` folder:
```bash
python treeiso.py
```
When prompted, enter the path to the directory containing your LAS/LAZ or CSV files. The resulting LAZ or CSV file will be saved in the same directory with the suffix *_treeiso.laz or *_treeiso.csv.

> **Note:**  
> The program first checks for the presence of the C++-based `cut_pursuit_py` module. If it is not installed, it falls back to the pure Python version (`cut_pursuit_L0.py`). The C++ interface is provided via pre-built wheels for Ubuntu, macOS, and Windows, so you typically do not need a C++ compiler. The C++ version offers a speed-up of over 5× compared to the pure Python implementation, though the results may differ slightly.

### 3.2 Pure Python Solution

#### Installation

1. Clone or download the repository:
   ```bash
   git clone https://github.com/artemislab/treeiso.git
   cd treeiso/Python
   ```
2. Install the dependencies:
   ```bash
   pip install -r requirements.txt
   ```

#### Usage

Run the `treeiso.py` script in the `Python` folder:
```bash
python treeiso.py
```
When prompted, enter the path to your directory containing LAS/LAZ files. The output file will be saved in the same input directory.

## 4. R Version

An R script suggested by Gergo Dioszegi uses the CloudCompare command-line interface to automate processing. This script remains untested.


## Example isolated trees: 
|Raw TLS Example1| After treeiso isolation| Top view|
|:---:|:---:|:---:|
|<img width="635" alt="demo1_crop" src="https://user-images.githubusercontent.com/8785889/182312969-7c81949f-67fa-409b-bb24-73b094917c52.png">|<img width="635" alt="demo1_treeiso_crop" src="https://user-images.githubusercontent.com/8785889/182313130-7d6ba091-b2cb-4482-ae21-3ee2ec5bb34e.png">|<img width="700" alt="demo1_treeiso2_crop" src="https://user-images.githubusercontent.com/8785889/182313150-a14e7a3e-79b0-4d40-a400-ce89520e5e4d.png">|


|Raw TLS Example2| After treeiso isolation| Top view|
|:---:|:---:|:---:|
|<img width="1354" alt="TAspen1_crop" src="https://user-images.githubusercontent.com/8785889/182315836-78bdd338-2678-4d29-bcce-1442777e6918.png">|<img width="1354" alt="TAspen1_treeiso_crop" src="https://user-images.githubusercontent.com/8785889/182315950-f44937e2-0ec8-4b43-a980-0da9b7925534.png">|<img width="724" alt="TAspen1_treeiso2_crop" src="https://user-images.githubusercontent.com/8785889/182315979-d8e9c5ea-7a00-4490-8bab-c43fd6af9498.png">|


|Raw TLS Example3| After treeiso isolation| Top view|
|:---:|:---:|:---:|
|<img width="1228" alt="NPoplar2_crop" src="https://user-images.githubusercontent.com/8785889/182318128-ddeac093-f48d-4e69-a9b2-560f5a25335e.png">|<img width="1228" alt="NPoplar2_treeiso_crop" src="https://user-images.githubusercontent.com/8785889/182318163-290e0663-4e0a-455c-85d5-5160f00a8e33.png">|<img width="1026" alt="NPoplar2_treeiso2_crop" src="https://user-images.githubusercontent.com/8785889/182318398-fc428e17-1fd7-415c-9594-e0b9e3f0c340.png">|





