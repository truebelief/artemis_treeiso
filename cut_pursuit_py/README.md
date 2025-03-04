# Cut Pursuit Segmentation

## Overview

Cut Pursuit is an efficient algorithm for segmenting point clouds by minimizing a functional over a graph. This package provides a Python interface to the C++ implementation of the Cut Pursuit algorithm.

## Installation

```bash
pip install cut-pursuit
```

## Usage Example

```python
import numpy as np
import cut_pursuit

# Assume pcd is a numpy array of 3D points (N x 3)
def segment_point_cloud(pcd, k=7, reg_strength=1.0):
    # Preprocess point cloud 
    pcd = pcd - np.mean(pcd, axis=0)
    
    # Compute k-nearest neighbors
    from scipy.spatial import cKDTree
    kdtree = cKDTree(pcd)
    _, nn_idx = kdtree.query(pcd, k=k)
    
    # Prepare graph structure
    indices = nn_idx[:, 1:]  # exclude self
    n_nodes = len(pcd)
    
    # Create edge lists
    eu = np.repeat(np.arange(n_nodes), k-1)
    ev = indices.ravel()
    
    # Edge weights 
    edge_weights = np.ones_like(eu, dtype=np.float32)
    
    # Perform cut pursuit
    segments = cut_pursuit.perform_cut_pursuit(
        K=k,              # Number of neighbors
        reg_strength=reg_strength,  # Regularization strength
        D=3,              # Dimension of points
        pc_vec=pcd,        # Point cloud 
        edge_weights=edge_weights,
        first_edge=np.cumsum(np.repeat(k-1, n_nodes+1))[:-1],
        adj_vertices=ev
    )
    
    return segments

# Example usage
point_cloud = np.random.rand(1000, 3)
segmentation = segment_point_cloud(point_cloud)
print(f"Number of segments: {len(np.unique(segmentation))}")
```

## Dependencies

- NumPy
- C++11 compatible compiler

## Citation

If you use this implementation, please cite the original paper:

Landrieu, L., & Obozinski, G. (2017). Cut Pursuit: Fast Algorithms to Learn Piecewise Constant Functions on General Weighted Graphs. SIAM Journal on Imaging Sciences, 10(4), 1724-1766.

## License

[MIT]
