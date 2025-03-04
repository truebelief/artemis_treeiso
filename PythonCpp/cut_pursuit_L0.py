import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import warnings
warnings.filterwarnings('ignore')

from scipy.spatial import cKDTree
import maxflow
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import time

@dataclass
class CPParameter:
    """Parameters for Cut Pursuit algorithm"""
    reg_strenth: float = 0.0
    flow_steps: int = 3
    max_ite_main: int = 6
    stopping_ratio: float = 0.0001


class CutPursuit:
    def __init__(self, n_vertices: int = 1, verbose: bool = True):
        self.verbose = verbose
        self.n_vertex = n_vertices
        self.dim = 1  # This will be updated in setup_cp
        self.parameter = CPParameter()

        # Graph structure
        self.n_total_vertices = n_vertices + 2  # including source and sink
        self.source = n_vertices
        self.sink = n_vertices + 1

        # Store vertex properties as NumPy arrays
        self.vertex_weights = np.zeros(self.n_total_vertices)
        # Don't initialize observation dimension yet, will do in setup_cp
        self.vertex_observations = None
        self.vertex_values = None
        self.vertex_colors = np.full(self.n_total_vertices, -1)
        self.vertex_components = np.zeros(self.n_total_vertices, dtype=int)

        # Edge structure using structured array
        self.edges_dtype = np.dtype([
            ('u', np.int32),
            ('v', np.int32),
            ('weight', np.float32),
            ('capacity', np.float32),
            ('is_active', np.bool_),
            ('real_edge', np.bool_)
        ])
        self.edges = np.zeros(0, dtype=self.edges_dtype)
        self.n_edge = 0

        # Component structure using arrays
        self.max_components = n_vertices
        self.component_indices = [[] for _ in range(self.max_components)]
        self.root_vertex = np.zeros(self.max_components, dtype=np.int32)
        self.saturated_vertices_logi = np.zeros(n_vertices, dtype=bool)
        self.n_active_components = 1
        self.saturated_components = np.zeros(1, dtype=bool)

        # Initialize source and sink
        self.vertex_weights[self.source] = 1.0
        self.vertex_weights[self.sink] = 1.0

        # Cache for optimization
        self.flow_graph = maxflow.Graph[int]()


    def set_parameters(self, flow_steps=4, max_ite_main=20, stopping_ratio=0.001, reg_strenth=0.0):
        """Set parameters for the Cut Pursuit algorithm"""
        self.parameter.flow_steps = flow_steps
        self.parameter.max_ite_main = max_ite_main
        self.parameter.stopping_ratio = stopping_ratio
        self.parameter.reg_strenth = reg_strenth

    def initialize(self):
        # Initialize first component efficiently using sparse matrix
        self.component_indices[0] = list(range(self.n_vertex))
        self.root_vertex[0] = 0
        self.vertex_components[:self.n_vertex] = 0


        # Rest of initialization remains the same...
        self.compute_value(0)

        # Create source/sink edges efficiently
        n_source_sink_edges = 4 * self.n_vertex
        source_sink_edges = np.zeros(n_source_sink_edges, dtype=self.edges_dtype)

        vertices = np.arange(self.n_vertex)
        idx = np.arange(0, n_source_sink_edges, 4)

        source_sink_edges['u'][idx] = self.source
        source_sink_edges['v'][idx] = vertices
        source_sink_edges['u'][idx + 1] = vertices
        source_sink_edges['v'][idx + 1] = self.source
        source_sink_edges['u'][idx + 2] = vertices
        source_sink_edges['v'][idx + 2] = self.sink
        source_sink_edges['u'][idx + 3] = self.sink
        source_sink_edges['v'][idx + 3] = vertices

        if len(self.edges) > 0:
            self.edges = np.concatenate([self.edges, source_sink_edges])
        else:
            self.edges = source_sink_edges
        self.n_edge = len(self.edges)

    def run(self) -> Tuple[List[float], List[float]]:
        """Optimized main optimization loop"""
        self.initialize()
        if self.verbose:
            print(f"Graph with {self.n_vertex} vertices and {len(self.edges)} edges "
                  f"and observation of dimension {self.dim}")

        # Initial energy computation
        fidelity_energy, penalty_energy = self.compute_energy()
        energy_zero = fidelity_energy
        old_energy = fidelity_energy + penalty_energy

        # Pre-allocate arrays for benchmarking
        energy_out = np.zeros(self.parameter.max_ite_main)
        time_out = np.zeros(self.parameter.max_ite_main)


        start_time = time.time()

        # Main optimization loop
        for ite_main in range(self.parameter.max_ite_main):
            # Compute optimal binary partition
            saturation = self.split()

            # Reduce graph
            self.reduce()

            # Compute new energy
            fidelity_energy, penalty_energy = self.compute_energy()
            current_total_energy = fidelity_energy + penalty_energy

            # Store benchmarking data
            energy_out[ite_main] = current_total_energy
            time_out[ite_main] = time.time() - start_time
            if self.verbose:
                print(f"Iteration {ite_main + 1:3} - {self.n_active_components:4} components - "
                      f"Saturation {100.0 * saturation / self.n_vertex:5.1f}% - "
                      f"Quadratic Energy {100 * fidelity_energy / energy_zero:6.3f}% - "
                      f"Timer {time_out[ite_main]:.2f}s")

            # Check stopping criteria
            if saturation == self.n_vertex:
                if self.verbose:
                    print("All components are saturated")
                break

            if ((old_energy - current_total_energy) / old_energy < self.parameter.stopping_ratio):
                if self.verbose:
                    print("Stopping criterion reached")
                break
            old_energy = current_total_energy

        return (
            energy_out[:ite_main + 1].tolist(),
            time_out[:ite_main + 1].tolist()
        )

    def compute_value(self, ind_com: int):
        comp_vertices = self.component_indices[ind_com]
        weights = self.vertex_weights[comp_vertices]
        total_weight = np.sum(weights)

        if total_weight > 0:
            comp_value = np.sum(weights[:, np.newaxis] * self.vertex_observations[comp_vertices], axis=0) / total_weight
        else:
            comp_value = np.zeros(self.dim)

        self.vertex_values[comp_vertices] = comp_value  # This updates the array in place
        self.vertex_components[comp_vertices] = ind_com

        return comp_value, total_weight



    def compute_energy(self) -> Tuple[float, float]:
        """Optimized energy computation using vectorized operations"""
        mask = self.vertex_weights[:self.n_vertex] > 0

        # Compute differences vectorized, only for actual vertices
        diff = self.vertex_observations[:self.n_vertex][mask] - self.vertex_values[:self.n_vertex][mask]
        weights = self.vertex_weights[:self.n_vertex][mask]

        # Compute energy using dot product
        fidelity_energy = 0.5 * np.sum(weights * np.sum(diff * diff, axis=1))

        # Regularization energy (penalty term)
        active_edges = self.edges[self.edges['is_active'] & self.edges['real_edge']]
        penalty_energy = 0.5 * self.parameter.reg_strenth * np.sum(active_edges['weight'])

        return fidelity_energy, penalty_energy

    def split(self) -> int:
        """Optimized split computation"""
        binary_label = np.zeros(self.n_vertex, dtype=bool)
        self.init_labels(binary_label)

        # Pre-allocate centers array
        centers = np.zeros((self.n_active_components, 2, self.dim))

        self.edge_mask=~self.edges['is_active'] & self.edges['real_edge']
        self.real_edges=self.edges[self.edge_mask]

        # Flow approximation loop
        for i_step in range(self.parameter.flow_steps):
            self.compute_centers(centers, binary_label)
            self.set_capacities(centers)
            binary_label.fill(False)
            source_idx=self.compute_max_flow()
            binary_label[source_idx] = True

        self.vertex_colors.fill(False)
        self.vertex_colors[source_idx]=True
        return self.activate_edges()

    def init_labels(self, binary_label: np.ndarray):
        """Initialize labels using Quickshift"""
        active_comps = np.where(~self.saturated_components)[0]

        for ind_com in active_comps:
            comp_vertices = self.component_indices[ind_com]
            observations = self.vertex_observations[comp_vertices]
            variances = np.var(observations, axis=0, ddof=0)
            var_dim = np.argmax(variances)
            median_value = np.median(observations[:, var_dim])
            binary_label[comp_vertices] = observations[:, var_dim] > median_value

    def compute_centers(self, centers: np.ndarray, binary_label: np.ndarray):
        # Get indices of unsaturated components
        active_comps = np.where(~self.saturated_components)[0]

        # Initialize arrays to store weights and observations
        total_weights_label0 = np.zeros(len(active_comps))
        total_weights_label1 = np.zeros(len(active_comps))
        sum_obs_label0 = np.zeros((len(active_comps), self.dim))
        sum_obs_label1 = np.zeros((len(active_comps), self.dim))

        for idx, ind_com in enumerate(active_comps):
            comp_vertices = self.component_indices[ind_com]
            weights = self.vertex_weights[comp_vertices]
            observations = self.vertex_observations[comp_vertices]
            labels = binary_label[comp_vertices]

            # Compute sums and weights for both labels
            weights_label0 = weights[~labels]
            weights_label1 = weights[labels]

            if len(weights_label0) == 0 or len(weights_label1) == 0:
                self.saturated_components[ind_com] = True
                centers[ind_com] = self.vertex_values[comp_vertices[0]]
                continue

            obs_label0 = observations[~labels]
            obs_label1 = observations[labels]

            total_weights_label0[idx] = np.sum(weights_label0)
            total_weights_label1[idx] = np.sum(weights_label1)

            sum_obs_label0[idx] = np.sum(obs_label0 * weights_label0[:, np.newaxis], axis=0)
            sum_obs_label1[idx] = np.sum(obs_label1 * weights_label1[:, np.newaxis], axis=0)

        # Compute centers for all components at once
        centers[active_comps, 0] = sum_obs_label0 / total_weights_label0[:, np.newaxis]
        centers[active_comps, 1] = sum_obs_label1 / total_weights_label1[:, np.newaxis]

    def set_capacities(self, centers: np.ndarray):
        """Optimized capacity setting using batch processing and vectorized operations"""
        SCALE_FACTOR = 1000

        # Initialize a new graph for this iteration
        self.flow_graph.reset()
        # Add nodes all at once
        node_ids = self.flow_graph.add_nodes(self.n_vertex)

        # Initialize arrays for terminal capacities
        source_caps = np.zeros(self.n_vertex, dtype=np.int64)
        sink_caps = np.zeros(self.n_vertex, dtype=np.int64)

        # Process each component
        for ind_com in range(self.n_active_components):
            if self.saturated_components[ind_com]:
                continue

            comp_vertices = np.array(self.component_indices[ind_com], dtype=np.int32)
            if len(comp_vertices) == 0:
                continue

            # Compute costs
            obs_diff0 = self.vertex_observations[comp_vertices] - centers[ind_com, 0]
            obs_diff1 = self.vertex_observations[comp_vertices] - centers[ind_com, 1]
            vertex_weights = self.vertex_weights[comp_vertices]

            cost_B = 0.5 * np.sum(obs_diff0 ** 2, axis=1) * vertex_weights
            cost_notB = 0.5 * np.sum(obs_diff1 ** 2, axis=1) * vertex_weights

            # Determine capacities for source and sink
            mask_to_sink = cost_B <= cost_notB
            cost_diff = np.abs(cost_B - cost_notB) * SCALE_FACTOR
            cost_diff = cost_diff.astype(np.int64)

            # Assign capacities
            source_caps[comp_vertices] = np.where(~mask_to_sink, cost_diff, 0)
            sink_caps[comp_vertices] = np.where(mask_to_sink, cost_diff, 0)

        # Add terminal edges in batch
        self.flow_graph.add_grid_tedges(node_ids, source_caps, sink_caps)

        # Add regular edges efficiently
        if len(self.real_edges) > 0:
            # Pre-compute all edge capacities
            edge_caps = (self.real_edges['weight'] * self.parameter.reg_strenth * SCALE_FACTOR).astype(np.int64)
            # Add edges with forward and backward capacities
            self.flow_graph.add_edges(
                self.real_edges['u'].astype(np.int32),  # source vertices
                self.real_edges['v'].astype(np.int32),  # target vertices
                edge_caps,  # forward capacities
                edge_caps  # backward capacities
            )

    def compute_max_flow(self):
        """Updated maximum flow computation using pymaxflow"""
        flow = self.flow_graph.maxflow()
        reachable = np.where(self.flow_graph.get_grid_segments(np.arange(self.n_vertex)))[0]
        return reachable

    def activate_edges(self) -> int:
        # Compute saturation directly
        saturation = np.sum([len(self.component_indices[i]) for i in range(self.n_active_components) if self.saturated_components[i]])

        # Find crossing edges efficiently
        edges_mask = self.edges['real_edge']
        u_colors = self.vertex_colors[self.edges['u'][edges_mask]]
        v_colors = self.vertex_colors[self.edges['v'][edges_mask]]
        crossing_edges = u_colors != v_colors
        crossing_indices = np.where(edges_mask)[0][crossing_edges]

        # Activate edges
        self.edges['is_active'][crossing_indices] = True
        self.edge_mask=~self.edges['is_active'] & self.edges['real_edge']
        self.real_edges=self.edges[self.edge_mask]
        return saturation

    def reduce(self):
        """Compute reduced graph and perform backward step if needed"""
        self.compute_connected_components()
        n_comp = self.n_active_components
        for ind_com in range(n_comp):
            self.compute_value(ind_com)

    def compute_connected_components(self):
        """Optimized connected components computation with efficient grouping"""
        if len(self.real_edges) == 0:
            # Create single arange array and reuse it
            vertex_indices = np.arange(self.n_vertex, dtype=np.int32)
            self.n_active_components = self.n_vertex
            self.component_indices = [[i] for i in vertex_indices]  # Still needs list comprehension for structure
            self.root_vertex = vertex_indices.copy()
            self.vertex_components[:self.n_vertex] = vertex_indices
            self.saturated_components = np.zeros(self.n_vertex, dtype=bool)
            return

        graph = csr_matrix(
            (np.ones(len(self.real_edges), dtype=bool), (self.real_edges['u'], self.real_edges['v'])),
            shape=(self.n_vertex, self.n_vertex),
            dtype=bool
        )

        # Compute connected components
        n_components, labels = connected_components(
            graph, directed=False, return_labels=True
        )

        # Initialize new component structures
        self.n_active_components = n_components

        # Use argsort to group vertices by component label efficiently
        sort_idx = np.argsort(labels)
        sorted_labels = labels[sort_idx]

        # Find the boundaries between different components
        boundaries = np.nonzero(np.diff(sorted_labels))[0] + 1
        boundaries = np.concatenate([[0], boundaries, [len(labels)]])

        # Create component indices using the boundaries
        self.component_indices = np.split(sort_idx, boundaries[1:-1])

        # Set root vertices (first vertex of each component)
        self.root_vertex = np.array([indices[0] if len(indices) else 0 for indices in self.component_indices], dtype=np.int32)

        # Update vertex components directly
        self.vertex_components[:self.n_vertex] = labels

        # Reset saturation status
        self.saturated_components = np.zeros(n_components, dtype=bool)



def setup_cp(n_nodes: int, n_edges: int, n_obs: int,
             observation: np.ndarray, eu: np.ndarray, ev: np.ndarray,
             edge_weight: np.ndarray, node_weight: np.ndarray, verbose:bool) -> CutPursuit:
    """Optimized setup for Cut Pursuit"""
    cp = CutPursuit(n_nodes,verbose=verbose)
    cp.dim = n_obs

    # Initialize observation arrays with correct dimensions
    cp.vertex_observations = np.zeros((cp.n_total_vertices, n_obs))
    cp.vertex_values = np.zeros((cp.n_total_vertices, n_obs))

    # Initialize vertex properties efficiently
    cp.vertex_weights[:n_nodes] = node_weight
    cp.vertex_observations[:n_nodes] = observation

    # Create edges array efficiently
    edges = np.zeros(2 * n_edges, dtype=cp.edges_dtype)

    # Forward edges
    edges['u'][:n_edges] = eu
    edges['v'][:n_edges] = ev
    edges['weight'][:n_edges] = edge_weight
    edges['capacity'][:n_edges] = edge_weight
    edges['real_edge'][:n_edges] = True

    # Reverse edges
    edges['u'][n_edges:] = ev
    edges['v'][n_edges:] = eu
    edges['weight'][n_edges:] = edge_weight
    edges['capacity'][n_edges:] = edge_weight
    edges['real_edge'][n_edges:] = True

    cp.edges = edges
    cp.n_edge = len(edges)

    return cp

# def cut_pursuit(n_nodes, n_edges, n_obs, y, eu, ev, edge_weight, node_weight, lambda_):
def cut_pursuit(n_nodes: int, n_edges: int, n_obs: int,
                observation: np.ndarray, eu: np.ndarray, ev: np.ndarray,
                edge_weight: np.ndarray, node_weight: np.ndarray,
                lambda_: float,verbose:bool) -> Tuple[np.ndarray, List[List[int]], np.ndarray, np.ndarray, np.ndarray]:
    """Main cut pursuit function with optimized setup and execution"""
    # Set random seed for reproducibility
    # np.random.seed(1)

    # print("L0-CUT PURSUIT")

    # Setup and run cut pursuit
    cp = setup_cp(n_nodes, n_edges, n_obs, observation, eu, ev, edge_weight, node_weight, verbose)

    # Set parameters
    cp.parameter.flow_steps = 4
    cp.parameter.max_ite_main = 20
    cp.parameter.stopping_ratio = 0.001
    cp.parameter.reg_strenth = lambda_

    # Run optimization
    energy_out, time_out = cp.run()

    return (
        cp.vertex_values[:n_nodes].copy(),  # Ensure we return a copy
        cp.component_indices[:cp.n_active_components],
        cp.vertex_components[:n_nodes].copy(),
        np.array(energy_out),
        np.array(time_out)
    )

def perform_cut_pursuit(reg_strength,D,pc_vec,edge_weights,Eu,Ev,verbose,progress_callback=None):
    solution, components, in_component, energy_out, time_out = cut_pursuit(
        n_nodes=len(pc_vec),
        n_edges=len(Eu),
        n_obs=D,
        observation=pc_vec,
        eu=Eu,
        ev=Ev,
        edge_weight=edge_weights,
        node_weight=np.ones(len(pc_vec)),
        lambda_=reg_strength,
        verbose=verbose
    )

    return in_component


def decimate_pcd(columns,min_res):
    _, block_idx_uidx, block_inverse_idx = np.unique(np.floor(columns[:,:3]/min_res).astype(np.int32),axis=0, return_index=True, return_inverse=True)
    return block_idx_uidx,block_inverse_idx

if __name__ == "__main__":
    import os
    import laspy

    import cProfile
    import pstats
    import io


    K = 5
    regStrength = 1.0
    min_res=0.05

    path_to_las = r"..\data\LPine1_demo.laz"
    pcd_basename = os.path.basename(path_to_las)[:-4]
    las = laspy.read(path_to_las)
    pcd=np.transpose([las.x,las.y,las.z])
    pcd=pcd-np.mean(pcd[:,:3],axis=0)

    dec_idx_uidx, dec_inverse_idx = decimate_pcd(pcd[:, :3], min_res)  # reduce points first
    pcd_dec = pcd[dec_idx_uidx]

    # Create profiler
    pr = cProfile.Profile()
    pr.enable()

    # Time the decimation
    t0 = time.time()

    point_count=len(pcd_dec)
    kdtree = cKDTree(pcd_dec[:,:3])
    nn_D, nn_idx = kdtree.query(pcd_dec, k=K + 1)

    # Remove self-connections
    # distances = nn_D[:, 1:]
    indices = nn_idx[:, 1:]

    # Create edge list
    n_nodes = len(pcd_dec)
    n_obs=3
    n_edges=n_nodes*K

    eu = np.repeat(np.arange(n_nodes), K)
    ev = indices.ravel()

    y=pcd_dec[:,:3]-np.mean(pcd_dec[:,:3],axis=0)

    # Edge weights following C++ implementation
    edge_weight = np.ones_like(eu)
    node_weight = np.ones(point_count)

    solution, components, in_component, energy_out, time_out = cut_pursuit(
        n_nodes=n_nodes,
        n_edges=n_edges,
        n_obs=n_obs,
        observation=y,
        eu=eu,
        ev=ev,
        edge_weight=edge_weight,
        node_weight=node_weight,
        lambda_=regStrength,
        verbose=False
    )


    main_algo_time = time.time() - t0
    pr.disable()

    # Print detailed stats
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
    ps.print_stats(30)  # Print top 30 time-consuming functions

    print("\n=== Overall Timing ===")
    print(f"Main algorithm time: {main_algo_time:.2f} seconds")
    print("\n=== Detailed Function Profiling ===")
    print(s.getvalue())

    las.add_extra_dim(laspy.ExtraBytesParams(name="cutpursuit", type="int32", description="cutpursuit"))
    las.cutpursuit = in_component[dec_inverse_idx]
    las.write(r"../output/LPine1_demo_result.laz")


