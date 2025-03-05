import numpy as np
import networkx as nx
from networkx.algorithms.flow import boykov_kolmogorov
import warnings
import random

import heapq
from scipy.spatial import cKDTree

from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import enum
import time
class FidelityType(enum.Enum):
    L2 = 0
    LINEAR = 1
    KL = 2
    SPG = 3

warnings.filterwarnings('ignore')

@dataclass
class CPParameter:
    """Parameters for Cut Pursuit algorithm"""
    reg_strenth: float = 0.0  # regularization strength, multiply the edge weight
    cutoff: int = 0  # minimal component size
    flow_steps: int = 3  # number of steps in the optimal binary cut computation
    kmeans_ite: int = 5  # number of iteration in the kmeans sampling
    kmeans_resampling: int = 3  # number of kmeans re-initialization
    verbose: int = 2  # verbosity
    # max_ite_main: int = 2  # max number of iterations in the main loop
    max_ite_main: int = 6  # max number of iterations in the main loop
    backward_step: bool = True  # indicates if a backward step should be performed
    stopping_ratio: float = 0.0001  # when (E(t-1) - E(t) / (E(0) - E(t)) is too small, algorithm stops
    fidelity: FidelityType = FidelityType.L2  # the fidelity function
    smoothing: float = 0.1  # smoothing term (for KL divergence only)
    parallel: bool = True  # enable/disable parallelism
    weight_decay: float = 0.7  # for continued optimization of the flow steps


class CutPursuit:
    def __init__(self, n_vertices: int = 1):
        # Initialize main graph as NetworkX DiGraph
        self.main_graph = nx.DiGraph()
        self.components: List[List[int]] = [[]]
        self.root_vertex: List[int] = [0]
        self.saturated_components: List[bool] = [False]
        self.source: int = None
        self.sink: int = None
        self.dim: int = 1
        self.n_vertex: int = n_vertices
        self.n_edge: int = 0
        self.parameter = CPParameter()

    def set_parameters(self, flow_steps=4, weight_decay=0.7, kmeans_ite=8,
                       kmeans_resampling=5, max_ite_main=20, backward_step=True,
                       stopping_ratio=0.001, reg_strenth=0.0, cutoff=0):
        """Set parameters for the Cut Pursuit algorithm"""
        self.parameter.flow_steps = flow_steps
        self.parameter.weight_decay = weight_decay
        self.parameter.kmeans_ite = kmeans_ite
        self.parameter.kmeans_resampling = kmeans_resampling
        self.parameter.max_ite_main = max_ite_main
        self.parameter.backward_step = backward_step
        self.parameter.stopping_ratio = stopping_ratio
        self.parameter.reg_strenth = reg_strenth
        self.parameter.cutoff = cutoff

    def initialize(self):
        """Initialize the cut pursuit structure"""
        # First component contains all vertices (excluding source/sink)
        self.components[0] = list(range(self.n_vertex))
        self.root_vertex[0] = 0

        # Compute the value for the initial component
        self.compute_value(0)

        # Note: source and sink nodes already exist from setup_cp
        # Just need to add edges from source to vertices and vertices to sink
        e_index = self.n_edge
        for vertex_idx in range(self.n_vertex):
            # Add edge from source to vertex
            self.main_graph.add_edge(self.source, vertex_idx,
                                     index=e_index,
                                     weight=0.0,
                                     capacity=0.0,
                                     residual_capacity=0.0,
                                     is_active=False,
                                     real_edge=False)
            self.main_graph.add_edge(vertex_idx, self.source,
                                     index=e_index + 1,
                                     weight=0.0,
                                     capacity=0.0,
                                     residual_capacity=0.0,
                                     is_active=False,
                                     real_edge=False)
            e_index += 2

            # Add edge from vertex to sink
            self.main_graph.add_edge(vertex_idx, self.sink,
                                     index=e_index,
                                     weight=0.0,
                                     capacity=0.0,
                                     residual_capacity=0.0,
                                     is_active=False,
                                     real_edge=False)
            self.main_graph.add_edge(self.sink, vertex_idx,
                                     index=e_index + 1,
                                     weight=0.0,
                                     capacity=0.0,
                                     residual_capacity=0.0,
                                     is_active=False,
                                     real_edge=False)
            e_index += 2

        ii=1

    def run(self) -> Tuple[List[float], List[float]]:
        """Main optimization loop"""
        # First initialize the structure
        self.initialize()

        if self.parameter.verbose > 0:
            print(f"Graph {len(self.main_graph.nodes)} vertices and {len(self.main_graph.edges)} edges and observation of dimension {self.dim}")

        # Compute initial energy (with 1 component)
        energy_pair = self.compute_energy()
        energy_zero = energy_pair[0]  # Energy with 1 component
        old_energy = energy_zero  # Energy at the previous iteration

        # Initialize vectors for time and energy benchmarking
        energy_out = [0.0] * self.parameter.max_ite_main
        time_out = [0.0] * self.parameter.max_ite_main
        start_time = time.time()

        # Main loop
        for ite_main in range(1, self.parameter.max_ite_main + 1):
            # Compute optimal binary partition
            saturation = self.split()

            # Compute the new reduced graph
            self.reduce()

            # Compute energy after the iteration
            energy = self.compute_energy()
            current_total_energy = energy[0] + energy[1]
            energy_out[ite_main - 1] = current_total_energy
            time_out[ite_main - 1] = time.time() - start_time

            if self.parameter.verbose > 1:
                print(f"Iteration {ite_main:3} - {len(self.components):4} components - ", end="")
                print(f"Saturation {100.0 * saturation / self.n_vertex:5.1f} % - ", end="")
                print(f"Quadratic Energy {100 * current_total_energy / energy_zero:4.3f} % - ", end="")
                print(f"Timer {time.time() - start_time}")

            # Stopping checks
            if saturation == self.n_vertex:
                # All components are saturated
                if self.parameter.verbose > 1:
                    print("All components are saturated")
                break

            # Check relative energy progress
            if ((old_energy - current_total_energy) / old_energy <
                    self.parameter.stopping_ratio):
                # Relative energy progress stopping criterion
                if self.parameter.verbose > 1:
                    print("Stopping criterion reached")
                break

            if ite_main >= self.parameter.max_ite_main:
                # Max number of iterations
                if self.parameter.verbose > 1:
                    print("Max number of iteration reached")
                break

            old_energy = current_total_energy

        # Apply cutoff if specified
        if self.parameter.cutoff > 0:
            self.cutoff()

        # Return energy and time vectors
        # Trim arrays to actual number of iterations
        energy_out = energy_out[:ite_main]
        time_out = time_out[:ite_main]

        return energy_out, time_out

    def compute_value(self, ind_com: int) -> Tuple[np.ndarray, float]:
        """Compute the value and total weight of a component"""
        total_weight = 0.0
        comp_value = np.zeros(self.dim)

        # For each vertex in the component
        for vertex_idx in self.components[ind_com]:
            node = self.main_graph.nodes[vertex_idx]
            weight = node['weight']
            total_weight += weight

            # Update component value (weighted average)
            comp_value += node['observation'] * weight

            # Update vertex's component index
            node['in_component'] = ind_com

        # Normalize by total weight
        if total_weight > 0:
            comp_value /= total_weight

        # Set the computed value for all vertices in the component
        for vertex_idx in self.components[ind_com]:
            self.main_graph.nodes[vertex_idx]['value'] = comp_value

        return comp_value, total_weight

    def compute_energy(self) -> Tuple[float, float]:
        """
        Compute the energy of the current solution.
        Returns:
            Tuple[float, float]: (fidelity_energy, penalty_energy)
        """
        # First compute the fidelity term (quadratic energy)
        fidelity_energy = 0.0

        # For each vertex, compute the fidelity energy
        for vertex_idx in range(self.n_vertex):
            node = self.main_graph.nodes[vertex_idx]
            # Skip if vertex has no weight
            if node['weight'] == 0:
                continue

            # Compute squared difference between observation and current value
            diff = node['observation'] - node['value']
            fidelity_energy += 0.5 * node['weight'] * np.sum(diff * diff)

        # Then compute the penalty term (regularization)
        penalty_energy = 0.0

        # For each edge, compute the penalty energy
        for u, v, edge_data in self.main_graph.edges(data=True):
            # Skip non-real edges (connections to source/sink)
            if not edge_data['real_edge']:
                continue

            # Add penalty if edge is active
            penalty_energy += 0.5 * edge_data['is_active'] * self.parameter.reg_strenth * edge_data['weight']

        return fidelity_energy, penalty_energy


    def split(self) -> int:
        """
        Compute the optimal binary partition.
        Returns:
            int: Number of saturated vertices
        """
        # Load number of components
        nb_comp = len(self.components)

        # Stores whether each vertex is in B or not
        binary_label = [False] * self.n_vertex

        # Initialize the binary partition with kmeans
        self.init_labels(binary_label)
        # binary_label=np.loadtxt("F:\\prj\\CC2\\comp\\TreeAIBox\\test\\testCutPursuit\\init_labels.txt")
        # Centers is the value of each binary component in the optimal partition
        # Shape: (nb_comp, 2, dim) - for each component, 2 centers (B and not B), each of dimension dim
        centers = np.zeros((nb_comp, 2, self.dim))

        # Main loop - the optimal flow is iteratively approximated
        for i_step in range(self.parameter.flow_steps):
            # Compute centers (h_1 and h_2)
            centers.fill(0)  # Reset centers
            self.compute_centers(centers, nb_comp, binary_label)

            # Update the capacities of the flow graph
            self.set_capacities(centers)

            # self.export_graph_state("graph_before_maxflow_pred")

            # Compute max flow
            self.compute_max_flow()

            # self.export_graph_state("graph_after_maxflow_pred")

            # Update binary labels based on sink/source association
            for ind_com in range(nb_comp):
                if self.saturated_components[ind_com]:
                    continue

                for vertex_idx in self.components[ind_com]:
                    # A vertex belongs to B if it's connected to sink after max flow
                    binary_label[vertex_idx] = (self.main_graph.nodes[vertex_idx]['color'] ==
                                                self.main_graph.nodes[self.sink]['color'])
            # binary_label1 = np.loadtxt(r"F:\prj\CC2\comp\TreeAIBox\test\testCutPursuit\binary_label_truth.txt")
            # ii=1
        # Activate edges based on the partition
        saturation = self.activate_edges()

        return saturation

    def init_labels(self, binary_label: List[bool]):
        """Initialize the labeling for each component with kmeans"""
        # For each component
        for ind_com in range(len(self.components)):
            if self.saturated_components[ind_com] or len(self.components[ind_com]) <= 1:
                continue

            comp_size = len(self.components[ind_com])
            best_energy = float('inf')
            potential_label = [False] * comp_size
            energy_array = np.zeros(comp_size)

            # Try multiple kmeans initializations
            for init_kmeans in range(self.parameter.kmeans_resampling):
                # Initialize kernels - shape: (2, dim) for 2 centers
                kernels = np.zeros((2, self.dim))
                total_weight = [0.0, 0.0]

                # KM++ initialization
                # First kernel is random
                first_kernel = random.randrange(comp_size)
                first_vertex = self.components[ind_com][first_kernel]
                kernels[0] = self.main_graph.nodes[first_vertex]['observation']

                # Compute distances to first kernel for second kernel selection
                current_energy = 0.0

                for i_ver in range(comp_size):
                    vertex = self.components[ind_com][i_ver]
                    vertex_data = self.main_graph.nodes[vertex]
                    diff = vertex_data['observation'] - kernels[0]
                    energy_array[i_ver] = np.sum(diff * diff) * vertex_data['weight']
                    current_energy += energy_array[i_ver]

                # Select second kernel proportional to distance
                random_sample = random.random() * current_energy
                second_kernel = comp_size - 1  # default value
                for i_ver in range(comp_size):
                    random_sample -= energy_array[i_ver]
                    if random_sample < 0:
                        second_kernel = i_ver
                        break

                second_vertex = self.components[ind_com][second_kernel]
                kernels[1] = self.main_graph.nodes[second_vertex]['observation']

                # Main kmeans loop
                for _ in range(self.parameter.kmeans_ite):
                    # Assignment step
                    for i_ver in range(comp_size):
                        vertex = self.components[ind_com][i_ver]
                        vertex_data = self.main_graph.nodes[vertex]
                        obs = vertex_data['observation']

                        # Compute distances to both kernels
                        dist0 = np.sum((obs - kernels[0]) ** 2)
                        dist1 = np.sum((obs - kernels[1]) ** 2)

                        potential_label[i_ver] = dist0 > dist1

                    # Update step
                    kernels.fill(0)
                    total_weight = [0.0, 0.0]

                    for i_ver in range(comp_size):
                        vertex = self.components[ind_com][i_ver]
                        vertex_data = self.main_graph.nodes[vertex]
                        weight = vertex_data['weight']

                        if weight == 0:
                            continue

                        label_idx = int(potential_label[i_ver])
                        total_weight[label_idx] += weight
                        kernels[label_idx] += vertex_data['observation'] * weight

                    # Check for empty clusters
                    if total_weight[0] == 0 or total_weight[1] == 0:
                        break

                    # Normalize
                    for i in range(2):
                        if total_weight[i] > 0:
                            kernels[i] /= total_weight[i]

                # Compute energy for this initialization
                current_energy = 0
                for i_ver in range(comp_size):
                    vertex = self.components[ind_com][i_ver]
                    vertex_data = self.main_graph.nodes[vertex]
                    obs = vertex_data['observation']
                    weight = vertex_data['weight']
                    kernel = kernels[int(potential_label[i_ver])]
                    current_energy += np.sum((obs - kernel) ** 2) * weight

                # Keep best labeling
                if current_energy < best_energy:
                    best_energy = current_energy
                    for i_ver in range(comp_size):
                        vertex_idx = self.components[ind_com][i_ver]
                        binary_label[vertex_idx] = potential_label[i_ver]

    def compute_centers(self, centers: np.ndarray, nb_comp: int, binary_label: List[bool]):
        """Compute for each component the values of h_1 and h_2"""
        for ind_com in range(nb_comp):
            if self.saturated_components[ind_com]:
                continue
            self.compute_center(centers[ind_com], ind_com, binary_label)

    def compute_center(self, center: np.ndarray, ind_com: int, binary_label: List[bool]):
        """
        Compute centers for a single component.
        Args:
            center: shape (2, dim) array for storing the two centers
            ind_com: index of the component
            binary_label: binary labels for all vertices
        """
        # Initialize weights and centers
        total_weight = [0.0, 0.0]  # Must use list, not numpy array to match C++
        center.fill(0.0)  # Reset the centers

        # For each vertex in the component
        for i_ver in range(len(self.components[ind_com])):
            vertex_idx = self.components[ind_com][i_ver]
            node = self.main_graph.nodes[vertex_idx]

            # Skip if weight is zero
            if node['weight'] == 0:
                continue

            # Check binary label
            if binary_label[vertex_idx]:
                total_weight[0] += node['weight']
                for i_dim in range(self.dim):
                    center[0][i_dim] += node['observation'][i_dim] * node['weight']
            else:
                total_weight[1] += node['weight']
                for i_dim in range(self.dim):
                    center[1][i_dim] += node['observation'][i_dim] * node['weight']

        # Check for empty clusters
        if total_weight[0] == 0 or total_weight[1] == 0:
            # The component is saturated
            self.saturate_component(ind_com)
            center[0] = self.main_graph.nodes[self.components[ind_com][0]]['value'].copy()
            center[1] = center[0].copy()
        else:
            # Normalize by weights
            for i_dim in range(self.dim):
                center[0][i_dim] /= total_weight[0]
                center[1][i_dim] /= total_weight[1]

    def saturate_component(self, ind_com: int):
        """Mark a component as uncuttable and remove it from further graph-cuts"""
        self.saturated_components[ind_com] = True

        # For each vertex in the component
        for vertex_idx in self.components[ind_com]:
            # Set the capacities of edges to source and sink to zero
            self.main_graph[self.source][vertex_idx]['capacity'] = 0.0
            self.main_graph[vertex_idx][self.sink]['capacity'] = 0.0

    def activate_edges(self, allows_saturation: bool = True) -> int:
        """
        Analyze the optimal binary partition to detect saturated components and new activated edges
        Returns:
            int: Number of vertices in saturated components
        """
        saturation = 0
        nb_comp = len(self.components)

        # First check if the components are saturated
        for ind_com in range(nb_comp):
            if self.saturated_components[ind_com]:
                saturation += len(self.components[ind_com])
                continue

            total_weight = [0.0, 0.0]  # [sink_weight, source_weight]

            for vertex_idx in self.components[ind_com]:
                node = self.main_graph.nodes[vertex_idx]
                is_sink = (node['color'] == self.main_graph.nodes[self.sink]['color'])
                if is_sink:
                    total_weight[0] += node['weight']
                else:
                    total_weight[1] += node['weight']

            if allows_saturation and (total_weight[0] == 0 or total_weight[1] == 0):
                self.saturate_component(ind_com)
                saturation += len(self.components[ind_com])

        # Check which edges have been activated
        for u, v, edge_data in self.main_graph.edges(data=True):
            if not edge_data['real_edge']:
                continue

            color_v1 = self.main_graph.nodes[u]['color']
            color_v2 = self.main_graph.nodes[v]['color']

            color_combination = color_v1 + color_v2
            if color_combination in {0, 2, 8}:  # Same color combinations
                continue

            # The edge is active!
            edge_data['is_active'] = True
            edge_data['capacity'] = 0

            # Mark vertices as border
            self.main_graph.nodes[u]['is_border'] = True
            self.main_graph.nodes[v]['is_border'] = True

        return saturation

    def set_capacities(self, centers: np.ndarray):
        """Set capacities for the graph cut"""
        nb_comp = len(self.components)

        # Set capacities for vertex-source and vertex-sink edges
        for ind_com in range(nb_comp):
            if self.saturated_components[ind_com]:
                continue

            for vertex_idx in self.components[ind_com]:
                node = self.main_graph.nodes[vertex_idx]
                if node['weight'] == 0:
                    # No observation - no cut
                    self.main_graph[self.source][vertex_idx]['capacity'] = 0
                    self.main_graph[vertex_idx][self.sink]['capacity'] = 0
                    continue

                # Compute cost of being in B or not B
                cost_B = 0
                cost_notB = 0

                for i_dim in range(self.dim):
                    cost_B += 0.5 * node['weight'] * (
                            centers[ind_com, 0, i_dim] ** 2 -
                            2 * centers[ind_com, 0, i_dim] * node['observation'][i_dim]
                    )
                    cost_notB += 0.5 * node['weight'] * (
                            centers[ind_com, 1, i_dim] ** 2 -
                            2 * centers[ind_com, 1, i_dim] * node['observation'][i_dim]
                    )

                # Set capacities based on costs
                if cost_B > cost_notB:
                    self.main_graph[self.source][vertex_idx]['capacity'] = cost_B - cost_notB
                    self.main_graph[vertex_idx][self.sink]['capacity'] = 0.0
                else:
                    self.main_graph[self.source][vertex_idx]['capacity'] = 0.0
                    self.main_graph[vertex_idx][self.sink]['capacity'] = cost_notB - cost_B

        # Set capacities for vertex-vertex edges
        for u, v, edge_data in self.main_graph.edges(data=True):
            if not edge_data['real_edge']:
                continue

            if not edge_data['is_active']:
                edge_data['capacity'] = edge_data['weight'] * self.parameter.reg_strenth
            else:
                edge_data['capacity'] = 0

    def export_graph_state(self, prefix: str):
        """
        Export graph state focusing on maxflow-relevant attributes
        Args:
            prefix: Prefix for output files
        """
        # Export vertices focusing on flow-related attributes
        with open(f"{prefix}_vertices.txt", 'w') as f:
            f.write('vertex_id,color\n')
            for vertex_idx in self.main_graph.nodes():
                node = self.main_graph.nodes[vertex_idx]
                f.write(f"{vertex_idx},{node['color']}\n")

        # Export edges focusing on flow-related attributes
        with open(f"{prefix}_edges.txt", 'w') as f:
            f.write('source,target,capacity,residual_capacity,real_edge\n')
            for u, v, data in self.main_graph.edges(data=True):
                row = [
                    u,
                    v,
                    data['capacity'],
                    data['residual_capacity'],
                    data['real_edge']
                ]
                f.write(','.join(map(str, row)) + '\n')


    def compute_max_flow(self):
        """Compute maximum flow using NetworkX's Boykov-Kolmogorov algorithm"""
        # Create flow graph (copy capacities from main graph)
        G = self.main_graph.copy()

        try:
            # Compute the maximum flow
            R = boykov_kolmogorov(G, self.source, self.sink)

            # print("\nNon-zero flows in residual graph:")
            # flow_count = 0
            # for u, v in R.edges():
            #     flow = R[u][v].get('flow', 0)
            #     if flow > 0:
            #         print(f"Edge {u}->{v}: flow = {flow}")
            #         flow_count += 1
            # print(f"Total edges with non-zero flow: {flow_count}")

            # Update residual capacities and colors based on the result
            # reachable = set(nx.descendants(R, self.source))
            # reachable.add(self.source)
            #
            # # Update colors
            # source_color = 0
            # sink_color = 4
            #
            # # Set colors for all nodes
            # nx.set_node_attributes(self.main_graph, sink_color, 'color')
            # for node in reachable:
            #     self.main_graph.nodes[node]['color'] = source_color
            # Get the correct partition using the trees information from the residual graph
            source_tree, target_tree = R.graph["trees"]
            source_partition = set(source_tree)

            # Update vertex colors
            source_color = 0
            sink_color = 4

            # Set all nodes to sink color first
            nx.set_node_attributes(self.main_graph, sink_color, 'color')
            # Then update nodes in source partition
            for node in source_partition:
                self.main_graph.nodes[node]['color'] = source_color

            # Update residual capacities
            for u, v in R.edges():
                if self.main_graph.has_edge(u, v):
                    self.main_graph[u][v]['residual_capacity'] = R[u][v].get('flow', 0)

        except nx.NetworkXUnbounded:
            print("Error: Graph has path of infinite capacity")
            raise
        except Exception as e:
            print(f"Error in max flow computation: {str(e)}")
            raise


    def reduce(self):
        """
        Compute the reduced graph, and if need be perform a backward check.
        """
        # Compute the connected components
        self.compute_connected_components()

        if self.parameter.backward_step:
            # Compute the structure of the reduced graph
            self.compute_reduced_graph()
            # Check for beneficial merges
            self.merge(False)
        else:
            # Compute only the value associated to each connected component
            self.compute_reduced_value()

    def compute_reduced_value(self):
        """
        Compute the reduced value of each component without building the full reduced graph structure.
        This is used when backward_step is False.
        """
        # For each component, compute its reduced value
        for ind_com in range(len(self.components)):
            # Compute value and update vertices in the component
            self.compute_value(ind_com)
    def compute_connected_components(self):
        """
        Compute the connected components of the graph with active edges removed.
        """
        # Boolean vectors indicating whether edges and vertices have been seen
        edges_seen = {(u, v): False for u, v in self.main_graph.edges()}
        vertices_seen = {v: False for v in self.main_graph.nodes()}

        # Mark source and sink as seen
        vertices_seen[self.source] = True
        vertices_seen[self.sink] = True

        # Start with known roots
        old_components_size = len(self.components)
        for ind_com in range(len(self.root_vertex)):
            root = self.root_vertex[ind_com]
            if self.saturated_components[ind_com]:
                # This component is saturated, we don't need to recompute it
                for vertex_idx in self.components[ind_com]:
                    vertices_seen[vertex_idx] = True
            else:
                # Compute the new content of this component
                self.components[ind_com] = self.connected_comp_from_root(
                    root, len(self.components[ind_com]), vertices_seen, edges_seen)

        # Look for components that did not already exist
        for vertex_idx in range(self.n_vertex):
            if vertices_seen[vertex_idx]:
                continue

            # This vertex becomes the root of a new component
            root = vertex_idx
            current_component_size = len(self.components[self.main_graph.nodes[root]['in_component']])

            new_component = self.connected_comp_from_root(
                root, current_component_size, vertices_seen, edges_seen)
            self.components.append(new_component)
            self.root_vertex.append(root)
            self.saturated_components.append(False)


    def connected_comp_from_root(self, root: int, size_comp: int,
                                 vertices_seen: Dict[int, bool],
                                 edges_seen: Dict[Tuple[int, int], bool]) -> List[int]:
        """
        Compute the connected component associated with root by performing depth-first search.
        """
        vertices_added = []  # vertices in current connected component
        vertices_to_add = [root]  # vertices to be added (stack)

        while vertices_to_add:
            vertex_current = vertices_to_add.pop()

            if vertices_seen[vertex_current]:
                continue

            vertices_added.append(vertex_current)
            vertices_seen[vertex_current] = True

            # Explore neighbors
            for _, neighbor, edge_data in self.main_graph.edges(vertex_current, data=True):
                if edge_data['is_active'] or edges_seen.get((vertex_current, neighbor), False):
                    continue

                edge_seen = edges_seen.get((vertex_current, neighbor), False)
                edge_reverse_seen = edges_seen.get((neighbor, vertex_current), False)

                if not (edge_seen or edge_reverse_seen):
                    edges_seen[(vertex_current, neighbor)] = True
                    edges_seen[(neighbor, vertex_current)] = True
                    vertices_to_add.append(neighbor)

        return vertices_added


    def compute_reduced_graph(self):
        """
        Compute the adjacency structure between components and
        the weight and value of each component.
        """
        # Create new reduced graph
        self.reduced_graph = nx.DiGraph()

        # Fill the values of the reduced graph
        for ind_com in range(len(self.components)):
            component_values_and_weight = self.compute_value(ind_com)

            # Add vertex to reduced graph with computed values
            self.reduced_graph.add_node(ind_com,
                                        weight=component_values_and_weight[1],
                                        value=component_values_and_weight[0])

        # Compute edges of reduced graph
        self.borders = []
        ind_border_edge = 0

        # Look through all edges in main graph
        for u, v, edge_data in self.main_graph.edges(data=True):
            if not edge_data['real_edge']:
                continue

            # Get components of source and target
            comp1 = self.main_graph.nodes[u]['in_component']
            comp2 = self.main_graph.nodes[v]['in_component']

            if comp1 == comp2:
                continue

            # By convention, component_source is smallest index
            component_source = min(comp1, comp2)
            component_target = max(comp1, comp2)

            # Add border-edge if it doesn't exist
            if not self.reduced_graph.has_edge(component_source, component_target):
                self.reduced_graph.add_edge(component_source, component_target,
                                            index=ind_border_edge,
                                            weight=0)
                ind_border_edge += 1
                self.borders.append([])

            # Add weight of current edge to border edge
            self.reduced_graph[component_source][component_target]['weight'] += 0.5 * edge_data['weight']
            border_index = self.reduced_graph[component_source][component_target]['index']
            self.borders[border_index].append((u, v))


    def cutoff(self):
        """Apply the cutoff operation."""
        i = 0
        while True:
            self.compute_reduced_graph()
            n_merged = self.merge(True)  # True indicates this is a cutoff merge
            i += 1
            if n_merged == 0 or i > 50:
                break


    def merge(self, is_cutoff: bool) -> int:
        """
        Check whether the energy can be decreased by removing edges from the reduced graph.
        Returns number of components merged.
        """
        # Priority queue for potential merges
        merge_queue = []

        # Go through all edges in reduced graph to compute potential gains
        for source_comp, target_comp, border_data in self.reduced_graph.edges(data=True):
            if is_cutoff:
                # Check cutoff condition
                if (self.reduced_graph.nodes[source_comp]['weight'] >= self.parameter.cutoff and
                        self.reduced_graph.nodes[target_comp]['weight'] >= self.parameter.cutoff):
                    continue

            # Compute gain from merging
            merge_gain = self.compute_merge_gain(source_comp, target_comp)

            if is_cutoff or merge_gain[1] > 0:
                # Create merge information
                heapq.heappush(merge_queue,
                               (-merge_gain[1],  # negative for max-heap
                                source_comp,
                                target_comp,
                                border_data['index'],
                                merge_gain[0]))  # merged_value

        # Process merges
        n_merged = 0
        is_merged = [False] * len(self.components)
        to_destroy = [False] * len(self.components)

        while merge_queue:
            gain, comp1, comp2, border_index, merged_value = heapq.heappop(merge_queue)
            gain = -gain  # convert back to positive

            if not is_cutoff and gain <= 0:
                break

            if is_merged[comp1] or is_merged[comp2]:
                continue

            n_merged += 1

            # Merge components
            self.components[comp1].extend(self.components[comp2])
            self.saturated_components[comp1] = False

            # Update reduced graph
            self.reduced_graph.nodes[comp1]['weight'] += self.reduced_graph.nodes[comp2]['weight']
            self.reduced_graph.nodes[comp1]['value'] = merged_value

            # Deactivate border between comp1 and comp2
            for edge in self.borders[border_index]:
                self.main_graph.edges[edge]['is_active'] = False

            is_merged[comp1] = True
            is_merged[comp2] = True
            to_destroy[comp2] = True

        # Rebuild components, root_vertex, and saturated_components
        new_components = []
        new_root_vertex = []
        new_saturated_components = []

        for ind_com in range(len(self.components)):
            if to_destroy[ind_com]:
                continue

            new_components.append(self.components[ind_com])
            new_root_vertex.append(self.root_vertex[ind_com])
            new_saturated_components.append(self.saturated_components[ind_com])

            # Update vertex information
            for vertex_idx in self.components[ind_com]:
                self.main_graph.nodes[vertex_idx]['value'] = self.reduced_graph.nodes[ind_com]['value']
                self.main_graph.nodes[vertex_idx]['in_component'] = len(new_components) - 1

        self.components = new_components
        self.root_vertex = new_root_vertex
        self.saturated_components = new_saturated_components

        return n_merged


    def compute_merge_gain(self, comp1: int, comp2: int) -> Tuple[np.ndarray, float]:
        """
        Compute the gain obtained by merging two components.
        Returns:
            Tuple of (merged_value, gain)
        """
        merge_value = np.zeros(self.dim)
        gain = 0.0

        # Get component attributes
        comp1_data = self.reduced_graph.nodes[comp1]
        comp2_data = self.reduced_graph.nodes[comp2]

        # Compute merged value
        total_weight = comp1_data['weight'] + comp2_data['weight']
        for i_dim in range(self.dim):
            merge_value[i_dim] = (comp1_data['weight'] * comp1_data['value'][i_dim] +
                                  comp2_data['weight'] * comp2_data['value'][i_dim]) / total_weight

            # Compute energy gain
            gain += 0.5 * (pow(merge_value[i_dim], 2) * total_weight -
                           pow(comp1_data['value'][i_dim], 2) * comp1_data['weight'] -
                           pow(comp2_data['value'][i_dim], 2) * comp2_data['weight'])

        return merge_value, gain

def setup_cp(n_nodes: int, n_edges: int, n_obs: int,
             observation: np.ndarray, eu: np.ndarray, ev: np.ndarray,
             edge_weight: np.ndarray, node_weight: np.ndarray) -> CutPursuit:
    """Set up the Cut Pursuit problem"""
    cp = CutPursuit(n_nodes)
    cp.dim = n_obs
    cp.n_vertex = n_nodes

    # Add all vertices including source and sink from the start
    cp.source = n_nodes  # source will be the last vertex
    cp.sink = n_nodes + 1  # sink will be the last vertex + 1

    # Add regular vertices with their attributes
    for ind_nod in range(n_nodes):
        cp.main_graph.add_node(ind_nod,
                               weight=node_weight[ind_nod],
                               observation=observation[ind_nod].copy(),
                               value=np.zeros(n_obs),
                               color=-1,
                               is_border=False,
                               in_component=0)

    # Add source and sink at initialization
    cp.main_graph.add_node(cp.source,
                           weight=1.0,
                           observation=np.zeros(n_obs),
                           value=np.zeros(n_obs),
                           color=-1,
                           is_border=False,
                           in_component=0)
    cp.main_graph.add_node(cp.sink,
                           weight=1.0,
                           observation=np.zeros(n_obs),
                           value=np.zeros(n_obs),
                           color=-1,
                           is_border=False,
                           in_component=0)

    # Add edges with their attributes
    true_ind_edg = 0
    for ind_edg in range(n_edges):
        source, target = eu[ind_edg], ev[ind_edg]
        weight = edge_weight[ind_edg]

        # Add forward edge
        cp.main_graph.add_edge(source, target,
                               index=true_ind_edg,
                               weight=weight,
                               capacity=weight,
                               residual_capacity=0.0,
                               is_active=False,
                               real_edge=True)

        # Add reverse edge
        cp.main_graph.add_edge(target, source,
                               index=true_ind_edg + 1,
                               weight=weight,
                               capacity=weight,
                               residual_capacity=0.0,
                               is_active=False,
                               real_edge=True)

        true_ind_edg += 2

    cp.n_edge = true_ind_edg
    return cp

# def cut_pursuit(n_nodes, n_edges, n_obs, y, eu, ev, edge_weight, node_weight, lambda_, cutoff, mode, speed, weight_decay, verbose):

def cut_pursuit(n_nodes: int, n_edges: int, n_obs: int,
                observation: np.ndarray, eu: np.ndarray, ev: np.ndarray,
                edge_weight: np.ndarray, node_weight: np.ndarray,
                lambda_: float, cutoff: int, mode: float,
                speed: float, weight_decay: float,
                verbose: float = 0) -> Tuple[np.ndarray, List[List[int]], np.ndarray, np.ndarray, np.ndarray]:
    """
    Main function to run Cut Pursuit algorithm.

    Args:
        n_nodes: Number of nodes
        n_edges: Number of edges
        n_obs: Dimension of observation
        observation: Node observations array (n_nodes x n_obs)
        eu: Edge source nodes
        ev: Edge target nodes
        edge_weight: Edge weights
        node_weight: Node weights
        lambda_: Regularization strength
        cutoff: Minimal component size
        mode: Not used in this version
        speed: Not used in this version
        weight_decay: Weight decay for flow steps
        verbose: Verbosity level

    Returns:
        Tuple of:
            - solution: ndarray of shape (n_nodes, n_obs) containing computed values
            - components: List[List[int]] containing component assignments
            - in_component: ndarray of shape (n_nodes,) containing component indices
            - energy_out: ndarray containing energy values at each iteration
            - time_out: ndarray containing computation times at each iteration
    """
    # Set random seed for reproducibility
    np.random.seed(1)

    if verbose > 0:
        print("L0-CUT PURSUIT")

    # Initialize and parameterize the Cut Pursuit solver
    cp = setup_cp(n_nodes, n_edges, n_obs, observation, eu, ev, edge_weight, node_weight)

    cp.parameter.flow_steps = 4
    cp.parameter.weight_decay = weight_decay
    cp.parameter.kmeans_ite = 8
    cp.parameter.kmeans_resampling = 5
    # cp.parameter.max_ite_main = 2
    cp.parameter.max_ite_main = 20
    cp.parameter.backward_step = True
    cp.parameter.stopping_ratio = 0.001
    cp.parameter.reg_strenth = lambda_
    cp.parameter.cutoff = cutoff

    # Run the optimization
    energy_out, time_out = cp.run()

    # Get the reduced graph structure
    cp.compute_reduced_graph()

    # Get number of components
    n_nodes_red = len(cp.components)

    # Initialize return arrays
    solution = np.zeros((n_nodes, n_obs))
    in_component = np.zeros(n_nodes, dtype=np.int32)
    components = [[] for _ in range(n_nodes_red)]

    # Fill solution array with computed values
    for vertex_idx in range(n_nodes):
        solution[vertex_idx] = cp.main_graph.nodes[vertex_idx]['value']

    # Fill components and in_component arrays
    for ind_com in range(n_nodes_red):
        components[ind_com].extend(cp.components[ind_com])

        for vertex_idx in cp.components[ind_com]:
            in_component[vertex_idx] = cp.main_graph.nodes[vertex_idx]['in_component']

    return solution, components, in_component, np.array(energy_out), np.array(time_out)

# def perform_cut_pursuit(K, regStrength, pc, edgeWeight,Eu,Ev,in_component,components):
def perform_cut_pursuit(K, lambda_, pc):
    point_count=len(pc)
    if point_count == 0:
        return False

    kdtree = cKDTree(pc[:,:3])
    nn_D, nn_idx = kdtree.query(pc, k=K + 1)

    # Remove self-connections
    # distances = nn_D[:, 1:]
    indices = nn_idx[:, 1:]

    # Create edge list
    n_nodes = len(pc)
    n_obs=3
    n_edges=n_nodes*K
    cutoff = 0
    mode = 1.0
    speed = 0
    weight_decay = 0
    verbose = 1

    eu = np.repeat(np.arange(n_nodes), K)
    ev = indices.ravel()

    y=pc[:,:3]-np.mean(pc[:,:3],axis=0)

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
        lambda_=lambda_,
        cutoff=cutoff,
        mode=mode,
        speed=speed,
        weight_decay=weight_decay,
        verbose=verbose
    )
    return in_component







def decimate_pcd(columns,min_res):
    xyz_min = np.min(columns[:, :3], axis=0)
    xyz_max = np.max(columns[:, :3], axis=0)

    block_shape = np.floor((xyz_max[:3] - xyz_min[:3]) / min_res).astype(np.int32) + 1
    block_shape = block_shape[[1, 0, 2]]
    block_x = xyz_max[1] - columns[:, 1]
    block_y = columns[:, 0] - xyz_min[0]
    block_z = columns[:, 2] - xyz_min[2]
    block_ijk = np.floor(np.concatenate([block_x[:, np.newaxis], block_y[:, np.newaxis], block_z[:, np.newaxis]],axis=1) / min_res).astype(np.int32)
    block_idx = np.ravel_multi_index((np.transpose(block_ijk)).astype(np.int32), block_shape)
    _, block_idx_uidx, block_inverse_idx = np.unique(block_idx, return_index=True, return_inverse=True)
    # columns_dec = columns[block_idx_uidx]
    return block_idx_uidx,block_inverse_idx

if __name__ == "__main__":
    # Example usage with mock point cloud data
    # pc = [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0], [2.0, 3.0, 4.0], [3.0, 4.0, 5.0]]  # Replace with actual point cloud data
    K = 4
    regStrength = 1.0


    import os
    import laspy
    # path_to_las = r"..\data\JP10_plot_2cm_test2.las"
    path_to_las = r"..\data\JP10_plot_2cm_test3.las"
    # path_to_las = r"F:\prj\CC2\comp\TreeAIBox\test\testCutPursuit\JP10_plot_2cm_test2_treeiso1.las"
    pcd_basename = os.path.basename(path_to_las)[:-4]
    las = laspy.read(path_to_las)
    pcd=np.transpose([las.x,las.y,las.z])

    min_res=0.05
    reg_strength = 1.0

    dec_idx_uidx, dec_inverse_idx = decimate_pcd(pcd[:, :3], min_res)  # reduce points first
    pcd_dec = pcd[dec_idx_uidx]

    in_component=perform_cut_pursuit(K, regStrength, pcd_dec)

    # labels = perform_fast_cut_pursuit(pcd_dec, k=10, reg_strength=1.0, cutoff=0)
    # reg_strength = 0.5  # Adjust this parameter to control segmentation granularity
    #
    # labels=shortestpath3D(pcd, pcd[:,-1],pcd[:,-2], min_res,max_isolated_distance)
    #
    # # pred=np.zeros(len(pcd))
    # # pred[stem_ind]=labels
    # # np.savetxt(r"C:\Users\ZXI\Downloads\tmp.txt",np.concatenate([pcd,labels[:,np.newaxis]],axis=-1))
    las.add_extra_dim(laspy.ExtraBytesParams(name="cutpursuit", type="int32", description="cutpursuit"))
    las.cutpursuit = in_component[dec_inverse_idx]
    las.write(r"C:\Users\ZXI\Downloads\JP10_plot_2cm_test3_cutpursuit.laz")
    # las.write(r"F:\prj\CC2\comp\TreeAIBox\test\testCutPursuit\JP10_plot_2cm_test2_cutpursuit_python.laz")



