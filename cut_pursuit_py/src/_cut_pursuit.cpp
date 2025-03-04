#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include <pybind11/stl.h>

#include "types.hpp"
#include "cp_d0_dist.hpp"

namespace py = pybind11;

// Wrapper function to perform cut pursuit segmentation
py::array_t<index_t> perform_cut_pursuit(
    const real_t reg_strength,
    size_t D,
    py::array_t<real_t, py::array::c_style> pc_vec,
    py::array_t<real_t, py::array::c_style> edge_weights,
    py::array_t<index_t, py::array::c_style> Eu,
    py::array_t<index_t, py::array::c_style> Ev,
	bool verbose
) {
    // Validate inputs
    auto pc_buf = pc_vec.request();
    auto ew_buf = edge_weights.request();
    auto eu_buf = Eu.request();
    auto ev_buf = Ev.request();

    if (pc_buf.ndim != 2 || pc_buf.shape[1] != D) {
        throw std::runtime_error("Point cloud must be 2D array with shape (N, D)");
    }

    const index_t pointCount = static_cast<index_t>(pc_buf.shape[0]);
    const index_t E = static_cast<index_t>(eu_buf.shape[0]);

    if (ew_buf.size != E || ev_buf.size != E) {
        throw std::runtime_error("Edge weights and edge vertex arrays must match");
    }

    // Convert point cloud to flat array for observations
    std::vector<real_t> Y(pointCount * D);
    real_t* pc_ptr = static_cast<real_t*>(pc_buf.ptr);
    for (size_t i = 0; i < pointCount; i++) {
        for (size_t d = 0; d < D; d++) {
            Y[i * D + d] = pc_ptr[i * D + d];
        }
    }

    // Compute first_edge from Eu
    std::vector<index_t> first_edge(pointCount + 1, 0);
    std::vector<index_t> adj_vertices(E);
    std::vector<real_t> edge_weights_vec(E);

    // Get pointers to buffer data
    index_t* eu_ptr = static_cast<index_t*>(eu_buf.ptr);
    index_t* ev_ptr = static_cast<index_t*>(ev_buf.ptr);
    real_t* ew_ptr = static_cast<real_t*>(ew_buf.ptr);

    // Count edges per vertex
    for (index_t i = 0; i < E; i++) {
        first_edge[eu_ptr[i] + 1]++;
        adj_vertices[i] = ev_ptr[i];
        edge_weights_vec[i] = ew_ptr[i];
    }

    // Compute cumulative sum for first_edge
    for (size_t i = 1; i <= pointCount; i++) {
        first_edge[i] += first_edge[i-1];
    }

    // Create cut pursuit instance
    using CP = Cp_d0_dist<real_t, index_t, comp_t>;
    
    CP* cp = new CP(pointCount, E, first_edge.data(), adj_vertices.data(), Y.data(), D);

    // Set edge weights with regularization parameter
    cp->set_edge_weights(edge_weights_vec.data(), reg_strength);

    // Set loss function to quadratic
    cp->set_loss(cp->quadratic_loss());

    // Set solver parameters
    cp->set_cp_param(1e-4, 20, verbose);  // tolerance, max iterations, verbose

    // Run the cut pursuit algorithm
    cp->cut_pursuit(true);

    // Get components assignment
    const comp_t* comp_assign;
    cp->get_components(&comp_assign);

    // Copy results
    std::vector<index_t> in_component(pointCount);
    for (size_t i = 0; i < pointCount; i++) {
        in_component[i] = static_cast<index_t>(comp_assign[i]);
    }

    delete cp;
	
	// Create a NumPy array of the right size
    py::array_t<index_t> result(pointCount);
    auto buf = result.request();
    index_t* ptr = static_cast<index_t*>(buf.ptr);
    
    // Copy data from the vector into the NumPy array
    std::copy(in_component.begin(), in_component.end(), ptr);
	
    return result;
}

PYBIND11_MODULE(_cut_pursuit, m) {
    m.doc() = "Cut Pursuit Implementation for Point Cloud Segmentation";
    
    m.def("perform_cut_pursuit", &perform_cut_pursuit, 
        "Perform cut pursuit segmentation",
        py::arg("reg_strength"),
        py::arg("D"),
        py::arg("pc_vec"),
        py::arg("edge_weights"),
        py::arg("Eu"),
        py::arg("Ev"),
		py::arg("verbose"));
}