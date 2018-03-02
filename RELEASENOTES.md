# Release Notes

## Version 1.1

### Major Features

 * `grasp_start_opt_3rd_party`: Starts optional toolboxes.
 * `grasp_bibliography`: Returns the papers to add to the references based on the started third party toolboxes.
 * `grasp_eigendecomposition`: Stable decomposition of directed adjacency using block diagonal Schur decomposition [1].
 * `grasp_gft_gui`: Port to GraSP of a GUI displaying the GFT of a graph [2].
 * `grasp_eigendecomposition`: irregularity aware GFT [4].
 * `grasp_voronoi_areas`: Voronoi cell area inner product matrix [4].
 * `grasp_gfm_l1`: first graph Fourier modes computation using an l1-norm graph variation [4].

### New Third Party Toolboxes

 * [GSPBox](https://github.com/epfl-lts2/gspbox/)
 * [SGWT](http://wiki.epfl.ch/sgwt/)
 * Toolboxes from [USC - STAC](https://github.com/STAC-USC/).

### Minor Features

 * Scripts to generate the toolbox online documentation.
 * `grasp_directed_torus`: one argument for equal dimensions.
 * `grasp_translation`: custom definition of frequencies, or equivalently builds an isometric graph operator.
 * `grasp_minnesota`: more adjacency matrix construction schemes.
 * `grasp_translation`: new translation operator based a symmetric graph shift [3].
 * `grasp_convolution`: convolutive operator output (matrix) if only two arguments.
 * `grasp_plane_rnd`: options to control the threshold on weights of displayed edges.
 * `grasp_show_fouriermodes`: pass along options to `grasp_show_graph`.
 * `grasp_adjacency_gaussian`, `grasp_adjacency_degreenorm`: degree-normalized weights.
 * `grasp_subgraph`: perform a subgraph computation from a subset of vertices.
 * `grasp_show_graph`: reworked design of the arrow for directed edges.
 * `grasp_show_graph`: allows for automatic computation of the graph boundaries.

### Bugfixes

 * grasp_adjacency_knn: incorrect weights fixed.
 * Latex `tikzgraph` command: append Latex/TikZ code before the boundaries and the scale are drawn.
 * `grasp_subaxis_matrix_dimensions`: background was not set with 2 arguments with a numerical 2nd argument.
 * `grasp_show_graph`: Text in 3D fixed.
 * `grasp_exportcsv_signal`: there was an error when using sparse signals.

### References

[1] https://hal.archives-ouvertes.fr/tel-01256044

[2] http://biron.usc.edu/wiki/index.php/Graph_Fourier_Transform_Interactive_GUI

[3] https://doi.org/10.1109/GlobalSIP.2016.7905858

[4] https://hal.archives-ouvertes.fr/hal-01708695
