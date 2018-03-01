# Release Notes

## Version 1.1

### Major Features

 * `grasp_start_opt_3rd_party`: Starts optional toolboxes.
 * `grasp_bibliography`: Returns the papers to add to the references based on the started third party toolboxes.
 * `grasp_eigendecomposition`: Stable decomposition of directed adjacency using block diagonal Schur decomposition [1].
 * `grasp_gft_gui`: Port to GraSP of a GUI displaying the GFT of a graph [2].

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

### Bugfixes

 * grasp_adjacency_knn: incorrect weights fixed.
 * Latex `tikzgraph` command: append Latex/TikZ code before the boundaries and the scale are drawn.
 * `grasp_subaxis_matrix_dimensions`: background was not set with 2 arguments with a numerical 2nd argument.
 * `grasp_show_graph`: Text in 3D fixed.

### References

[1] https://hal.archives-ouvertes.fr/tel-01256044
[2] http://biron.usc.edu/wiki/index.php/Graph_Fourier_Transform_Interactive_GUI
[3] https://doi.org/10.1109/GlobalSIP.2016.7905858