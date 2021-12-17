# GraSP: Available third party toolboxes

This page describes the available toolboxes from GraSP.
These are intalled in the subfolder [3rdParty/](3rdParty/) when [`grasp_install()`](grasp_install.m) is executed.

For the optional toolboxes, these can be started using [`grasp_start_opt_3rd_party('name')`](grasp_start_opt_3rd_party.m) where `name` is the name of the toolbox.

Details can be found either in the link to the toolbox (header), or to the associated scientific communication.

## Graph Signal Processing

### [GSPbox](https://github.com/epfl-lts2/gspbox/)

Name | Optional | Paper(s)
---|:---:|---
`gspbox` | Yes | [https://arxiv.org/abs/1408.5781]()

Unmaintained Matlab toolbox for graph signal processing.
Some third party toolboxes, or codes from the community, may need access to its functions.
Note that Matlab structures can be converted between GSPbox and GraSP using [Util/grasp_to_gspbox](Util/grasp_to_gspbox.m) and [Util/grasp_from_gspbox](Util/grasp_from_gspbox.m).

### [Spectral Graph Wavelet Transform](https://wiki.epfl.ch/sgwt/)

Name | Optional | Paper(s)
---|:---:|---
`sgwt` | Yes | [http://dx.doi.org/10.1016/j.acha.2010.04.005]()

The classical graph wavelet transform from David Hammond, Pierre Vandergheynst, Rémi Gribonval.

### [Active Semi-Supervised Learning with Sampling Theory](https://github.com/STAC-USC/Active_SSL_with_Sampling_Theory/)

Name | Optional | Paper(s)
---|:---:|---
`usc_ssl_sampling` | Yes | [https://doi.org/10.1145/2623330.2623760]()

Toolbox to perform active sampling (sampled vertex set selection) and signal reconstruction from a subset of samples (semi-supervised learning).

### [Graph Laplacian Learning](https://github.com/STAC-USC/Graph_Learning/)

Name | Optional | Paper(s)
---|:---:|---
`usc_graph_learning` | Yes | [https://doi.org/10.1109/JSTSP.2017.2726975]() <br /> [https://doi.org/10.1109/TSIPN.2018.2872157]()

Graph Learning using Gaussian Markov Random Fields models.

### [Discriminative Graph Learning with Sparsity Regularization](https://github.com/STAC-USC/Disc-GLasso/)

Name | Optional | Paper(s)
---|:---:|---
`usc_disc_glasso` | Yes | [https://doi.org/10.1109/ICASSP.2017.7952698]()

Perform graph learning using Graphical LASSO.

### [Fast Implementation for Graph Fourier Transforms based on Grid Symmetries](https://github.com/STAC-USC/symmetric_grid/)

Name | Optional | Paper(s)
---|:---:|---
`usc_symmetric_grid` | Yes | [https://doi.org/10.1109/ICASSP.2017.7952929]()

Exploits symmetries in the graph to implement a fast graph Fourier transform.

### [FastGSSS (Fast Graph Sampling Set Selection) Toolbox](https://github.com/ychtanaka/FastGSSS/)

Name | Optional | Paper(s)
---|:---:|---
`fast_gsss` | Yes | [https://doi.org/10.1109/TSP.2019.2908129]()

Matlab toolbox for eigendecomposition-free sampling set selection.

### [Graph Learning with Monotone Topology Properties](https://github.com/STAC-USC/graph_learning_properties/)

Name | Optional | Paper(s)
---|:---:|---
`usc_graph_learning_properties` | Yes | [https://doi.org/10.1109/JSTSP.2017.2726975]() <br /> [https://doi.org/10.1109/TSP.2018.2813337]()

Implementation of a graph Learning approach based on the `usc_graph_learning` toolbox above.

## Graphs

### [GrTheory](https://www.mathworks.com/matlabcentral/fileexchange/4266-grtheory-graph-theory-toolbox)

Name | Optional | Paper(s)
---|:---:|---
`GrTheory` | Yes |

This toolbox implements a number of useful functions to handle and transform graphs.

### [MatlabBGL](https://www.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl)

Name | Optional | Paper(s)
---|:---:|---
`MatlabBGL` | Yes | 

This is another toolbox implementing useful function for graphs.
GraSP is however only using the Minnesota graph data from this toolbox (the dataset provided by Matlab does not include node labels).


### [USC STAC Graph Structures](https://github.com/STAC-USC/GraphStructures/)

Name | Optional | Paper(s)
---|:---:|---
`usc_graphs` | Yes | 

Several graph structures for GraSP used the USC's lab STAC.

### [Molène Dataset](https://github.com/bgirault-usc/Molene-Dataset/)

Name | Optional | Paper(s)
---|:---:|---
`bgirault_molene_dataset` | Yes | [https://doi.org/10.1109/EUSIPCO.2015.7362637]()

Molène dataset, preprocessed for Matlab from [data.gouv.fr](https://www.data.gouv.fr/fr/datasets/donnees-horaires-des-55-stations-terrestres-de-la-zone-large-molene-sur-un-mois/) by Benjamin Girault.

## Miscellaneous tools

### [MyPatcher](https://github.com/bgirault-usc/MyPatcher)

Name | Optional | Paper(s)
---|:---:|---
| No |

This is a simple Matlab function to apply a patch (obtained through `diff -u`).
In the context of GraSP, this is most useful to apply a patch to a third party toolbox after installing it.
    
### [Anneal](https://www.mathworks.com/matlabcentral/fileexchange/10548-general-simulated-annealing-algorithm)

Name | Optional | Paper(s)
---|:---:|---
`anneal` | Yes | 

This is a simple library to perform optimization using simulated annealing.
Used in grasp to implement the approach described in [Duality between Time Series and Networks](https://doi.org/10.1371/journal.pone.0023378).

### [Chebfun](https://github.com/chebfun/chebfun/)

Name | Optional | Paper(s)
---|:---:|---
`chebfun` | Yes | [PDF Guide](https://www.chebfun.org/docs/guide/chebfun_guide.pdf)

Toolbox for Matlab to perform Chebyshev polynomial approximations.
