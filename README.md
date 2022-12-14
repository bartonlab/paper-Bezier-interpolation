# Overview 

This repository contains data and scripts for reproducing the results accompanying the manuscript

### Bézier interpolation improves the inference of dynamical models from data ###

Kai Shimagaki<sup>1</sup> and John P. Barton<sup>1,2,#</sup>

<sup>1</sup> Department of Physics and Astronomy, University of California, Riverside<br>
<sup>2</sup> Department of Computational and Systems Biology, University of Pittsburgh School of Medicine

<sup>#</sup> correspondence to [jpbarton@pitt.edu](mailto:jpbarton@pitt.edu)

# Contents 

Scripts for the analysis and creating figures are organized in `note` a directory, in a Jupyter Notebook as shown below:


1. Scripts for generating sequences from Wright-Fisher models can be found in the `Wright-Fisher_model_binary.ipynb` notebook. Since this process takes around 1 hours, we also provide pre-generated Wright-Fisher data sets that can be found in the [Zenodo record](https://doi.org/10.5281/zenodo.7153464). After downloading, need to extract the contents of the archive `data.zip` directly under `paper-Bezier-interpolation`.
2. Scripts for analyzing and inferring genetic selection coefficients for Wright-Fisher models are contained in the `Inference_Wright-Fisher.ipynb` notebook. 
3. Scripts for analysis and inferring interaction parameters of Ornstein-Uhlenbeck process are organized in the `Inference_Ornstein-Uhlenbeck.ipynb` notebook. 
4. Scripts for analyzing and inferring genetic selection coefficients for HIV-1 data are located in the `Inference_HIV.ipynb notebook`. The required input data set can be found in the previous study, [Sohail, et.al](https://www.nature.com/articles/s41587-020-0737-3). 
5. Finally, scripts for analysis and figures contained in the manuscript are located in the `figures.ipynb` notebook. 

Our preprint is available at https://arxiv.org/abs/2209.11285 .

# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).
