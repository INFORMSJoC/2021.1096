[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The data in this repository is an archive of the data that were used in the research reported on in the paper 
[Distributionally Robust Optimization under Decision-Dependent Ambiguity Set with Applications to Machine Scheduling and Humanitarian Logistics](https://doi.org/10.1287/ijoc.2021.1096) 
by N. Noyan and G. Rudolf and M.A. Lejeune.

M. Lejeune acknowledges the support of the Office of Naval Research [Grant N000141712420] for this study.

## Cite

To cite this software, please cite the [paper](https://doi.org/10.1287/ijoc.2020.1022) and the software

Below is the BibTex for citing this version of the data.

```
@article{Data.IJOC.2021.1096,
  author =        {N. Noyan and G. Rudolf and M.A. Lejeune},
  publisher =     {INFORMS Journal on Computing},
  title =         {Data for Distributionally Robust Optimization under Decision-Dependent Ambiguity Set with Applications to Machine Scheduling and Humanitarian Logistics},
  year =          {2020},
  doi =           {10.5281/zenodo.xxxxxxx},
  url =           {https://github.com/INFORMSJoC/2021.1096},
}  
```

## Content

This repository includes
1. **AMPL Data Files** are in the folder [FinalDataFiles](AMPLFiles/FinalDataFiles). We consider ten instances, and index the data folders from 1 to 10, accordingly. 
1. An **AMPL Mod File** [here](AMPLFiles/FinalModFile.mod)
1. **AMPL Run Files** created to solve a particular set of MIP formulations, allowing efficient batch runs. For example, [FinalAlternativeFormulations_COM.run]  (AMPLFiles/FinalAlternativeFormulations_COM.run) was used to solve the following four MIP formulations under a particular set of parameter setting: 
    * Cost minimizing basic completely comonotone MILP (CCM)
    * Budget-constrained basic completely comonotone MILP (CCM)
    * Cost minimizing completely comonotone MILP with RLT (CCM-RLT)
    * Budget-Constrained completely comonotone MILP with (CCM-RLT). 
    The comments on the mod file along with the problem definitions in the main run files provide insights about the alternative MIP formulations considered in our study.
1. We created additional run files to get results for a batch of main run files under different parameter settings. These files are provided in the folder [Batch-RunFiles](AMPLFiles/Batch-RunFiles) For illustrative purposes. Each type of file is given for a particular instance index such as `DataSet5`; we basically modify the index information to get results for the other instances.
1. All the output files obtained in our computational study are available in the folder [OutputFiles](AMPLFiles/OutputFiles). Since there is a very large number of output files (almost two thousand problem instances), for convenience, we also provided the excel files including the key results retrieved from the output files. 

We next outline how to associate the output files with the tables and figures presented in the paper.

1. Figure 1: Optimal objective function value (robustified CVaR of TWCT) for varying radius and budget.
The corresponding output files are available under the following folders: 
    * [DataSet1](AMPLFiles/OutputFiles/ModelAnalysis/DataSet1), where the combined key results are summarized in the excel file [summary_outputs_modelanalysis.xlsx](AmplFiles/OutputFiles/summary_outputs_modelanalysis.xlsx) 
    * [Fixing](AMPLFiles/OutputFiles/ModelAnalysis/DataSet1/Fixing) where “fixing” refers to setting all the control decisions to 1, i.e., the setting with “no compression decisions.” The combined key results are summarized in the excel file [summary_outputs_fixing.xlsx](AmplFiles/OutputFiles/summary_outputs_fixing.xlsx)
1. Table 1: Computational performance – comonotone data, Table 2: Impact of modeling parameters on performance of CCM-RLT, and Table 5: Computational performance on non-	  comonotone instances.
    * The corresponding output files are available in the folder [ComputationalPerformance](AMPLFiles/OutputFiles/ComputationalPerformance), where the combined key results are summarized in the excel file (summary_outputs_computational_performance.xslx)[AMPLFiles/OutputFiles/summary_outputs_computational_performance.xslx].
1. Figure 3: Robustified CVaR of TWCT versus Total Compression Cost for increasing   (trade-off coefficient) values in {0.1, 0.2,…0.9,1}.
    * The corresponding output files are available in the folder (ParetoAnalysis)[AMPLFiles/OutputFiles/ModelAnalysis/DataSet1/ParetoAnalysis], where the combined key results are summarized in the excel file [summary_outputs_pareto.xslx](AMPLFiles/OutputFiles/summary_outputs_pareto.xslx).
1. Figure 4: Optimal objective function values and solutions for illustrative example.
    * The corresponding output files are available in the folder [ToyExample](AMPLFiles/OutputFiles/ToyExample)
1. Table 4: Impact of radius ( ) on performance of CCM-RLT and Figure 5: Solution times for varying radius.
    * The corresponding output files are available in the folder [KappaAnalysis](AMPLFiles/OutputFiles/ComputationalPerformance/KappaAnalysis), where the combined key results are summarized in the excel file [summary_outputs_computational_performance_kappaimpact](AMPLFiles/OutputFiles/summary_outputs_computational_performance_kappaimpact.xslx)
