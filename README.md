# PC-SAFT-Crystallization-Workflow
Associated Content: Nethrue Pramuditha Mendis, Jiayuan Wang, Richard Lakerveld, "A Workflow for Crystallization Process Design with Simultaneous Process Optimization and Solvent Selection based on the Perturbed-Chain Statistical Associating Fluid Theory", Chemie Ingenieur Technik, 2022

Prerequisites:
MATLAB (The MathWorks, Inc.)
GAMS (GAMS Development Corporation) with CONOPT solver license
MATLAB function file *vert2con* (can be downloaded at https://www.mathworks.com/matlabcentral/fileexchange/7895-vert2con-vertices-to-constraints)

Provided in the git repository:
The MATLAB function files *paraest*, *relxopt*, *fixedopt*, and *mischeck* with associated GAMS files
PC-SAFT parameters of the solvents used in the manuscript

Setting up:
1. Interface GAMS and MATLAB with GDXMRW procedure (details can be found at https://www.gams.com/latest/docs/T_GDXMRW.html)
2. Copy all MATLAB function files (including *vert2con*) and associated GAMS files to a single folder
3. Run MATLAB and set the folder containing all MATLAB function files and associated GAMS files to the current folder

Details on the inputs and outputs for each MATLAB function is provided within the function file itself (in comments).
