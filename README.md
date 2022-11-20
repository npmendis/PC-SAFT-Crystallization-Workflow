# PC-SAFT-Crystallization-Workflow
Associated Content: Nethrue Pramuditha Mendis, Jiayuan Wang, Richard Lakerveld, "A Workflow for Crystallization Process Design with Simultaneous Process Optimization and Solvent Selection based on the Perturbed-Chain Statistical Associating Fluid Theory", Chemie Ingenieur Technik, 2022

Prerequisites:
1. MATLAB (The MathWorks, Inc.)
2. GAMS (GAMS Development Corporation) with CONOPT solver license
3. MATLAB function file *vert2con* (can be downloaded at https://www.mathworks.com/matlabcentral/fileexchange/7895-vert2con-vertices-to-constraints)

Provided in the git repository:
1. The MATLAB function files *paraest*, *relxopt*, *fixedopt*, and *mischeck* with associated GAMS files
2. PC-SAFT parameters of the polar and nonpolar solvents used in the manuscript (MAT files)

Setting up:
1. Interface GAMS and MATLAB with GDXMRW procedure (details can be found at https://www.gams.com/latest/docs/T_GDXMRW.html)
2. Copy all MATLAB function files (including *vert2con*) and associated GAMS files in the git repository to a single folder
3. Run MATLAB and set the folder containing all MATLAB function files and associated GAMS files to the current folder

Details on the inputs and outputs for each MATLAB function is provided within the function file itself (in comments).
