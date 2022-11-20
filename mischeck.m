function [solution] = mischeck(pcsaft_para_sol1,pcsaft_para_sol2,M_site,P,T,F,x_sol1,x_sol2,rho_guess_sol1,rho_guess_sol2,multistart_count)
%
% Associated Content:
% Nethrue Pramuditha Mendis, Jiayuan Wang, Richard Lakerveld, "A Workflow 
% for Crystallization Process Design with Simultaneous Process Optimization 
% and Solvent Selection based on the Perturbed-Chain Statistical Associating 
% Fluid Theory", Chemie Ingenieur Technik, 2022
%
%--------------------------------------------------------------------------
% Function Overview
%--------------------------------------------------------------------------
% This function checks the miscibility of two solvents in a stream of known
% composition
%
% Note that all the inputs to this function has to be consistent with the
% previous steps of the workflow.
%
% The PC-SAFT model equations are based on the following sources:
% 1. J. Gross and G. Sadowski, “Perturbed-Chain SAFT:  An Equation of State
% Based on a Perturbation Theory for Chain Molecules,” Ind. Eng. Chem. Res.
% , vol. 40, no. 4, pp. 1244–1260, 2001.
% 2. J. Gross and G. Sadowski, “Application of the Perturbed-Chain SAFT 
% Equation of State to Associating Systems,” Ind. Eng. Chem. Res., vol. 41,
% no. 22, pp. 5510–5515, 2002.
% 3. W. G. Chapman, K. E. Gubbins, G. Jackson, and M. Radosz, “New 
% Reference Equation of State for Associating Liquids,” Ind. Eng. Chem. 
% Res., vol. 29, no. 8, pp. 1709–1721, 1990.
% 4, Stanley H. Huang and Maciej Radosz, “Equation of State for Small, 
% Large, Polydisperse, and Associating Molecules: Extension to Fluid 
% Mixtures,” Ind. Eng. Chem. Res., vol. 30, pp. 1994–2005, 1991.
%--------------------------------------------------------------------------
% Function Inputs
%-------------------------------------------------------------------------- 
% pcsaft_para_sol1/pcsaft_para_sol2: PC-SAFT pure component parameters of
% Solvent 1/Solvent 2
%
%   Example 1:
%   pcsaft_para_sol1 =  [3.53	3.48	316.94	1822.33	0.009]
%
% *Note 2: The columns of the vector are in the following order: 
% Column 1-segment number, Column 2-segment diameter, 
% Column 3-dispersion energy
% paramter, Column 4-association energy, and Column 5-association volume.
%--------------------------------------------------------------------------  
% M_site: The total number of association sites (the sum of donors and 
% acceptors) of a compound (e.g., if the association scheme is 1/1, 
% M_site = 2, if the association scheme is 2/2, M_site =4, etc.). The input
% is given in vector form. See the example below.
%
%   Example 1: if Solvent 1 is polar and Solvent 2 is nonpolar,
%   M_site = [2 0]
%-------------------------------------------------------------------------- 
% T/P: pressure/temperature of the stream in K/bar
%
%   Example 1: T = 298.15, P = 1
%-------------------------------------------------------------------------- 
% F: Total flow rate of the stream in mols-1
%
% Example 1: F = 7.23
%-------------------------------------------------------------------------- 
% x_sol1/x_sol2: mole fractions of Solvent 1/Solvent 2 in the stream
%
% Example 1: x_sol1 = 0.85, x_sol2 =0.1
%
% *Note 1: The mole fractions here are the actuals ones (not the solute (X)
% free values)
%--------------------------------------------------------------------------
% rho_guess_sol1/rho_guess_sol2: An initial guesses for the number density 
% of pure Solvent 1/Solvent 2 in number of molecules per cubic angstrom. 
% The typical range is 0.001-0.01. We recommend using the values obtained 
% for 'rho_pure' in 'fixedopt'.
%--------------------------------------------------------------------------
% multistart_count: The number of multistart attempts. A typical range is 
% 100-1000.
%
%   Example 1. multistart_count=1000
%--------------------------------------------------------------------------
% Function Outputs
%--------------------------------------------------------------------------
% This function returns the output 'solution' in the structure form. Each 
% row corresponds to a locally optimal solution, and the fields (columns) 
% represent some selected variables.
%--------------------------------------------------------------------------
% TPDF: Tangent plane distance function value. Solvents are immiscible if
% this value is smaller than zero.
%--------------------------------------------------------------------------
% rho_pure: pure component number density (number of molecules per cubic 
% angstrom) of Solvent 1/Solvent 2
%--------------------------------------------------------------------------
% rho: number density (number of molecules per cubic angstrom) of the
% solute (X) free mixture of Solvent 1 and Solvent 2
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

assoc_site.name='assoc_site';
assoc_site.type='parameter';
assoc_site.val=M_site;
assoc_site.form='full';

pressure.name='pressure';
pressure.type='parameter';
pressure.val=P;
pressure.form='full';

temperature.name='temperature';
temperature.type='parameter';
temperature.val=T;
temperature.form='full';

mole_frac_sol1.name='mole_frac_sol1';
mole_frac_sol1.type='parameter';
mole_frac_sol1.val=F*x_sol1/(F*x_sol1+F*x_sol2);
mole_frac_sol1.form='full';   

para_sol.name='para_sol';
para_sol.type='parameter';
para_sol.val=[pcsaft_para_sol1; pcsaft_para_sol2];
para_sol.form='full';

density_guess1.name='density_guess1';
density_guess1.type='parameter';
density_guess1.val=rho_guess_sol1;
density_guess1.form='full';

density_guess2.name='density_guess2';
density_guess2.type='parameter';
density_guess2.val=rho_guess_sol2;
density_guess2.form='full';

%--------------------------------------------------------------------------

j=0;
solution=[];

for i=1:1:multistart_count
    
    wgdx('misc_check_in', assoc_site, pressure, temperature, para_sol, mole_frac_sol1, density_guess1, density_guess2);
    [out1, out2, out3, out4, out5]=gams('miscibility_check');

    if (out1.val==1 && out2.val==0)
        j=j+1;
        solution(j).TPDF=out3.val;
        solution(j).rho_pure=out4.val;
        solution(j).rho=out5.val;
    end
       
    if (j==0)
        fprintf(1,'multistart iterations = %i, solutions = %i\n',i,j);    
    else
        fprintf(1,'multistart iterations = %i, solutions = %i, last solution = %10.5f\n',i,j,solution(j).TPDF); 
    end

end


end

