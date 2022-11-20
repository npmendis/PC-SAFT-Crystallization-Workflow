function [solution] = paraest(pcsaft_para_sol,M_site,Tm,DHf,P,T,x_sol,para_X_lo,para_X_up,multistart_count)
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
% This function returns a list of locally optimal solutions to the PC-SAFT
% pure component parameter estimation problem of X.
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
% pcsaft_para_sol: PC-SAFT pure component parameters of the 
% reference solvents (in matrix form).
%
%   Example 1:
%                       [4.38	3.68	256.56	2578.77	0.003;
%                       3.43	3.53	261.59	2493.54	0.002;
%   pcsaft_para_sol =   3.53	3.48	316.94	1822.33	0.009;
%                       3.00	4.03	312.58	  0  	  0;
%                       3.27	3.88	287.12	  0	      0]
%
% *Note 1: Rows and columns of the matrix represent reference solvents and
% their PC-SAFT pure component parameters, respectively.
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
%   Example 1: for reference solvents ethanol, ethyl acetate, acetone, 
%   acetic acid, toluene, cyclohexane, and 1,4-dioxane, whose association
%   schemes are 1/1, 1/1, 1/1, 1/1, 0/0, 0/0, and 1/1, respectively. Then,
%   M_site = [2 2 2 2 0 0 2]
%-------------------------------------------------------------------------- 
% Tm/DHf: The value of melting temperature/enthalpy of fusion of X in 
% K/Jmol-1.
%
%   Example 1: Tm = 416.15, DHf = 29800
%-------------------------------------------------------------------------- 
% T/P: pressure/temperature at which the solubility is measured/predicted
% in K/bar
%
%   Example 1: T = 298.15, P = 1
%--------------------------------------------------------------------------
% x_sol: Mole fraction solubility of X in each reference 
%
%   Example 1: for reference solvents ethanol, ethyl acetate, acetone, 
%   acetic acid, toluene, cyclohexane, and 1,4-dioxane
%   x_sol = [0.0855,0.0448,0.0828,0.0435,0.00129,2.34e-05,0.0516]
%--------------------------------------------------------------------------
% para_X_lo: lower bounds of the PC-SFAT pure component paramters of X
%
%   Example 1: para_X_lo = [1 1 150 500 0.010]
%
% *Note 1: The columns of the vector are in the following order: 
% Column 1-segment number, Column 2-segment diameter, 
% Column 3-dispersion energy
% paramter, Column 4-association energy, and Column 5-association volume.
%--------------------------------------------------------------------------
% para_X_up: upper bounds of the PC-SFAT pure component paramters of X
%
%   Example 1: para_X_up = [7 7 500 4000 0.010]
%
% *Note 1: The columns of the vector are in the following order: 
% Column 1-segment number, Column 2-segment diameter, 
% Column 3-dispersion energy
% paramter, Column 4-association energy, and Column 5-association volume.
% *Note 2: If a certain parameter needs to be fixed, choose the same lower
% and upper bounds for that parameter.
%--------------------------------------------------------------------------
% multistart_count: The number of multistart attempts. A typical range is 
% 100-1000.
%
%   Example 1. multistart_count=100
%--------------------------------------------------------------------------
% Function Outputs
%--------------------------------------------------------------------------
% This function returns the output 'solution' in the structure form. Each 
% row corresponds to a locally optimal solution, and the fields (columns) 
% represent some selected variables.
%--------------------------------------------------------------------------
% objective: objective function value of the fitting problem. 
%--------------------------------------------------------------------------
% rho/ln_acti_coeff_exp/ln_acti_coeff: number density (number of molecules 
% per cubic angstrom) of each saturated solution of X/experimental 
% activity coefficient (logarithm) at the experimental saturation point/
% predicted activity coefficient (logarithm) at the experimental saturation
% point
%
%   Example 1:  for reference solvents ethanol, ethyl acetate, acetone, 
%   acetic acid, toluene, cyclohexane, and 1,4-dioxane, if the filed 
%   'ln_acti_coeff_exp' contains following data,
%   [1,-0.9494
%    2,-0.3031
%    3,-0.9173
%    4,-0.2736
%    5,3.24449
%    6,7.25415
%    7,-0.4444]
%    Column 1 is the reference solvent (in the given order) and Column 2 is
%    the logarithm of the activity coefficient of X.
%
% *Note 1: Compare ln_acti_coeff_exp with ln_acti_coeff to assess the
% quality of the fit
%--------------------------------------------------------------------------
% m/sigma/epsilon/bi_eps/kappa/rho_X: stand for the segment number/segment
% diameter/dispersion energy parameter/association energy/association
% volume/pure component number density (number of molecules 
% per cubic angstrom) of X.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

comp.name='comp';
comp.uels={1:1:length(pcsaft_para_sol)+1};

mixture.name='mixture';
mixture.uels={1:1:length(pcsaft_para_sol)};

assoc_site.name='assoc_site';
assoc_site.type='parameter';
assoc_site.val=M_site;
assoc_site.form='full';

melting_prop1.name='melting_prop1';
melting_prop1.type='parameter';
melting_prop1.val=Tm;
melting_prop1.form='full';

melting_prop2.name='melting_prop2';
melting_prop2.type='parameter';
melting_prop2.val=DHf;
melting_prop2.form='full';

pressure.name='pressure';
pressure.type='parameter';
pressure.val=P;
pressure.form='full';

temperature.name='temperature';
temperature.type='parameter';
temperature.val=T;
temperature.form='full';

xl=zeros(length(pcsaft_para_sol)+1,length(pcsaft_para_sol));
xl(1,:)=x_sol;
for i=1:1:length(pcsaft_para_sol)
    xl(i+1,i)=1-xl(1,i);
end
mf_l.name='mf_l';
mf_l.type='parameter';
mf_l.val=xl;
mf_l.form='full';

R=8.31446;
lngamma=DHf/R*(1/Tm-1/T)-log(x_sol);
ln_act_coeff.name='ln_act_coeff';
ln_act_coeff.type='parameter';
ln_act_coeff.val=lngamma;
ln_act_coeff.form='full';

pcsaft_para_lo=[para_X_lo; pcsaft_para_sol];
para_lo.name='para_lo';
para_lo.type='parameter';
para_lo.val=pcsaft_para_lo;
para_lo.form='full';

pcsaft_para_up=[para_X_up; pcsaft_para_sol];
para_up.name='para_up';
para_up.type='parameter';
para_up.val=pcsaft_para_up;
para_up.form='full';

%--------------------------------------------------------------------------

j=0;
solution=[];

for i=1:1:multistart_count

    pcsaft_para_l=[[1+(3-1)*rand, 1+(3-1)*rand, 150+(500-150)*rand, 1000+(3000-1000)*rand, 0.01]; pcsaft_para_sol];
    para_l.name='para_l';
    para_l.type='parameter';
    para_l.val=pcsaft_para_l;
    para_l.form='full';

    rho_guess=rand;
    while rho_guess<=0.1 || rho_guess>=0.7
        rho_guess=rand;
    end
    density_guess.name='density_guess';
    density_guess.type='parameter';
    density_guess.val=rho_guess;
    density_guess.form='full';

    wgdx('para_est', comp, mixture, assoc_site, melting_prop1, melting_prop2, pressure, temperature, mf_l, ln_act_coeff, para_l, para_lo, para_up, density_guess);

    [out1, out2, out3, out4, out5, out6, out7, out8, out9, out10, out11, out12]=gams('parameter_estimation2');
    
    if (out1.val==1 && out2.val==0)
       j=j+1;
       solution(j).objective=out3.val;
       solution(j).rho_X=out4.val;
       solution(j).rho=out5.val;
       solution(j).m=out6.val(1,2);
       solution(j).sigma=out7.val(1,2);
       solution(j).epsilon=out8.val(1,2);
       solution(j).bi_eps=out9.val(1,2);
       solution(j).kappa=out10.val(1,2);
       solution(j).ln_acti_coeff_exp=out11.val;
       solution(j).lnacti_coeff=out12.val;
    end
       
    if (j==0)
        fprintf(1,'multistart iterations = %i, solutions = %i\n',i,j);    
    else
        fprintf(1,'multistart iterations = %i, solutions = %i, last solution = %10.5f\n',i,j,solution(j).objective); 
    end

end



end

