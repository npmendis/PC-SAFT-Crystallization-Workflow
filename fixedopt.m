function [solution] = fixedopt(process_type,obj_type,para_X,para_sol1,sol1_database,para_sol2,sol2_database,M_site,Tm,DHf,P,P_lb,P_ub,T,T_lb,T_ub,molefrac_X_in,yield_lb,rho_guess,multistart_count)
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
% This function returns a list of locally optimal solutions to the 
% crystallization process optimization problem,  where the process 
% operating conditions are optimized for known solvent(s) 
% (fixed-solvent optimization problem).
%
% This process involves the crystallization of compound X. Refer to the 
% main manuscript for the general process configuration. In this function 
% file, the compounds and the streams in the process are numbered as 
% follows:
% For the compounds,
% X-'1', Solvent 1-'2’, Solvent 2-‘3’;
% For the process streams, 
% Stream S1-'1', Stream S2-'2', Stream S3-'3’, Stream S4-‘4’, Stream
% S5-'5'.
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
% process_type: a value that indicates the crystallization type and nature
% of the inputs
%   process_type = 1, a single-stage cooling crystallizer with a filter,
%   free variables: Solvent 1 type, Stream S1 composition; 
%   process_type = 21, a single-stage antisolvent crystallizer with a
%   filter, free variables: Solvent 1 type, Solvent 2 type, Stream S1
%   composition, Stream S2 flow rate;
%   process_type = 22, a single-stage antisolvent crystallizer with a
%   filter, free variables: Solvent 2 type, S2 Stream flow rate, (Solvent 1
%   type and Stream S1 composition is known);
%   process type = 3, a single-stage evaporative crystallizer with a
%   filter, free variables: Solvent 1 type, Stream S1 composition, 
%   crystallizer temperature, crystallizer pressure;.
%-------------------------------------------------------------------------- 
% obj_type: a number that indicates the objective function type
%   obj_type = 1, maximize crystal yield;
%   obj_type = 2, miximize the crystal yield normalized by the total
%   solvent consumption;
%   obj_type = 3, minimize the energy consumption normalized by the crystal
%   yield.
%
% *Note 1: 'obj_type = 3' option is only available when process_type = 3.   
%-------------------------------------------------------------------------- 
% para_X: PC-SAFT pure component paramter-vector of X
%
%   Example 1: para_X = [1.79 4.60 485.04 1938.26 0.010]
%
% *Note 1: The columns of the vector are in the following order: 
% Column 1-segment number, Column 2-segment diameter, 
% Column 3-dispersion energy
% paramter, Column 4-association energy, and Column 5-association volume.
%--------------------------------------------------------------------------
% para_sol1/para_sol2: PC-SAFT pure component paramter-vector of Solvent
% 1/Solvent 2 (which may be hypothetical)
%
%   Example 1: for a polar solvent,
%   para_sol1 = [4.38 3.68 256.56 2578.77 0.003]
%   para_sol2 = [4.38 3.68 256.56 2578.77 0.003]
%
%   Example 2: for a nonpolar solvent,
%   para_sol1 = [3.65 3.87 228.92]
%   para_sol2 = [3.65 3.87 228.92]
%
% *Note 1: The columns of the vector are in the following order: 
% Column 1-segment number, Column 2-segment diameter, 
% Column 3-dispersion energy
% paramter, Column 4-association energy, and Column 5-association volume.
%-------------------------------------------------------------------------- 
% sol1_database/sol2_database: PC-SAFT pure component parameters of the 
% Solvent 1/Solvent 2 candidates (in matrix form).
%
%   Example 1: for a database of polar solvents,
%                   [4.38	3.68	256.56	2578.77	0.003;
%                   3.43	3.53	261.59	2493.54	0.002;
%   sol1_database = 3.53	3.48	316.94	1822.33	0.009;
%                   ----	----	------	------	-----;
%                   ----	----	------	------	-----]
%                   [4.38	3.68	256.56	2578.77	0.003;
%                   3.43	3.53	261.59	2493.54	0.002;
%   sol2_database = 3.53	3.48	316.94	1822.33	0.009;
%                   ----	----	------	------	-----;
%                   ----	----	------	------	-----]
%
%   Example 2: for a database of nonpolar solvents,
%                   [3.00	4.03	312.58;
%                   3.27	3.88	287.12;
%   sol1_database = 2.81	3.72	285.69;
%                   ----	----	------;
%                   ----	----	------]
%                   [3.00	4.03	312.58;
%                   3.27	3.88	287.12;
%   sol2_database = 2.81	3.72	285.69;
%                   ----	----	------;
%                   ----	----	------]
%
% *Note 1: Rows and columns of the matrix represent solvent candidates and
% their PC-SAFT pure component parameters, respectively. Pure component
% parameters are in the same order as described under 'para_X'.
%
% *Note 2: When process_type = 1, 21, or 3, sol1_database should include 
% all the solvent candidates under consideration. When process_type = 22, 
% sol1_database can be left empty (e.g., sol1_database = []).
%
% *Note 3: The input 'sol2_database' is only required when 
% process_type = 21 or 22. For other process types, this input
% can be left empty (e.g., sol2_database = []);
%--------------------------------------------------------------------------  
% M_site: The total number of association sites (the sum of donors and 
% acceptors) of a compound (e.g., if the association scheme is 1/1, 
% M_site = 2, if the association scheme is 2/2, M_site =4, etc.). The input
% is given in vector form. See the example below.
%
%   Example 1: consider a case where X, Solvent 1, and Solvent 2 have 
%   association schemes of 2/2, 1/1, and 0/0 (nonpolar). Then,
%   M_site = [4 2 0]
%
% *Note 1: The first column represents X (Compound '1'), the second column
% represents Solvent 1 (Compound '2'), and third column represents 
% Solvent 2 (Compound '3'). If Solvent 2 is not present, leave the 
% corresponding column empty or assign any random value.
%-------------------------------------------------------------------------- 
% Tm/DHf: The value of melting temperature/enthalpy of fusion of X in 
% K/Jmol-1.
%
%   Example 1: Tm = 416.15, DHf = 29800
%-------------------------------------------------------------------------- 
% T/P: The stream and the crystallizer pressure/temperature values in 
% K/bar.
%
%   Example 1: for process_type 21, T = [298.15 298.15 298.15], P = [1 1 1]
%
% *Note 1: The first column represents Stream S1, the second column
% represents Stream S2, and third column represents the crystallizer. If
% Stream S2 is not present, assign any random value to the corresponding 
% column. See Example 2.
%
%   Example 2: for process_type 1, T = [298.15 0 298.15], P = [1 0 1]
%
% *Note 2: the stream and the crystallizer pressure/temperature values are 
% always fixed except for the case process_type = 3. When process_type = 3, 
% these inputs are the initial guesses for the stream and the crystallizer 
% pressure/temperature values.
%--------------------------------------------------------------------------
% T_lb/T_ub/P_lb/P_ub: Lower and upper bounds of stream and 
% crystallizer pressure/temperature values in bar.
%
%   Example 1. T_lb = [298.15 0 298.15], T_ub = [298.15 0 323,15], 
%   P_lb = [1 0 0.1], P_ub = [1 0 1]
%
% *Note 1: These inputs are only needed when process_type = 3. For all the 
% other instances, these inputs can be left empty (e.g., T_lb=[]).
%
% *Note 2: Lower bound = upper bound implies that particular variable is
% fixed.
%
% *Note 3: To get a meaningful result when process_type = 3, have 1) the 
% temperature and the pressure of Stream S1 fixed, 2) crystallizer 
% temperature higher than that of Stream S1, and 3) crystallizer pressure 
% lower than that of Stream S1 by setting lower and upper bounds
% appropriately as shown in the Example 1. In Example 1 bounds are imposed
% in such a way that the crystallizer temperature is always higher than
% that of Stream S1, and the crystallizer pressure is always lower than
% that of Stream S1.
%--------------------------------------------------------------------------
% molefrac_X_in: Mole fraction composition of X in Stream S1.
%
%   Example 1. molefrac_X_in=0.1
%
% *Note 1:This input is only needed when process_type = 22. 
% On the other occations, it can be left empty (e.g., molefrac_X_in=[]).
%--------------------------------------------------------------------------
% yield_lb: Lower bound for the crystal yield in mol/s (Note that the feed
% rate of X to the process is 1 mol/s).
%
%   Example 1. yield_lb=0.8
%
% *Note 1:This input is only needed when process_type = 3. 
% On the other occations, it can be left empty (e.g., yield_lb=[]).
%--------------------------------------------------------------------------
% rho_guess: An initial guess for the number density of pure X in number of
% molecules per cubic angstrom. The typical range is 0.001-0.01. We
% recommend using X's number density value obtained from the function 
% 'paraest' (in the parameter estimation step).
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
% m/sigma/epsilon/bi_eps/kappa: stand for segment number/segment
% diameter/dispersion energy parameter/association energy/association
% volume. These variables are properties of the compound type.
%
%   Example 1: if the field 'm' contains following data, 
%   [1,1.79;2,3.46;3,2.66],
%   the segment number values of X (Compound '1'), Solvent 1 (Compound
%   '2'), and Solvent 2 (Compound '3') are 1.79, 3.46, and 2.66,
%   respectively.
%--------------------------------------------------------------------------
% temperature/pressure/rho/F: stand for temperature (K)/pressure (bar)/
% number density (number of molecules per cubic angstrom)/flow rate (mol/s)
% of each stream. These variables are properties of the stream type.
%
%   Example 1: if the field 'F' contains following data, 
%   [1,5.28;3,4.28;4,0.10]
%   the flow rates of Stream S1 (Stream '1'), Stream S2 (Stream '2'),
%   Stream S3 (Stream '3'), Stream S4 (Stream '4'), Stream S5 (Stream '5')
%   are 5.28, 0, 4.28, 4.10, and 0 mol/s, respectively.
%
%   *Note 1: The flow rates of Stream S1 and Stream S5 cannot be seen in
%   the mtrix, indicating those values are zero.
%--------------------------------------------------------------------------
% rho_pure/x: pure component number density (number of molecules per cubic 
% angstrom)/mole fraction of each compound in each stream. These variables
% are properties of both the compound type and the stream type.
%
%   Example 1: if the filed 'x' contains following data,
%   [1,1,0.19;
%    1,3,2.38e-07;
%    1,4,1;
%    2,1,0.81;
%    2,3,1.00;
%    3,2,1],
%   Column 1, Column 2, and Column 3 indicates the compound type, stream 
%   type, and the mole fraction value, respectively. For instance, the
%   first row tells that the mole fraction of X (Compound '1') in Stream S1
%   (Stream '1') is 0.19 (see the numbering of compounds and streams in
%   'Function Overview')
%--------------------------------------------------------------------------
% Q: crystallizer heating duty in J/s 
% (calculated only when process_type =3).
%--------------------------------------------------------------------------
% deviation_sol1/deviation_sol1: Approximated objective function-deviations 
% for each solvent candidate in the Solvent 1/Solvent 2 database. This
% output is essential when para_sol1/para_sol2 is hypothetical.
%
%   Example 1: if the filed 'deviation_sol1' contains following data,
%   [1,0.18;
%    2,0.15;
%    3,0.13;
%    -,----;
%    -,----;]
%   Column 1 and Column 2 indicates the solvent candidate and its
%   deviation, repectively. Solvent candidates are in the same order as in
%   the input 'sol1_database'. For instance, if the first solvent candidate
%   in 'sol1_database' is benzyl alcohol, the approaximated objective
%   function-deviation for benzyl alcohol is 0.18.
%
%--------------------------------------------------------------------------
% objective: objective function value of the optimization problem. Units
% are mol/s, s-1, and MJ/mol when obj_type = 1, 2, and 3, respectively.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

solute_para.name='solute_para';
solute_para.type='parameter';
solute_para.val=para_X;
solute_para.form='full';

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

if process_type==3
    pressure_lb.name='pressure_lb';
    pressure_lb.type='parameter';
    pressure_lb.val=P_lb;
    pressure_lb.form='full';

    pressure_ub.name='pressure_ub';
    pressure_ub.type='parameter';
    pressure_ub.val=P_ub;
    pressure_ub.form='full';

    temperature_lb.name='temperature_lb';
    temperature_lb.type='parameter';
    temperature_lb.val=T_lb;
    temperature_lb.form='full';

    temperature_ub.name='temperature_ub';
    temperature_ub.type='parameter';
    temperature_ub.val=T_ub;
    temperature_ub.form='full';

    y_lb.name='y_lb';
    y_lb.type='parameter';
    y_lb.val=yield_lb;
    y_lb.form='full';   
end

if process_type==22
    mole_frac_in.name='mole_frac_in';
    mole_frac_in.type='parameter';
    mole_frac_in.val=molefrac_X_in;
    mole_frac_in.form='full';   
end

if process_type==21 || process_type==1 || process_type==3
    no_sol1.name='no_sol1';
    no_sol1.uels={1:1:length(sol1_database)};

    sol1s.name='sol1s';
    sol1s.type='parameter';
    sol1s.val=sol1_database;
    sol1s.form='full';
end

if process_type==21 || process_type==22
    no_sol2.name='no_sol2';
    no_sol2.uels={1:1:length(sol2_database)};

    sol2s.name='sol2s';
    sol2s.type='parameter';
    sol2s.val=sol2_database;
    sol2s.form='full';
end

density_guess.name='density_guess';
density_guess.type='parameter';
density_guess.val=rho_guess;
density_guess.form='full';

obj.name='obj';
obj.type='parameter';
obj.val=obj_type;
obj.form='full';

%--------------------------------------------------------------------------

j=0;
solution=[];

for i=1:1:multistart_count
    
    if process_type==21
        pcsaft_para_sol=[0 0 0 0 0; para_sol1; para_sol2];
    elseif process_type==22
        pcsaft_para_sol=[0 0 0 0 0;para_sol1; para_sol2];
    else
        pcsaft_para_sol=[0 0 0 0 0; para_sol1; 0 0 0 0 0];   
    end
    
    para_sol.name='para_sol';
    para_sol.type='parameter';
    para_sol.val=pcsaft_para_sol;
    para_sol.form='full';
    
    if process_type==21
        wgdx('fixedopt_in', solute_para, assoc_site, melting_prop1, melting_prop2, pressure, temperature, para_sol, no_sol1, sol1s, no_sol2, sol2s, density_guess,obj);
        [out1, out2, out3, out4, out5, out6, out7, out8, out9, out10, out11, out12]=gams('fixedopt-type21');

        if (out1.val==1 && out2.val==0)
            j=j+1;
            solution(j).objective=out3.val;
            solution(j).temperature=out4.val;
            solution(j).pressure=out5.val;
            solution(j).rho_pure=out6.val;
            solution(j).rho=out7.val;
            solution(j).F=out8.val;
            solution(j).x=out9.val;
            solution(j).Q=out10.val;
            solution(j).deviation_sol1=out11.val;
            solution(j).deviation_sol2=out12.val;
        end
       
    elseif process_type==22
        wgdx('fixedopt_in', solute_para, assoc_site, melting_prop1, melting_prop2, pressure, temperature, mole_frac_in, para_sol, no_sol2, sol2s, density_guess,obj);
        [out1, out2, out3, out4, out5, out6, out7, out8]=gams('fixedopt-type22');

        if (out1.val==1 && out2.val==0)
            j=j+1;
            solution(j).objective=out3.val;
            solution(j).rho_pure=out4.val;
            solution(j).rho=out5.val;
            solution(j).F=out6.val;
            solution(j).x=out7.val;
            solution(j).deviation_sol2=out8.val;
        end

    elseif process_type==1
        wgdx('fixedopt_in', solute_para, assoc_site, melting_prop1, melting_prop2, pressure, temperature, para_sol, no_sol1, sol1s, density_guess,obj);
        [out1, out2, out3, out4, out5, out6, out7, out8, out9, out10, out11]=gams('fixedopt-type1');

        if (out1.val==1 && out2.val==0)
            j=j+1;
            solution(j).objective=out3.val;
            solution(j).temperature=out4.val;
            solution(j).pressure=out5.val;
            solution(j).rho_pure=out6.val;
            solution(j).rho=out7.val;
            solution(j).F=out8.val;
            solution(j).x=out9.val;
            solution(j).Q=out10.val;
            solution(j).deviation_sol1=out11.val;
        end

    elseif process_type==3
        wgdx('fixedopt_in', solute_para, assoc_site, melting_prop1, melting_prop2, pressure, pressure_lb, pressure_ub, temperature, temperature_lb, temperature_ub, y_lb, para_sol, no_sol1, sol1s, density_guess,obj);
        [out1, out2, out3, out4, out5, out6, out7, out8 out9 out10 out11]=gams('fixedopt-type3');

        if (out1.val==1 && out2.val==0)
            j=j+1;
            solution(j).objective=out3.val;
            solution(j).temperature=out4.val;
            solution(j).pressure=out5.val;
            solution(j).rho_pure=out6.val;
            solution(j).rho=out7.val;
            solution(j).F=out8.val;
            solution(j).x=out9.val;
            solution(j).Q=out10.val;
            solution(j).deviation_sol1=out11.val;
        end
    end
       
    if (j==0)
        fprintf(1,'multistart iterations = %i, solutions = %i\n',i,j);    
    else
        fprintf(1,'multistart iterations = %i, solutions = %i, last solution = %10.5f\n',i,j,solution(j).objective); 
    end

end


end

