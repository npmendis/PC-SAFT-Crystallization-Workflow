$set matout "'misc_check_out.gdx', solver_status, infes, objval, rho_pure, rho"

* Associated Content:
* Nethrue Pramuditha Mendis, Jiayuan Wang, Richard Lakerveld, "A Workflow for Crystallization Process Design with Simultaneous Process Optimization and Solvent Selection based on the Perturbed-Chain Statistical Associating Fluid Theory", Chemie Ingenieur Technik, 2022

*GAMS model for the miscibility check of two solvents

*The PC-SAFT part of the model is based on:
* 1. J. Gross and G. Sadowski, �Perturbed-Chain SAFT:  An Equation of State Based on a Perturbation Theory for Chain Molecules,� Ind. Eng. Chem. Res., vol. 40, no. 4, pp. 1244�1260, 2001.
* 2. J. Gross and G. Sadowski, �Application of the Perturbed-Chain SAFT Equation of State to Associating Systems,� Ind. Eng. Chem. Res., vol. 41, no. 22, pp. 5510�5515, 2002.
* 3. W. G. Chapman, K. E. Gubbins, G. Jackson, and M. Radosz, �New Reference Equation of State for Associating Liquids,� Ind. Eng. Chem. Res., vol. 29, no. 8, pp. 1709�1721, 1990.
* 4, Stanley H. Huang and Maciej Radosz, �Equation of State for Small, Large, Polydisperse, and Associating Molecules: Extension to Fluid Mixtures,� Ind. Eng. Chem. Res., vol. 30, pp. 1994�2005, 1991.

*-------------------------------------------------------------------------------
*Definitions
*-------------------------------------------------------------------------------

Set
   i components /1*2/
   state streams /1*2/
   para pure component parameters /1*5/;
   alias (i,j,k);

Parameters

   M_site(i) Association site on component i

   P_init initial guess for operating pressure in bar
   T_init initial guess for operating temperature in K

   x_sol1_in solute inlet mole fraction

   pcsaft_para_sol_ini(i,para) initial guess for pure component parameters of the solvents

   rho_ini1 solvent1 pure component sumber density guess
   rho_ini2 solvent2 pure component sumber density guess
   ;

$gdxin misc_check_in
$load M_site=assoc_site
$load P_init=pressure
$load T_init=temperature
$load x_sol1_in=mole_frac_sol1
$load pcsaft_para_sol_ini=para_sol
$load rho_ini1=density_guess1
$load rho_ini2=density_guess2
$gdxin

Positive Variables
   P(state) operating pressure in bar
   T(state) operating temperature in K
   F(state) flow rate mols-1
   x(i,state) composition (molar fraction)
   ;

Variables
   obj
   ;

*PC-SAFT-related
*-------------------------------------------------------------------------------

Set
   n /0,1,2,3/
   es_i /0,1,2/
   es_j /0,1,2,3,4,5,6/

Parameters
   R universal constant /8.31446/
   kb boltzmann constant   /1.38064852e2/
   np(n) /'0' 0,'1' 1,'2' 2,'3' 3/
   es_jp(es_j) /'0' 0,'1' 1,'2' 2,'3' 3,'4' 4,'5' 5,'6' 6/;

Table ap(es_j,es_i)
                      0                       1                       2
       0         0.9105631445           -0.3084016918           -0.0906148351
       1         0.6361281449           0.1860531159            0.4527842806
       2         2.6861347891           -2.5030047259           0.5962700728
       3        -26.547362491           21.419793629            -1.7241829131
       4         97.759208784           -65.255885330           -4.1302112531
       5        -159.59154087           83.318680481            13.776631870
       6         91.297774084           -33.746922930           -8.6728470368 ;
Table bp(es_j,es_i)
                      0                       1                       2
       0         0.7240946941          -0.5755498075            0.0976883116
       1         2.2382791861           0.6995095521            -0.2557574982
       2         -4.0025849485          3.8925673390            -9.1558561530
       3        -21.003576815           -17.215471648           20.642075974
       4         26.855641363           192.67226447            -38.804430052
       5         206.55133841          -161.82646165            93.626774077
       6        -355.60235612           -165.20769346           -29.666905585;

Positive Variables

   bi_eps(i)  association energy for pure component
   kappa(i)   association volume for pure component
   epsilon(i) reduced potential depth for pure component
   sigma(i)   segment diameter for pure component
   m(i)       segement number for pure component

   Mw(i)       molar mass for pure component

   bi_epsij(i,j) association energy for mixtures
   kappaij(i,j)  association volume for mixtures
   d(i,state) temperature dependent segment diameter
   z(n,state)  Zeta(crystallizer)
   z_pure(n,i,state)       Zeta for the pure component(crystallizer)
   rho(state)      total number density of mixtures
   rho_pure(i,state)     total number density of pure component(crystallizer)
   dij(i,j,state)   intermedia variable

   XA(i,state) mole fraction of i not bounded at site A assume that XA=XB
   XA_pure(i,state)
   mm(state)  mean segment number in the mixture
   order1(state)      first order perturbation term
   order1_pure(i,state)     first order perturbation term for pure component(crystallizer)
   order2(state)     second order perturbation term
   order2_pure(i,state)  second order perturbation term for pure component(crystallizer)
   sigij(i,j)       parameters for a pair of unlike segments
   epsij(i,j,state) parameters for a pair of unlike segments
   Ztot(state)      total compressibility factor(crystallizer)
   Ztot_pure(i,state)       total compressibility factor for pure component(crystallizer)
   ghs(i,j,state)   hard-sphere radial distribution function
   ghs_pure(i,state)   hard-sphere radial distribution function for pure component(crystallizer)
   assoc(i,j,state) association strength
   assoc_pure(i,state)     association strength for pure component(crystallizer)

Variables

   aassoc(state)     Helmholtz free energy due to association
   aassoc_pure(i,state)    Helmholtz free energy due to association for pure component
   ahs(state)    Helmholtz free energy due to hard-sphere
   ahs_pure(i,state)        Helmholtz free energy due to hard-sphere for pure component
   ahc(state)       Helmholtz free energy due to hard-chain
   ahc_pure(i,state)         Helmholtz free energy due to hard-chain for pure component
   a(es_j,state) power series coefficient 'a'
   a_pure(es_j,i,state)
   b(es_j,state) power series coefficient 'b'
   b_pure(es_j,i,state)
   I1(state)
   I1_pure(i,state)
   I2(state)
   I2_pure(i,state)
   C1(state)
   C1_pure(i,state)
   adisp(state)       Helmholtz free energy due to dispertion
   adisp_pure(i,state)     Helmholtz free energy due to dispertion for pure component
   Zhs(state)   compressibility factor due to hard-sphere
   Zhs_pure(i,state)      compressibility factor due to hard-sphere for pure component
   rhogij(i,j,state)
   rhog_pure(i,state)
   Zhc(state)      compressibility factor due to hard-chain
   Zhc_pure(i,state)       compressibility factor due to hard-chain for pure component
   I1z(state) I1 into eta derivative
   I1z_pure(i,state)
   I2z(state) I2 into eta derivative
   I2z_pure(i,state)
   C2(state)
   C2_pure(i,state)
   Zdisp(state)    compressibility factor due to dispersion
   Zdisp_pure(i,state)     compressibility factor due to dispersion for pure component
   Zassoc(state)    comporessibility factor due to association
   Zassoc_pure(i,state)     comporessibility factor due to association for pure component
   dzdx(n,i,state) derivative of zeta
   dahsdx(i,state)
   dghsdx(i,j,k,state)
   dassocdx(i,j,k,state)
   daassocdx(k,state)
   rhoassoc(i,j,state)
   rhoassoc_pure(i,state)
   dXAdx(i,k,state)
   rhoXA(i,state)
   rhoXA_pure(i,state)
   Zassoc(state)
   Zassoc_pure(i,state)
   dahcdx(i,state)
   ord1dx(i,state)
   ord2dx(i,state)
   I1dx(i,state)
   I2dx(i,state)
   adx(es_j,i,state) derivative of power series coefficient 'a'
   bdx(es_j,i,state) derivative of power series coefficient 'b'
   C1dx(i,state)
   dadispdx(i,state)
   daresdx(i,state)
   lnphi(i,state)
   lnphi_pure(i,state)
   ares(state)
   ares_pure(i,state)
   binary_interact(i,j,state)
   a_binary(i,j)
   b_binary(i,j)
   TPDF tangent place distance function
;

*-------------------------------------------------------------------------------
*Variable initialization
*-------------------------------------------------------------------------------

execseed = 1e8*(frac(jnow));

P.fx('1')=P_init;
P.fx('2')=P_init;

T.fx('1')=T_init;
T.fx('2')=T_init;

*Flow rates
F.fx('1')=1;
F.fx('2')=1;

*Compositions
x.fx('1','1')=x_sol1_in;
x.l('2','1')=1-x.l('1','1');

x.l('1','2')=1*uniform(0,1);
x.l('2','2')=1-x.l('1','2');

*PC-SAFT-related
*-------------------------------------------------------------------------------

m.fx(i)=pcsaft_para_sol_ini(i,'1');
sigma.fx(i)=pcsaft_para_sol_ini(i,'2');
epsilon.fx(i)=pcsaft_para_sol_ini(i,'3');
bi_eps.fx(i)=pcsaft_para_sol_ini(i,'4');
kappa.fx(i)=pcsaft_para_sol_ini(i,'5');

a_binary.fx(i,j)=0;
b_binary.fx(i,j)=0;

d.l(i,state)=sigma.l(i)*(1-0.12*exp(-3*epsilon.l(i)/T.l(state)));

rho.l(state)=6/pi*0.5/sum(i,x.l(i,state)*m.l(i)*power(d.l(i,state),3));
rho.lo(state)=0.001;
rho.up(state)=0.035;

rho_pure.l(i,state)=6/pi*0.5/(m.l(i)*power(d.l(i,state),3));
rho_pure.l('1',state)=rho_ini1;
rho_pure.l('2',state)=rho_ini2;
rho_pure.lo(i,state)=0.003;
rho_pure.up(i,state)=0.035;

XA.l(i,state)=0.228;
XA.lo(i,state)=0.000001;

XA_pure.l(i,state)=0.117;
XA_pure.lo(i,state)=0.000001;

Ztot.l(state)=P.l(state)/(kb*T.l(state)*rho.l(state));
Ztot_pure.l(i,state)=P.l(state)/(kb*T.l(state)*rho_pure.l(i,state));

z.l(n,state)=pi/6*rho.l(state)*sum(i,x.l(i,state)*m.l(i)*power(d.l(i,state),np(n)));
z.up('3',state)=0.74048;

z_pure.l(n,i,state)=pi/6*rho_pure.l(i,state)*m.l(i)*power(d.l(i,state),np(n));
z_pure.up('3',i,state)=0.74048;

dij.l(i,j,state)=d.l(i,state)*d.l(j,state)/(d.l(i,state)+d.l(j,state));

ghs.l(i,j,state)=1/(1-z.l('3',state))+dij.l(i,j,state)*3*z.l('2',state)/(1-z.l('3',state))**2+dij.l(i,j,state)*dij.l(i,j,state)
                            *2*z.l('2',state)**2/(1-z.l('3',state))**3;
ghs_pure.l(i,state)=1/(1-z_pure.l('3',i,state))+d.l(i,state)*3/2*z_pure.l('2',i,state)/(1-z_pure.l('3',i,state))**2+d.l(i,state)*d.l(i,state)
                            /2*z_pure.l('2',i,state)**2/(1-z_pure.l('3',i,state))**3;

sigij.l(i,j)=0.5*(sigma.l(i)+sigma.l(j));
epsij.l(i,j,state)=(epsilon.l(i)*epsilon.l(j))**0.5;

bi_epsij.l(i,j)=bi_eps.l(i)/2+bi_eps.l(j)/2;
kappaij.l(i,j)=(kappa.l(i)*kappa.l(j))**0.5*power((2*(sigma.l(i)*sigma.l(j))**0.5/(sigma.l(i)+sigma.l(j))),3);

assoc.l(i,j,state)=(d.l(i,state)+d.l(j,state))**3/8*ghs.l(i,j,state)*kappaij.l(i,j)*(exp(bi_epsij.l(i,j)/T.l(state))-1);
assoc_pure.l(i,state)=d.l(i,state)**3*ghs_pure.l(i,state)*kappa.l(i)*(exp(bi_eps.l(i)/T.l(state))-1);

mm.l(state)=sum(i,x.l(i,state)*m.l(i));

ahs.l(state)=1/z.l('0',state)*(3*z.l('1',state)*z.l('2',state)/(1-z.l('3',state))+z.l('2',state)**3/z.l('3',state)/(1-z.l('3',state))**2+
                  (z.l('2',state)**3/z.l('3',state)**2-z.l('0',state))*log(1-z.l('3',state)));
ahs_pure.l(i,state)=1/z_pure.l('0',i,state)*(3*z_pure.l('1',i,state)*z_pure.l('2',i,state)/(1-z_pure.l('3',i,state))+z_pure.l('2',i,state)**3/z_pure.l('3',i,state)/(1-z_pure.l('3',i,state))**2+
                  (z_pure.l('2',i,state)**3/z_pure.l('3',i,state)**2-z_pure.l('0',i,state))*log(1-z_pure.l('3',i,state)));

ahc.l(state)=mm.l(state)*ahs.l(state)-sum(i,x.l(i,state)*(m.l(i)-1)*log(ghs.l(i,i,state)));
ahc_pure.l(i,state)=m.l(i)*ahs_pure.l(i,state)-(m.l(i)-1)*log(ghs_pure.l(i,state));

order1.l(state)=sum((i,j),x.l(i,state)*x.l(j,state)*m.l(i)*m.l(j)*epsij.l(i,j,state)/T.l(state)*sigij.l(i,j)**3);
order1_pure.l(i,state)=m.l(i)*m.l(i)*epsilon.l(i)/T.l(state)*sigma.l(i)**3;

order2.l(state)=sum((i,j),x.l(i,state)*x.l(j,state)*m.l(i)*m.l(j)*epsij.l(i,j,state)**2/T.l(state)**2*sigij.l(i,j)**3);
order2_pure.l(i,state)=m.l(i)*m.l(i)*epsilon.l(i)*epsilon.l(i)/T.l(state)/T.l(state)*sigma.l(i)**3;

a.l(es_j,state)=ap(es_j,'0')+(mm.l(state)-1)/mm.l(state)*ap(es_j,'1')+(mm.l(state)-1)*(mm.l(state)-2)/mm.l(state)/mm.l(state)*ap(es_j,'2');
a_pure.l(es_j,i,state)=ap(es_j,'0')+(m.l(i)-1)/m.l(i)*ap(es_j,'1')+(m.l(i)-1)*(m.l(i)-2)/m.l(i)/m.l(i)*ap(es_j,'2');

b.l(es_j,state)=bp(es_j,'0')+(mm.l(state)-1)/mm.l(state)*bp(es_j,'1')+(mm.l(state)-1)*(mm.l(state)-2)/mm.l(state)/mm.l(state)*bp(es_j,'2');
b_pure.l(es_j,i,state)=bp(es_j,'0')+(m.l(i)-1)/m.l(i)*bp(es_j,'1')+(m.l(i)-1)*(m.l(i)-2)/m.l(i)/m.l(i)*bp(es_j,'2');

I1.l(state)=sum(es_j,a.l(es_j,state)*power(z.l('3',state),es_jp(es_j)));
I1_pure.l(i,state)=sum(es_j,a_pure.l(es_j,i,state)*power(z_pure.l('3',i,state),es_jp(es_j)));

I2.l(state)=sum(es_j,b.l(es_j,state)*power(z.l('3',state),es_jp(es_j)));
I2_pure.l(i,state)=sum(es_j,b_pure.l(es_j,i,state)*power(z_pure.l('3',i,state),es_jp(es_j)));

C1.l(state)=1/(1+mm.l(state)*(8*z.l('3',state)-2*z.l('3',state)**2)/(1-z.l('3',state))**4+(1-mm.l(state))*(20*z.l('3',state)-27*z.l('3',state)**2+12*
                z.l('3',state)**3-2*z.l('3',state)**4)/(1-z.l('3',state))**2/(2-z.l('3',state))**2);
C1_pure.l(i,state)=1/(1+m.l(i)*(8*z_pure.l('3',i,state)-2*z_pure.l('3',i,state)**2)/(1-z_pure.l('3',i,state))**4+(1-m.l(i))*(20*z_pure.l('3',i,state)-27*z_pure.l('3',i,state)**2+12*
                z_pure.l('3',i,state)**3-2*z_pure.l('3',i,state)**4)/(1-z_pure.l('3',i,state))**2/(2-z_pure.l('3',i,state))**2);

adisp.l(state)=-2*pi*rho.l(state)*I1.l(state)*order1.l(state)-pi*rho.l(state)*mm.l(state)*C1.l(state)*I2.l(state)*order2.l(state);
adisp_pure.l(i,state)=-2*pi*rho_pure.l(i,state)*I1_pure.l(i,state)*order1_pure.l(i,state)-pi*rho_pure.l(i,state)*m.l(i)*C1_pure.l(i,state)*I2_pure.l(i,state)*order2_pure.l(i,state);

Zhs.l(state)=z.l('3',state)/(1-z.l('3',state))+3*z.l('1',state)*z.l('2',state)/z.l('0',state)/(1-z.l('3',state))**2+(3*z.l('2',state)**3-
                  z.l('3',state)*z.l('2',state)**3)/z.l('0',state)/(1-z.l('3',state))**3;
Zhs_pure.l(i,state)=z_pure.l('3',i,state)/(1-z_pure.l('3',i,state))+3*z_pure.l('1',i,state)*z_pure.l('2',i,state)/z_pure.l('0',i,state)/(1-z_pure.l('3',i,state))**2+(3*z_pure.l('2',i,state)**3-
                  z_pure.l('3',i,state)*z_pure.l('2',i,state)**3)/z_pure.l('0',i,state)/(1-z_pure.l('3',i,state))**3;

rhogij.l(i,j,state)=z.l('3',state)/(1-z.l('3',state))**2+dij.l(i,j,state)*(3*z.l('2',state)/(1-z.l('3',state))**2+6*
                  z.l('2',state)*z.l('3',state)/(1-z.l('3',state))**3)+dij.l(i,j,state)*dij.l(i,j,state)*(4*z.l('2',state)**2/(1-z.l('3',state))**3+6*
                  z.l('2',state)**2*z.l('3',state)/(1-z.l('3',state))**4);
rhog_pure.l(i,state)=z_pure.l('3',i,state)/(1-z_pure.l('3',i,state))**2+d.l(i,state)/2*(3*z_pure.l('2',i,state)/(1-z_pure.l('3',i,state))**2+6*
                  z_pure.l('2',i,state)*z_pure.l('3',i,state)/(1-z_pure.l('3',i,state))**3)+d.l(i,state)*d.l(i,state)/4*(4*z_pure.l('2',i,state)**2/(1-z_pure.l('3',i,state))**3+6*
                  z_pure.l('2',i,state)**2*z_pure.l('3',i,state)/(1-z_pure.l('3',i,state))**4);

Zhc.l(state)=mm.l(state)*Zhs.l(state)-sum(i,x.l(i,state)*(m.l(i)-1)/ghs.l(i,i,state)*rhogij.l(i,i,state));
Zhc_pure.l(i,state)=m.l(i)*Zhs_pure.l(i,state)-(m.l(i)-1)/ghs_pure.l(i,state)*rhog_pure.l(i,state);

I1z.l(state)=sum(es_j,a.l(es_j,state)*(es_jp(es_j)+1)*power(z.l('3',state),es_jp(es_j)));
I1z_pure.l(i,state)=sum(es_j,a_pure.l(es_j,i,state)*(es_jp(es_j)+1)*power(z_pure.l('3',i,state),es_jp(es_j)));

I2z.l(state)=sum(es_j,b.l(es_j,state)*(es_jp(es_j)+1)*power(z.l('3',state),es_jp(es_j)));
I2z_pure.l(i,state)=sum(es_j,b_pure.l(es_j,i,state)*(es_jp(es_j)+1)*power(z_pure.l('3',i,state),es_jp(es_j)));

C2.l(state)=-C1.l(state)*C1.l(state)*(mm.l(state)*(-4*z.l('3',state)**2+20*z.l('3',state)+8)/(1-z.l('3',state))**5+(1-mm.l(state))*(2*z.l('3',state)**3
                +12*z.l('3',state)**2-48*z.l('3',state)+40)/(1-z.l('3',state))**3/(2-z.l('3',state))**3);
C2_pure.l(i,state)=-C1_pure.l(i,state)*C1_pure.l(i,state)*(m.l(i)*(-4*z_pure.l('3',i,state)**2+20*z_pure.l('3',i,state)+8)/(1-z_pure.l('3',i,state))**5+(1-m.l(i))*(2*z_pure.l('3',i,state)**3
                +12*z_pure.l('3',i,state)**2-48*z_pure.l('3',i,state)+40)/(1-z_pure.l('3',i,state))**3/(2-z_pure.l('3',i,state))**3);

Zdisp.l(state)=-2*pi*rho.l(state)*I1z.l(state)*order1.l(state)-pi*rho.l(state)*mm.l(state)*(C1.l(state)*I2z.l(state)
                                  +C2.l(state)*z.l('3',state)*I2.l(state))*order2.l(state);
Zdisp_pure.l(i,state)=-2*pi*rho_pure.l(i,state)*I1z_pure.l(i,state)*order1_pure.l(i,state)-pi*rho_pure.l(i,state)*m.l(i)*(C1_pure.l(i,state)*I2z_pure.l(i,state)
                                  +C2_pure.l(i,state)*z_pure.l('3',i,state)*I2_pure.l(i,state))*order2_pure.l(i,state);

*-------------------------------------------------------------------------------
*Equations define
*-------------------------------------------------------------------------------

Equations

mole_frac_bound(state)

obj_function
;

*PC-SAFT-related
*-------------------------------------------------------------------------------

Equations

d_define(i,state)

z_define(n,state)
z_pure_define(n,i,state)

dij_define(i,j,state)

ghs_define(i,j,state)
ghs_pure_define(i,state)

binary_interact_define(i,j,state)
binary_interact_constraint1
binary_interact_constraint2

sigij_define(i,j)
epsij_define(i,j,state)
bi_epsij_define(i,j)
kappaij_define(i,j)

assoc_define(i,j,state)
assoc_pure_define(i,state)

mm_define(state)

XA_define(i,state)
XA_pure_define(i,state)

aassoc_define(state)
aassoc_pure_define(i,state)

ahs_define(state)
ahs_pure_define(i,state)

ahc_define(state)
ahc_pure_define(i,state)

order1_define(state)
order1_pure_define(i,state)

order2_define(state)
order2_pure_define(i,state)

a_define(es_j,state)
a_pure_define(es_j,i,state)

b_define(es_j,state)
b_pure_define(es_j,i,state)

ares_define(state)
ares_pure_define(i,state)

I1_define(state)
I1_pure_define(i,state)

I2_define(state)
I2_pure_define(i,state)

C1_define(state)
C1_pure_define(i,state)

adisp_define(state)
adisp_pure_define(i,state)

Zhs_define(state)
Zhs_pure_define(i,state)

rhogij_define(i,j,state)
rhog_pure_define(i,state)

Zhc_define(state)
Zhc_pure_define(i,state)

I1z_define(state)
I1z_pure_define(i,state)

I2z_define(state)
I2z_pure_define(i,state)

C2_define(state)
C2_pure_define(i,state)

Zdisp_define(state)
Zdisp_pure_define(i,state)

Ztot_define(state)
Ztot_pure_define(i,state)

PVNRT(state)
PVNRT_pure(i,state)

dzdx_define(n,i,state)

dahsdx_define(i,state)

dghsdx_define(i,j,k,state)

dassocdx_define(i,j,k,state)

rhoassoc_define(i,j,state)
rhoassoc_pure_define(i,state)

dXAdx_define(i,k,state)

rhoXA_define(i,state)
rhoXA_pure_define(i,state)

daassocdx_define(k,state)

Zassoc_define(state)
Zassoc_pure_define(i,state)

dahcdx_define(i,state)

ord1dx_define(i,state)

ord2dx_define(i,state)

adx_define(es_j,i,state)

bdx_define(es_j,i,state)

I1dx_define(i,state)

I2dx_define(i,state)

C1dx_define(i,state)

dadispdx_define(i,state)

daresdx_define(i,state)

lnphi_define(i,state)
lnphi_pure_define(i,state)

TPDF_define
;

*-------------------------------------------------------------------------------
*Equations
*-------------------------------------------------------------------------------

mole_frac_bound(state)..sum(i,x(i,state))=e=1;

TPDF_define..TPDF=e=sum(i,x(i,'2')*(lnphi(i,'2')+log(x(i,'2'))-(lnphi(i,'1')+log(x(i,'1')))));

obj_function.. obj=e=TPDF;


*PC-SAFT-related
* ------------------------------------------------------------------------------

d_define(i,state)..d(i,state)=e=sigma(i)*(1-0.12*exp(-3*epsilon(i)/T(state)));

z_define(n,state)..z(n,state)=e=pi/6*rho(state)*sum(i,x(i,state)*m(i)*power(d(i,state),np(n)));
z_pure_define(n,i,state)..z_pure(n,i,state)=e=pi/6*rho_pure(i,state)*m(i)*power(d(i,state),np(n));

dij_define(i,j,state)..dij(i,j,state)=e=d(i,state)*d(j,state)/(d(i,state)+d(j,state));

ghs_define(i,j,state)..ghs(i,j,state)=e=1/(1-z('3',state))+dij(i,j,state)*3*z('2',state)/(1-z('3',state))**2+dij(i,j,state)*dij(i,j,state)
                            *2*z('2',state)**2/(1-z('3',state))**3;
ghs_pure_define(i,state)..ghs_pure(i,state)=e=1/(1-z_pure('3',i,state))+d(i,state)*3/2*z_pure('2',i,state)/(1-z_pure('3',i,state))**2+d(i,state)*d(i,state)
                            /2*z_pure('2',i,state)**2/(1-z_pure('3',i,state))**3;

binary_interact_define(i,j,state)..binary_interact(i,j,state)=e=a_binary(i,j)+b_binary(i,j)*T(state);
binary_interact_constraint1..a_binary('1','2')=e=a_binary('2','1');
binary_interact_constraint2..b_binary('1','2')=e=b_binary('2','1');

sigij_define(i,j)..sigij(i,j)=e=0.5*(sigma(i)+sigma(j));
epsij_define(i,j,state)..epsij(i,j,state)=e=(1-binary_interact(i,j,state))*(epsilon(i)*epsilon(j))**0.5;

bi_epsij_define(i,j)..bi_epsij(i,j)=e=bi_eps(i)/2+bi_eps(j)/2;
kappaij_define(i,j)..kappaij(i,j)=e=(kappa(i)*kappa(j))**0.5*power((2*(sigma(i)*sigma(j))**0.5/(sigma(i)+sigma(j))),3);

assoc_define(i,j,state)..assoc(i,j,state)=e=(d(i,state)+d(j,state))**3/8*ghs(i,j,state)*kappaij(i,j)*(exp(bi_epsij(i,j)/T(state))-1);
assoc_pure_define(i,state)..assoc_pure(i,state)=e=d(i,state)**3*ghs_pure(i,state)*kappa(i)*(exp(bi_eps(i)/T(state))-1);

mm_define(state)..mm(state)=e=sum(i,x(i,state)*m(i));

XA_define(i,state)..1/XA(i,state)=e=1+sum(j,x(j,state)*rho(state)*XA(j,state)*M_site(j)/2*assoc(i,j,state));
XA_pure_define(i,state)..1=e=XA_pure(i,state) +rho_pure(i,state)*XA_pure(i,state) *XA_pure(i,state) *M_site(i)/2*assoc_pure(i,state);

aassoc_define(state)..aassoc(state)=e=sum(i,x(i,state)*((M_site(i)*log(XA(i,state))-M_site(i)*XA(i,state)/2)+M_site(i)/2));
aassoc_pure_define(i,state)..aassoc_pure(i,state)=e=M_site(i)*log(XA_pure(i,state) )-M_site(i)*XA_pure(i,state) /2+M_site(i)/2;

ahs_define(state)..ahs(state)=e=1/z('0',state)*(3*z('1',state)*z('2',state)/(1-z('3',state))+z('2',state)**3/z('3',state)/(1-z('3',state))**2+
                  (z('2',state)**3/z('3',state)**2-z('0',state))*log(1-z('3',state)));
ahs_pure_define(i,state)..ahs_pure(i,state)=e=1/z_pure('0',i,state)*(3*z_pure('1',i,state)*z_pure('2',i,state)/(1-z_pure('3',i,state))+z_pure('2',i,state)**3/z_pure('3',i,state)/(1-z_pure('3',i,state))**2+
                  (z_pure('2',i,state)**3/z_pure('3',i,state)**2-z_pure('0',i,state))*log(1-z_pure('3',i,state)));

ahc_define(state)..ahc(state)=e=mm(state)*ahs(state)-sum(i,x(i,state)*(m(i)-1)*log(ghs(i,i,state)));
ahc_pure_define(i,state)..ahc_pure(i,state)=e=m(i)*ahs_pure(i,state)-(m(i)-1)*log(ghs_pure(i,state));

order1_define(state)..order1(state)=e=sum((i,j),x(i,state)*x(j,state)*m(i)*m(j)*epsij(i,j,state)/T(state)*sigij(i,j)**3);
order1_pure_define(i,state)..order1_pure(i,state)=e=m(i)*m(i)*epsilon(i)/T(state)*sigma(i)**3;

order2_define(state)..order2(state)=e=sum((i,j),x(i,state)*x(j,state)*m(i)*m(j)*epsij(i,j,state)**2/T(state)**2*sigij(i,j)**3);
order2_pure_define(i,state)..order2_pure(i,state)=e=m(i)*m(i)*epsilon(i)*epsilon(i)/T(state)/T(state)*sigma(i)**3;

a_define(es_j,state)..a(es_j,state)=e=ap(es_j,'0')+(mm(state)-1)/mm(state)*ap(es_j,'1')+(mm(state)-1)*(mm(state)-2)/mm(state)/mm(state)*ap(es_j,'2');
a_pure_define(es_j,i,state)..a_pure(es_j,i,state)=e=ap(es_j,'0')+(m(i)-1)/m(i)*ap(es_j,'1')+(m(i)-1)*(m(i)-2)/m(i)/m(i)*ap(es_j,'2');

b_define(es_j,state)..b(es_j,state)=e=bp(es_j,'0')+(mm(state)-1)/mm(state)*bp(es_j,'1')+(mm(state)-1)*(mm(state)-2)/mm(state)/mm(state)*bp(es_j,'2');
b_pure_define(es_j,i,state)..b_pure(es_j,i,state)=e=bp(es_j,'0')+(m(i)-1)/m(i)*bp(es_j,'1')+(m(i)-1)*(m(i)-2)/m(i)/m(i)*bp(es_j,'2');

I1_define(state)..I1(state)=e=sum(es_j,a(es_j,state)*power(z('3',state),es_jp(es_j)));
I1_pure_define(i,state)..I1_pure(i,state)=e=sum(es_j,a_pure(es_j,i,state)*power(z_pure('3',i,state),es_jp(es_j)));

I2_define(state)..I2(state)=e=sum(es_j,b(es_j,state)*power(z('3',state),es_jp(es_j)));
I2_pure_define(i,state)..I2_pure(i,state) =e=sum(es_j,b_pure(es_j,i,state)*power(z_pure('3',i,state),es_jp(es_j)));

C1_define(state)..C1(state)=e=1/(1+mm(state)*(8*z('3',state)-2*z('3',state)**2)/(1-z('3',state))**4+(1-mm(state))*(20*z('3',state)-27*z('3',state)**2+12*
                z('3',state)**3-2*z('3',state)**4)/(1-z('3',state))**2/(2-z('3',state))**2);
C1_pure_define(i,state)..C1_pure(i,state)=e=1/(1+m(i)*(8*z_pure('3',i,state)-2*z_pure('3',i,state)**2)/(1-z_pure('3',i,state))**4+(1-m(i))*(20*z_pure('3',i,state)-27*z_pure('3',i,state)**2+12*
                z_pure('3',i,state)**3-2*z_pure('3',i,state)**4)/(1-z_pure('3',i,state))**2/(2-z_pure('3',i,state))**2);

adisp_define(state)..adisp(state)=e=-2*pi*rho(state)*I1(state)*order1(state)-pi*rho(state)*mm(state)*C1(state)*I2(state)*order2(state);
adisp_pure_define(i,state)..adisp_pure(i,state)=e=-2*pi*rho_pure(i,state)*I1_pure(i,state)*order1_pure(i,state)-pi*rho_pure(i,state)*m(i)*C1_pure(i,state)*I2_pure(i,state) *order2_pure(i,state);

ares_define(state)..ares(state)=e=adisp(state)+ahc(state)+aassoc(state);
ares_pure_define(i,state)..ares_pure(i,state)=e=adisp_pure(i,state)+ahc_pure(i,state)+aassoc_pure(i,state);

Zhs_define(state)..Zhs(state)=e=z('3',state)/(1-z('3',state))+3*z('1',state)*z('2',state)/z('0',state)/(1-z('3',state))**2+(3*z('2',state)**3-
                  z('3',state)*z('2',state)**3)/z('0',state)/(1-z('3',state))**3;
Zhs_pure_define(i,state)..Zhs_pure(i,state)=e=z_pure('3',i,state)/(1-z_pure('3',i,state))+3*z_pure('1',i,state)*z_pure('2',i,state)/z_pure('0',i,state)/(1-z_pure('3',i,state))**2+(3*z_pure('2',i,state)**3-
                  z_pure('3',i,state)*z_pure('2',i,state)**3)/z_pure('0',i,state)/(1-z_pure('3',i,state))**3;

rhogij_define(i,j,state)..rhogij(i,j,state)=e=z('3',state)/(1-z('3',state))**2+dij(i,j,state)*(3*z('2',state)/(1-z('3',state))**2+6*
                  z('2',state)*z('3',state)/(1-z('3',state))**3)+dij(i,j,state)*dij(i,j,state)*(4*z('2',state)**2/(1-z('3',state))**3+6*
                  z('2',state)**2*z('3',state)/(1-z('3',state))**4);
rhog_pure_define(i,state)..rhog_pure(i,state)=e=z_pure('3',i,state)/(1-z_pure('3',i,state))**2+d(i,state)/2*(3*z_pure('2',i,state)/(1-z_pure('3',i,state))**2+6*
                  z_pure('2',i,state)*z_pure('3',i,state)/(1-z_pure('3',i,state))**3)+d(i,state)*d(i,state)/4*(4*z_pure('2',i,state)**2/(1-z_pure('3',i,state))**3+6*
                  z_pure('2',i,state)**2*z_pure('3',i,state)/(1-z_pure('3',i,state))**4);

Zhc_define(state)..Zhc(state)=e=mm(state)*Zhs(state)-sum(i,x(i,state)*(m(i)-1)/ghs(i,i,state)*rhogij(i,i,state));
Zhc_pure_define(i,state)..Zhc_pure(i,state)=e=m(i)*Zhs_pure(i,state)-(m(i)-1)/ghs_pure(i,state)*rhog_pure(i,state);

I1z_define(state)..I1z(state)=e=sum(es_j,a(es_j,state)*(es_jp(es_j)+1)*power(z('3',state),es_jp(es_j)));
I1z_pure_define(i,state)..I1z_pure(i,state)=e=sum(es_j,a_pure(es_j,i,state)*(es_jp(es_j)+1)*power(z_pure('3',i,state),es_jp(es_j)));

I2z_define(state)..I2z(state)=e=sum(es_j,b(es_j,state)*(es_jp(es_j)+1)*power(z('3',state),es_jp(es_j)));
I2z_pure_define(i,state)..I2z_pure(i,state)=e=sum(es_j,b_pure(es_j,i,state)*(es_jp(es_j)+1)*power(z_pure('3',i,state),es_jp(es_j)));

C2_define(state)..C2(state)=e=-C1(state)*C1(state)*(mm(state)*(-4*z('3',state)**2+20*z('3',state)+8)/(1-z('3',state))**5+(1-mm(state))*(2*z('3',state)**3
                +12*z('3',state)**2-48*z('3',state)+40)/(1-z('3',state))**3/(2-z('3',state))**3);
C2_pure_define(i,state)..C2_pure(i,state)=e=-C1_pure(i,state)*C1_pure(i,state)*(m(i)*(-4*z_pure('3',i,state)**2+20*z_pure('3',i,state)+8)/(1-z_pure('3',i,state))**5+(1-m(i))*(2*z_pure('3',i,state)**3
                +12*z_pure('3',i,state)**2-48*z_pure('3',i,state)+40)/(1-z_pure('3',i,state))**3/(2-z_pure('3',i,state))**3);

Zdisp_define(state)..Zdisp(state)=e=-2*pi*rho(state)*I1z(state)*order1(state)-pi*rho(state)*mm(state)*(C1(state)*I2z(state)
                                  +C2(state)*z('3',state)*I2(state))*order2(state);
Zdisp_pure_define(i,state)..Zdisp_pure(i,state)=e=-2*pi*rho_pure(i,state)*I1z_pure(i,state)*order1_pure(i,state)-pi*rho_pure(i,state)*m(i)*(C1_pure(i,state)*I2z_pure(i,state)
                                  +C2_pure(i,state)*z_pure('3',i,state)*I2_pure(i,state) )*order2_pure(i,state);

Ztot_define(state)..Ztot(state)=e=1+Zhc(state)+Zdisp(state)+Zassoc(state);
Ztot_pure_define(i,state)..Ztot_pure(i,state)=e=1+Zhc_pure(i,state)+Zdisp_pure(i,state)+Zassoc_pure(i,state);

PVNRT(state)..0=e=P(state)-Ztot(state)*kb*T(state)*rho(state);
PVNRT_pure(i,state)..0=e=P(state)-Ztot_pure(i,state)*kb*T(state)*rho_pure(i,state);

dzdx_define(n,i,state)..dzdx(n,i,state)=e=pi/6*rho(state)*m(i)*power(d(i,state),np(n));

dahsdx_define(i,state)..dahsdx(i,state)=e=-dzdx('0',i,state)/z('0',state)*ahs(state)+1/z('0',state)*(3*(dzdx('1',i,state)*z('2',state)+z('1',state)*dzdx('2',i,state))/(1-z('3',state))
                        +3*z('1',state)*z('2',state)*dzdx('3',i,state)/(1-z('3',state))**2+3*z('2',state)**2*dzdx('2',i,state)/z('3',state)/(1-z('3',state))**2+
                        z('2',state)**3*dzdx('3',i,state)*(3*z('3',state)-1)/z('3',state)**2/(1-z('3',state))**3+((3*z('2',state)**2*z('3',state)*dzdx('2',i,state)-
                        2*z('2',state)**3*dzdx('3',i,state))/z('3',state)**3-dzdx('0',i,state))*log(1-z('3',state))+(z('0',state)-z('2',state)**3/z('3',state)**2)*
                        dzdx('3',i,state)/(1-z('3',state)));

dghsdx_define(i,j,k,state)..dghsdx(i,j,k,state)=e=dzdx('3',k,state)/(1-z('3',state))**2+dij(i,j,state)*(3*dzdx('2',k,state)/(1-z('3',state))**2+6*z('2',state)*dzdx('3',k,state)/
                        (1-z('3',state))**3)+dij(i,j,state)*dij(i,j,state)*(4*z('2',state)*dzdx('2',k,state)/(1-z('3',state))**3+6*z('2',state)**2*dzdx('3',k,state)/(1-z('3',state))**4);

dassocdx_define(i,j,k,state)..dassocdx(i,j,k,state)=e=(d(i,state)+d(j,state))**3/8*dghsdx(i,j,k,state)*kappaij(i,j)*(exp(bi_epsij(i,j)/T(state))-1);

rhoassoc_define(i,j,state)..rhoassoc(i,j,state)=e=(d(i,state)+d(j,state))**3/8*rhogij(i,j,state)*kappaij(i,j)*(exp(bi_epsij(i,j)/T(state))-1);
rhoassoc_pure_define(i,state)..rhoassoc_pure(i,state)=e=d(i,state)**3*rhog_pure(i,state)*kappa(i)*(exp(bi_eps(i)/T(state))-1);

dXAdx_define(i,k,state)..dXAdx(i,k,state)=e=-(XA(i,state)*XA(i,state))*rho(state)*(M_site(k)/2*XA(k,state)*assoc(i,k,state)+sum(j,x(j,state)*M_site(j)/2*dXAdx(j,k,state)*assoc(i,j,state)
                     +x(j,state)*M_site(j)/2*XA(j,state)*dassocdx(i,j,k,state)));

rhoXA_define(i,state)..rhoXA(i,state)=e=-(XA(i,state)*XA(i,state))*rho(state)*(sum(j,x(j,state)*M_site(j)/2*XA(j,state)*assoc(i,j,state)+x(j,state)*M_site(j)/2*(rhoXA(j,state)*assoc(i,j,state)+XA(j,state)*rhoassoc(i,j,state))));
rhoXA_pure_define(i,state)..rhoXA_pure(i,state) =e=-(XA_pure(i,state) *XA_pure(i,state) )*rho_pure(i,state)*(M_site(i)/2*XA_pure(i,state) *assoc_pure(i,state)+M_site(i)/2*(rhoXA_pure(i,state) *assoc_pure(i,state)+XA_pure(i,state) *rhoassoc_pure(i,state)));

daassocdx_define(k,state)..daassocdx(k,state)=e=M_site(k)*log(XA(k,state))-M_site(k)/2*XA(k,state)+M_site(k)/2+sum(i,x(i,state)*M_site(i)*dXAdx(i,k,state)*(1/XA(i,state)-0.5));

Zassoc_define(state)..Zassoc(state)=e=sum(i,x(i,state)*M_site(i)*rhoXA(i,state)*(1/XA(i,state)-0.5));
Zassoc_pure_define(i,state)..Zassoc_pure(i,state)=e=M_site(i)*rhoXA_pure(i,state) *(1/XA_pure(i,state) -0.5);

dahcdx_define(i,state)..dahcdx(i,state)=e=m(i)*ahs(state)+mm(state)*dahsdx(i,state)-sum(j,x(j,state)*(m(j)-1)/ghs(j,j,state)*dghsdx(j,j,i,state))-(m(i)-1)*log(ghs(i,i,state));

ord1dx_define(i,state)..ord1dx(i,state)=e=2*m(i)*sum(j,x(j,state)*m(j)*epsij(i,j,state)/T(state)*sigij(i,j)**3);

ord2dx_define(i,state)..ord2dx(i,state)=e=2*m(i)*sum(j,x(j,state)*m(j)*epsij(i,j,state)*epsij(i,j,state)/T(state)/T(state)*sigij(i,j)**3);

adx_define(es_j,i,state)..adx(es_j,i,state)=e=m(i)/mm(state)**2*ap(es_j,'1')+m(i)/mm(state)**2*(3-4/mm(state))*ap(es_j,'2');

bdx_define(es_j,i,state)..bdx(es_j,i,state)=e=m(i)/mm(state)**2*bp(es_j,'1')+m(i)/mm(state)**2*(3-4/mm(state))*bp(es_j,'2');

I1dx_define(i,state)..I1dx(i,state)=e=sum(es_j,(a(es_j,state)*es_jp(es_j)*dzdx('3',i,state)+adx(es_j,i,state)*z('3',state))*power(z('3',state),es_jp(es_j)-1));

I2dx_define(i,state)..I2dx(i,state)=e=sum(es_j,(b(es_j,state)*es_jp(es_j)*dzdx('3',i,state)+bdx(es_j,i,state)*z('3',state))*power(z('3',state),es_jp(es_j)-1));

C1dx_define(i,state)..C1dx(i,state)=e=C2(state)*dzdx('3',i,state)-C1(state)*C1(state)*(m(i)*(8*z('3',state)-2*z('3',state)**2)/(1-z('3',state))**4-m(i)*(20*z('3',state)-27*z('3',state)**2+
                  12*z('3',state)**3-2*z('3',state)**4)/(1-z('3',state))**2/(2-z('3',state))**2);

dadispdx_define(i,state)..dadispdx(i,state)=e=-2*pi*rho(state)*(I1dx(i,state)*order1(state)+I1(state)*ord1dx(i,state))-pi*rho(state)*((m(i)*C1(state)*I2(state)+mm(state)*C1dx(i,state)*I2(state)+mm(state)*C1(state)*I2dx(i,state))*
                  order2(state)+mm(state)*C1(state)*I2(state)*ord2dx(i,state));

daresdx_define(i,state)..daresdx(i,state)=e=dadispdx(i,state)+dahcdx(i,state)+daassocdx(i,state);

lnphi_define(i,state)..lnphi(i,state)=e=ares(state)+(Ztot(state)-1)-log(Ztot(state))+daresdx(i,state)-sum(j,x(j,state)*daresdx(j,state));
lnphi_pure_define(i,state)..lnphi_pure(i,state)=e=ares_pure(i,state)+(Ztot_pure(i,state)-1)-log(Ztot_pure(i,state));

model misc_check /all/;
misc_check.reslim=300;

Option
nlp=conopt;
OPTION DOMLIM =99999;

$if exist matdata.gms $include matdata.gms

solve misc_check using nlp minimizing obj;
option decimals = 8;
display rho.l, rho_pure.l, TPDF.l, x.l, obj.l;

set stat /modelStat, solveStat/;
parameter returnStat(stat),infes, objval, solver_status;

returnStat('modelStat') = misc_check.modelstat;
returnStat('solveStat') = misc_check.solvestat;
infes = misc_check.suminfes;
objval = misc_check.objval;
solver_status = misc_check.solvestat;

execute_unload %matout%;

