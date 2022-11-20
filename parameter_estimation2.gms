$set matout "'para_est_out.gdx', solver_status, infes, objval, rho_pure, rho, m, sigma, epsilon, bi_eps, kappa, ln_acti_coeff_exp, lnacti_coeff"

* Associated Content:
* Nethrue Pramuditha Mendis, Jiayuan Wang, Richard Lakerveld, "A Workflow for Crystallization Process Design with Simultaneous Process Optimization and Solvent Selection based on the Perturbed-Chain Statistical Associating Fluid Theory", Chemie Ingenieur Technik, 2022

*GAMS model for the PC-SAFT pure component parameter estimation

*The PC-SAFT part of the model is based on:
* 1. J. Gross and G. Sadowski, “Perturbed-Chain SAFT:  An Equation of State Based on a Perturbation Theory for Chain Molecules,” Ind. Eng. Chem. Res., vol. 40, no. 4, pp. 1244–1260, 2001.
* 2. J. Gross and G. Sadowski, “Application of the Perturbed-Chain SAFT Equation of State to Associating Systems,” Ind. Eng. Chem. Res., vol. 41, no. 22, pp. 5510–5515, 2002.
* 3. W. G. Chapman, K. E. Gubbins, G. Jackson, and M. Radosz, “New Reference Equation of State for Associating Liquids,” Ind. Eng. Chem. Res., vol. 29, no. 8, pp. 1709–1721, 1990.
* 4, Stanley H. Huang and Maciej Radosz, “Equation of State for Small, Large, Polydisperse, and Associating Molecules: Extension to Fluid Mixtures,” Ind. Eng. Chem. Res., vol. 30, pp. 1994–2005, 1991.

*-------------------------------------------------------------------------------
*Definitions
*-------------------------------------------------------------------------------

Set
   i(*) components
   state(*) mixtures
   para /1*5/
   alias (i,j,k);

Parameters

   M_site(i) Association site on component i

   Tm
   DHf

   P pressure: Unit in [bar]
   T temperature: Unit in [K]
   xl_ini(i,state)
   xup_ini(i,state)

   ln_acti_coeff_exp(state)

   pcsaft_para_ini(i,para)
   pcsaft_para_lo(i,para)
   pcsaft_para_up(i,para)
   rho_ini
   ;

$gdxin para_est
$load i=comp
$load state=mixture
$load M_site=assoc_site
$load Tm=melting_prop1
$load DHf=melting_prop2
$load P=pressure
$load T=temperature
$load xl_ini=mf_l
$load ln_acti_coeff_exp=ln_act_coeff
$load rho_ini=density_guess
$load pcsaft_para_ini=para_l
$load pcsaft_para_lo=para_lo
$load pcsaft_para_up=para_up
$gdxin

Positive Variables
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
   d(i) temperature dependent segment diameter
   z(n,state)  Zeta
   z_pure(n)       Zeta for the pure component
   rho(state)      total number density of mixtures
   rho_pure     total number density of pure component
   dij(i,j)   intermedia variable
   XA(i,state) mole fraction of i not bounded at site A assume that XA=XB
   XA_pure
   mm(state)  mean segment number in the mixture
   order1(state)      first order perturbation term
   order1_pure     first order perturbation term for pure component
   order2(state)     second order perturbation term
   order2_pure  second order perturbation term for pure component
   sigij(i,j)       parameters for a pair of unlike segments
   epsij(i,j,state) parameters for a pair of unlike segments
   Ztot(state)      total compressibility factor
   Ztot_pure       total compressibility factor for pure component
   ghs(i,j,state)   hard-sphere radial distribution function
   ghs_pure   hard-sphere radial distribution function for pure component
   assoc(i,j,state) association strength
   assoc_pure     association strength for pure component

Variables

   aassoc(state)     Helmholtz free energy due to association
   aassoc_pure    Helmholtz free energy due to association for pure component
   ahs(state)    Helmholtz free energy due to hard-sphere
   ahs_pure        Helmholtz free energy due to hard-sphere for pure component
   ahc(state)       Helmholtz free energy due to hard-chain
   ahc_pure         Helmholtz free energy due to hard-chain for pure component
   a(es_j,state) power series coefficient 'a'
   a_pure(es_j)
   b(es_j,state) power series coefficient 'b'
   b_pure(es_j)
   I1(state)
   I1_pure
   I2(state)
   I2_pure
   C1(state)
   C1_pure
   adisp(state)       Helmholtz free energy due to dispertion
   adisp_pure     Helmholtz free energy due to dispertion for pure component
   Zhs(state)   compressibility factor due to hard-sphere
   Zhs_pure      compressibility factor due to hard-sphere for pure component
   rhogij(i,j,state)
   rhog_pure
   Zhc(state)      compressibility factor due to hard-chain
   Zhc_pure       compressibility factor due to hard-chain for pure component
   I1z(state) I1 into eta derivative
   I1z_pure
   I2z(state) I2 into eta derivative
   I2z_pure
   C2(state)
   C2_pure
   Zdisp(state)    compressibility factor due to dispersion
   Zdisp_pure     compressibility factor due to dispersion for pure component
   Zassoc(state)    comporessibility factor due to association
   Zassoc_pure     comporessibility factor due to association for pure component
   dzdx(n,i,state) derivative of zeta
   dahsdx(i,state)
   dghsdx(i,j,k,state)
   dassocdx(i,j,k,state)
   daassocdx(k,state)
   rhoassoc(i,j,state)
   rhoassoc_pure
   dXAdx(i,k,state)
   rhoXA(i,state)
   rhoXA_pure
   Zassoc(state)
   Zassoc_pure
   dahcdx(i,state)
   ord1dx(i,state)
   ord2dx(i,state)
   I1dx(i,state)
   I2dx(i,state)
   adx(es_j,i,state) derivative of power series coefficient 'a'
   bdx(es_j,i,state) derivative of power series coefficient 'a'
   C1dx(i,state)
   dadispdx(i,state)
   daresdx(i,state)
   lnphi(i,state)
   lnphi_pure
   ares(state)
   ares_pure

   binary_interact(i,j,state)
   a_binary(i,j)
   b_binary(i,j)
   lnacti_coeff(state)
;

*-------------------------------------------------------------------------------
*Variable initialization
*-------------------------------------------------------------------------------

execseed = 1e8*(frac(jnow));

*pure component parameters

m.l(i)=pcsaft_para_ini(i,'1');
m.lo(i)=pcsaft_para_lo(i,'1');
m.up(i)=pcsaft_para_up(i,'1');

sigma.l(i)=pcsaft_para_ini(i,'2');
sigma.lo(i)=pcsaft_para_lo(i,'2');
sigma.up(i)=pcsaft_para_up(i,'2');

sigma.l(i)=pcsaft_para_ini(i,'2');
sigma.lo(i)=pcsaft_para_lo(i,'2');
sigma.up(i)=pcsaft_para_up(i,'2');

epsilon.l(i)=pcsaft_para_ini(i,'3');
epsilon.lo(i)=pcsaft_para_lo(i,'3');
epsilon.up(i)=pcsaft_para_up(i,'3');

bi_eps.l(i)=pcsaft_para_ini(i,'4');
bi_eps.lo(i)=pcsaft_para_lo(i,'4');
bi_eps.up(i)=pcsaft_para_up(i,'4');

kappa.l(i)=pcsaft_para_ini(i,'5');
kappa.lo(i)=pcsaft_para_lo(i,'5');
kappa.up(i)=pcsaft_para_up(i,'5');

*binary parameters

a_binary.fx(i,j)=0;
b_binary.fx(i,j)=0;

x.fx(i,state)=xl_ini(i,state);

*PC-SAFT-related
*-------------------------------------------------------------------------------

d.l(i)=sigma.l(i)*(1-0.12*exp(-3*epsilon.l(i)/T));

rho.l(state)=6/pi*0.5/sum(i,x.l(i,state)*m.l(i)*power(d.l(i),3));
rho.lo(state)=0.003;
rho.up(state)=0.04;

rho_pure.l=6/pi*rho_ini/(m.l('1')*power(d.l('1'),3));
rho_pure.lo=0.003;
rho_pure.up=0.04;

XA.l(i,state)=0.228;
XA.lo(i,state)=0.000001;

XA_pure.l=0.117;
XA_pure.lo=0.000001;

Ztot.l(state)=P/(kb*T*rho.l(state));
Ztot_pure.l=P/(kb*T*rho_pure.l);

z.up('3',state)=0.74048;
z_pure.up('3')=0.74048;

z.l(n,state)=pi/6*rho.l(state)*sum(i,x.l(i,state)*m.l(i)*power(d.l(i),np(n)));
z_pure.l(n)=pi/6*rho_pure.l*m.l('1')*power(d.l('1'),np(n));

dij.l(i,j)=d.l(i)*d.l(j)/(d.l(i)+d.l(j));

ghs.l(i,j,state)=1/(1-z.l('3',state))+dij.l(i,j)*3*z.l('2',state)/(1-z.l('3',state))**2+dij.l(i,j)*dij.l(i,j)
                            *2*z.l('2',state)**2/(1-z.l('3',state))**3;
ghs_pure.l=1/(1-z_pure.l('3'))+d.l('1')*3/2*z_pure.l('2')/(1-z_pure.l('3'))**2+d.l('1')*d.l('1')
                            /2*z_pure.l('2')**2/(1-z_pure.l('3'))**3;

sigij.l(i,j)=0.5*(sigma.l(i)+sigma.l(j));
epsij.l(i,j,state)=(epsilon.l(i)*epsilon.l(j))**0.5;

bi_epsij.l(i,j)=bi_eps.l(i)/2+bi_eps.l(j)/2;
kappaij.l(i,j)=(kappa.l(i)*kappa.l(j))**0.5*power((2*(sigma.l(i)*sigma.l(j))**0.5/(sigma.l(i)+sigma.l(j))),3);

assoc.l(i,j,state)=(d.l(i)+d.l(j))**3/8*ghs.l(i,j,state)*kappaij.l(i,j)*(exp(bi_epsij.l(i,j)/T)-1);
assoc_pure.l=d.l('1')**3*ghs_pure.l*kappa.l('1')*(exp(bi_eps.l('1')/T)-1);

mm.l(state)=sum(i,x.l(i,state)*m.l(i));

ahs.l(state)=1/z.l('0',state)*(3*z.l('1',state)*z.l('2',state)/(1-z.l('3',state))+z.l('2',state)**3/z.l('3',state)/(1-z.l('3',state))**2+
                  (z.l('2',state)**3/z.l('3',state)**2-z.l('0',state))*log(1-z.l('3',state)));
ahs_pure.l=1/z_pure.l('0')*(3*z_pure.l('1')*z_pure.l('2')/(1-z_pure.l('3'))+z_pure.l('2')**3/z_pure.l('3')/(1-z_pure.l('3'))**2+
                  (z_pure.l('2')**3/z_pure.l('3')**2-z_pure.l('0'))*log(1-z_pure.l('3')));

ahc.l(state)=mm.l(state)*ahs.l(state)-sum(i,x.l(i,state)*(m.l(i)-1)*log(ghs.l(i,i,state)));
ahc_pure.l=m.l('1')*ahs_pure.l-(m.l('1')-1)*log(ghs_pure.l);

order1.l(state)=sum((i,j),x.l(i,state)*x.l(j,state)*m.l(i)*m.l(j)*epsij.l(i,j,state)/T*sigij.l(i,j)**3);
order1_pure.l=m.l('1')*m.l('1')*epsilon.l('1')/T*sigma.l('1')**3;

order2.l(state)=sum((i,j),x.l(i,state)*x.l(j,state)*m.l(i)*m.l(j)*epsij.l(i,j,state)**2/T**2*sigij.l(i,j)**3);
order2_pure.l=m.l('1')*m.l('1')*epsilon.l('1')*epsilon.l('1')/T/T*sigma.l('1')**3;

a.l(es_j,state)=ap(es_j,'0')+(mm.l(state)-1)/mm.l(state)*ap(es_j,'1')+(mm.l(state)-1)*(mm.l(state)-2)/mm.l(state)/mm.l(state)*ap(es_j,'2');
a_pure.l(es_j)=ap(es_j,'0')+(m.l('1')-1)/m.l('1')*ap(es_j,'1')+(m.l('1')-1)*(m.l('1')-2)/m.l('1')/m.l('1')*ap(es_j,'2');

b.l(es_j,state)=bp(es_j,'0')+(mm.l(state)-1)/mm.l(state)*bp(es_j,'1')+(mm.l(state)-1)*(mm.l(state)-2)/mm.l(state)/mm.l(state)*bp(es_j,'2');
b_pure.l(es_j)=bp(es_j,'0')+(m.l('1')-1)/m.l('1')*bp(es_j,'1')+(m.l('1')-1)*(m.l('1')-2)/m.l('1')/m.l('1')*bp(es_j,'2');

I1.l(state)=sum(es_j,a.l(es_j,state)*power(z.l('3',state),es_jp(es_j)));
I1_pure.l=sum(es_j,a_pure.l(es_j)*power(z_pure.l('3'),es_jp(es_j)));

I2.l(state)=sum(es_j,b.l(es_j,state)*power(z.l('3',state),es_jp(es_j)));
I2_pure.l=sum(es_j,b_pure.l(es_j)*power(z_pure.l('3'),es_jp(es_j)));

C1.l(state)=1/(1+mm.l(state)*(8*z.l('3',state)-2*z.l('3',state)**2)/(1-z.l('3',state))**4+(1-mm.l(state))*(20*z.l('3',state)-27*z.l('3',state)**2+12*
                z.l('3',state)**3-2*z.l('3',state)**4)/(1-z.l('3',state))**2/(2-z.l('3',state))**2);
C1_pure.l=1/(1+m.l('1')*(8*z_pure.l('3')-2*z_pure.l('3')**2)/(1-z_pure.l('3'))**4+(1-m.l('1'))*(20*z_pure.l('3')-27*z_pure.l('3')**2+12*
                z_pure.l('3')**3-2*z_pure.l('3')**4)/(1-z_pure.l('3'))**2/(2-z_pure.l('3'))**2);

adisp.l(state)=-2*pi*rho.l(state)*I1.l(state)*order1.l(state)-pi*rho.l(state)*mm.l(state)*C1.l(state)*I2.l(state)*order2.l(state);
adisp_pure.l=-2*pi*rho_pure.l*I1_pure.l*order1_pure.l-pi*rho_pure.l*m.l('1')*C1_pure.l*I2_pure.l*order2_pure.l;

Zhs.l(state)=z.l('3',state)/(1-z.l('3',state))+3*z.l('1',state)*z.l('2',state)/z.l('0',state)/(1-z.l('3',state))**2+(3*z.l('2',state)**3-
                  z.l('3',state)*z.l('2',state)**3)/z.l('0',state)/(1-z.l('3',state))**3;
Zhs_pure.l=z_pure.l('3')/(1-z_pure.l('3'))+3*z_pure.l('1')*z_pure.l('2')/z_pure.l('0')/(1-z_pure.l('3'))**2+(3*z_pure.l('2')**3-
                  z_pure.l('3')*z_pure.l('2')**3)/z_pure.l('0')/(1-z_pure.l('3'))**3;

rhogij.l(i,j,state)=z.l('3',state)/(1-z.l('3',state))**2+dij.l(i,j)*(3*z.l('2',state)/(1-z.l('3',state))**2+6*
                  z.l('2',state)*z.l('3',state)/(1-z.l('3',state))**3)+dij.l(i,j)*dij.l(i,j)*(4*z.l('2',state)**2/(1-z.l('3',state))**3+6*
                  z.l('2',state)**2*z.l('3',state)/(1-z.l('3',state))**4);
rhog_pure.l=z_pure.l('3')/(1-z_pure.l('3'))**2+d.l('1')/2*(3*z_pure.l('2')/(1-z_pure.l('3'))**2+6*
                  z_pure.l('2')*z_pure.l('3')/(1-z_pure.l('3'))**3)+d.l('1')*d.l('1')/4*(4*z_pure.l('2')**2/(1-z_pure.l('3'))**3+6*
                  z_pure.l('2')**2*z_pure.l('3')/(1-z_pure.l('3'))**4);

Zhc.l(state)=mm.l(state)*Zhs.l(state)-sum(i,x.l(i,state)*(m.l(i)-1)/ghs.l(i,i,state)*rhogij.l(i,i,state));
Zhc_pure.l=m.l('1')*Zhs_pure.l-(m.l('1')-1)/ghs_pure.l*rhog_pure.l;

I1z.l(state)=sum(es_j,a.l(es_j,state)*(es_jp(es_j)+1)*power(z.l('3',state),es_jp(es_j)));
I1z_pure.l=sum(es_j,a_pure.l(es_j)*(es_jp(es_j)+1)*power(z_pure.l('3'),es_jp(es_j)));

I2z.l(state)=sum(es_j,b.l(es_j,state)*(es_jp(es_j)+1)*power(z.l('3',state),es_jp(es_j)));
I2z_pure.l=sum(es_j,b_pure.l(es_j)*(es_jp(es_j)+1)*power(z_pure.l('3'),es_jp(es_j)));

C2.l(state)=-C1.l(state)*C1.l(state)*(mm.l(state)*(-4*z.l('3',state)**2+20*z.l('3',state)+8)/(1-z.l('3',state))**5+(1-mm.l(state))*(2*z.l('3',state)**3
                +12*z.l('3',state)**2-48*z.l('3',state)+40)/(1-z.l('3',state))**3/(2-z.l('3',state))**3);
C2_pure.l=-C1_pure.l*C1_pure.l*(m.l('1')*(-4*z_pure.l('3')**2+20*z_pure.l('3')+8)/(1-z_pure.l('3'))**5+(1-m.l('1'))*(2*z_pure.l('3')**3
                +12*z_pure.l('3')**2-48*z_pure.l('3')+40)/(1-z_pure.l('3'))**3/(2-z_pure.l('3'))**3);

Zdisp.l(state)=-2*pi*rho.l(state)*I1z.l(state)*order1.l(state)-pi*rho.l(state)*mm.l(state)*(C1.l(state)*I2z.l(state)
                                  +C2.l(state)*z.l('3',state)*I2.l(state))*order2.l(state);
Zdisp_pure.l=-2*pi*rho_pure.l*I1z_pure.l*order1_pure.l-pi*rho_pure.l*m.l('1')*(C1_pure.l*I2z_pure.l
                                  +C2_pure.l*z_pure.l('3')*I2_pure.l)*order2_pure.l;

lnacti_coeff.l(state)=ln_acti_coeff_exp(state);

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

d_define(i)

z_define(n,state)
z_pure_define(n)

dij_define(i,j)

ghs_define(i,j,state)
ghs_pure_define

binary_interact_define(i,j,state)
binary_interact_constraint1
binary_interact_constraint2

sigij_define(i,j)
epsij_define(i,j,state)
bi_epsij_define(i,j)
kappaij_define(i,j)

assoc_define(i,j,state)
assoc_pure_define

mm_define(state)

XA_define(i,state)
XA_pure_define

aassoc_define(state)
aassoc_pure_define

ahs_define(state)
ahs_pure_define

ahc_define(state)
ahc_pure_define

order1_define(state)
order1_pure_define

order2_define(state)
order2_pure_define

a_define(es_j,state)
a_pure_define(es_j)

b_define(es_j,state)
b_pure_define(es_j)

ares_define(state)
ares_pure_define

I1_define(state)
I1_pure_define

I2_define(state)
I2_pure_define

C1_define(state)
C1_pure_define

adisp_define(state)
adisp_pure_define

Zhs_define(state)
Zhs_pure_define

rhogij_define(i,j,state)
rhog_pure_define

Zhc_define(state)
Zhc_pure_define

I1z_define(state)
I1z_pure_define

I2z_define(state)
I2z_pure_define

C2_define(state)
C2_pure_define

Zdisp_define(state)
Zdisp_pure_define

Ztot_define(state)
Ztot_pure_define

PVNRT(state)
PVNRT_pure

dzdx_define(n,i,state)

dahsdx_define(i,state)

dghsdx_define(i,j,k,state)

dassocdx_define(i,j,k,state)

rhoassoc_define(i,j,state)
rhoassoc_pure_define

dXAdx_define(i,k,state)

rhoXA_define(i,state)
rhoXA_pure_define

daassocdx_define(k,state)

Zassoc_define(state)
Zassoc_pure_define

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
lnphi_pure_define

acti_coeff_define(state)
;

*-------------------------------------------------------------------------------
*Equations
*-------------------------------------------------------------------------------

mole_frac_bound(state)..sum(i,x(i,state))=e=1;

obj_function.. obj=e=sum(state, (ln_acti_coeff_exp(state)-lnacti_coeff(state))*(ln_acti_coeff_exp(state)-lnacti_coeff(state)));

* ------------------------------------------------------------------------------

d_define(i)..d(i)=e=sigma(i)*(1-0.12*exp(-3*epsilon(i)/T));

z_define(n,state)..z(n,state)=e=pi/6*rho(state)*sum(i,x(i,state)*m(i)*power(d(i),np(n)));
z_pure_define(n)..z_pure(n)=e=pi/6*rho_pure*m('1')*power(d('1'),np(n));

dij_define(i,j)..dij(i,j)=e=d(i)*d(j)/(d(i)+d(j));

ghs_define(i,j,state)..ghs(i,j,state)=e=1/(1-z('3',state))+dij(i,j)*3*z('2',state)/(1-z('3',state))**2+dij(i,j)*dij(i,j)
                            *2*z('2',state)**2/(1-z('3',state))**3;
ghs_pure_define..ghs_pure=e=1/(1-z_pure('3'))+d('1')*3/2*z_pure('2')/(1-z_pure('3'))**2+d('1')*d('1')
                            /2*z_pure('2')**2/(1-z_pure('3'))**3;

binary_interact_define(i,j,state)..binary_interact(i,j,state)=e=a_binary(i,j)+b_binary(i,j)*T;
binary_interact_constraint1..a_binary('1','2')=e=a_binary('2','1');
binary_interact_constraint2..b_binary('1','2')=e=b_binary('2','1');

sigij_define(i,j)..sigij(i,j)=e=0.5*(sigma(i)+sigma(j));
epsij_define(i,j,state)..epsij(i,j,state)=e=(1-binary_interact(i,j,state))*(epsilon(i)*epsilon(j))**0.5;

bi_epsij_define(i,j)..bi_epsij(i,j)=e=bi_eps(i)/2+bi_eps(j)/2;
kappaij_define(i,j)..kappaij(i,j)=e=(kappa(i)*kappa(j))**0.5*power((2*(sigma(i)*sigma(j))**0.5/(sigma(i)+sigma(j))),3);

assoc_define(i,j,state)..assoc(i,j,state)=e=(d(i)+d(j))**3/8*ghs(i,j,state)*kappaij(i,j)*(exp(bi_epsij(i,j)/T)-1);
assoc_pure_define..assoc_pure=e=d('1')**3*ghs_pure*kappa('1')*(exp(bi_eps('1')/T)-1);

mm_define(state)..mm(state)=e=sum(i,x(i,state)*m(i));

XA_define(i,state)..1/XA(i,state)=e=1+sum(j,x(j,state)*rho(state)*XA(j,state)*M_site(j)/2*assoc(i,j,state));
XA_pure_define..1=e=XA_pure +rho_pure*XA_pure *XA_pure *M_site('1')/2*assoc_pure;

aassoc_define(state)..aassoc(state)=e=sum(i,x(i,state)*((M_site(i)*log(XA(i,state))-M_site(i)*XA(i,state)/2)+M_site(i)/2));
aassoc_pure_define..aassoc_pure=e=M_site('1')*log(XA_pure )-M_site('1')*XA_pure /2+M_site('1')/2;

ahs_define(state)..ahs(state)=e=1/z('0',state)*(3*z('1',state)*z('2',state)/(1-z('3',state))+z('2',state)**3/z('3',state)/(1-z('3',state))**2+
                  (z('2',state)**3/z('3',state)**2-z('0',state))*log(1-z('3',state)));
ahs_pure_define..ahs_pure=e=1/z_pure('0')*(3*z_pure('1')*z_pure('2')/(1-z_pure('3'))+z_pure('2')**3/z_pure('3')/(1-z_pure('3'))**2+
                  (z_pure('2')**3/z_pure('3')**2-z_pure('0'))*log(1-z_pure('3')));

ahc_define(state)..ahc(state)=e=mm(state)*ahs(state)-sum(i,x(i,state)*(m(i)-1)*log(ghs(i,i,state)));
ahc_pure_define..ahc_pure=e=m('1')*ahs_pure-(m('1')-1)*log(ghs_pure);

order1_define(state)..order1(state)=e=sum((i,j),x(i,state)*x(j,state)*m(i)*m(j)*epsij(i,j,state)/T*sigij(i,j)**3);
order1_pure_define..order1_pure=e=m('1')*m('1')*epsilon('1')/T*sigma('1')**3;

order2_define(state)..order2(state)=e=sum((i,j),x(i,state)*x(j,state)*m(i)*m(j)*epsij(i,j,state)**2/T**2*sigij(i,j)**3);
order2_pure_define..order2_pure=e=m('1')*m('1')*epsilon('1')*epsilon('1')/T/T*sigma('1')**3;

a_define(es_j,state)..a(es_j,state)=e=ap(es_j,'0')+(mm(state)-1)/mm(state)*ap(es_j,'1')+(mm(state)-1)*(mm(state)-2)/mm(state)/mm(state)*ap(es_j,'2');
a_pure_define(es_j)..a_pure(es_j)=e=ap(es_j,'0')+(m('1')-1)/m('1')*ap(es_j,'1')+(m('1')-1)*(m('1')-2)/m('1')/m('1')*ap(es_j,'2');

b_define(es_j,state)..b(es_j,state)=e=bp(es_j,'0')+(mm(state)-1)/mm(state)*bp(es_j,'1')+(mm(state)-1)*(mm(state)-2)/mm(state)/mm(state)*bp(es_j,'2');
b_pure_define(es_j)..b_pure(es_j)=e=bp(es_j,'0')+(m('1')-1)/m('1')*bp(es_j,'1')+(m('1')-1)*(m('1')-2)/m('1')/m('1')*bp(es_j,'2');

I1_define(state)..I1(state)=e=sum(es_j,a(es_j,state)*power(z('3',state),es_jp(es_j)));
I1_pure_define..I1_pure=e=sum(es_j,a_pure(es_j)*power(z_pure('3'),es_jp(es_j)));

I2_define(state)..I2(state)=e=sum(es_j,b(es_j,state)*power(z('3',state),es_jp(es_j)));
I2_pure_define..I2_pure =e=sum(es_j,b_pure(es_j)*power(z_pure('3'),es_jp(es_j)));

C1_define(state)..C1(state)=e=1/(1+mm(state)*(8*z('3',state)-2*z('3',state)**2)/(1-z('3',state))**4+(1-mm(state))*(20*z('3',state)-27*z('3',state)**2+12*
                z('3',state)**3-2*z('3',state)**4)/(1-z('3',state))**2/(2-z('3',state))**2);
C1_pure_define..C1_pure=e=1/(1+m('1')*(8*z_pure('3')-2*z_pure('3')**2)/(1-z_pure('3'))**4+(1-m('1'))*(20*z_pure('3')-27*z_pure('3')**2+12*
                z_pure('3')**3-2*z_pure('3')**4)/(1-z_pure('3'))**2/(2-z_pure('3'))**2);

adisp_define(state)..adisp(state)=e=-2*pi*rho(state)*I1(state)*order1(state)-pi*rho(state)*mm(state)*C1(state)*I2(state)*order2(state);
adisp_pure_define..adisp_pure=e=-2*pi*rho_pure*I1_pure*order1_pure-pi*rho_pure*m('1')*C1_pure*I2_pure *order2_pure;

ares_define(state)..ares(state)=e=adisp(state)+ahc(state)+aassoc(state);
ares_pure_define..ares_pure=e=adisp_pure+ahc_pure+aassoc_pure;

Zhs_define(state)..Zhs(state)=e=z('3',state)/(1-z('3',state))+3*z('1',state)*z('2',state)/z('0',state)/(1-z('3',state))**2+(3*z('2',state)**3-
                  z('3',state)*z('2',state)**3)/z('0',state)/(1-z('3',state))**3;
Zhs_pure_define..Zhs_pure=e=z_pure('3')/(1-z_pure('3'))+3*z_pure('1')*z_pure('2')/z_pure('0')/(1-z_pure('3'))**2+(3*z_pure('2')**3-
                  z_pure('3')*z_pure('2')**3)/z_pure('0')/(1-z_pure('3'))**3;

rhogij_define(i,j,state)..rhogij(i,j,state)=e=z('3',state)/(1-z('3',state))**2+dij(i,j)*(3*z('2',state)/(1-z('3',state))**2+6*
                  z('2',state)*z('3',state)/(1-z('3',state))**3)+dij(i,j)*dij(i,j)*(4*z('2',state)**2/(1-z('3',state))**3+6*
                  z('2',state)**2*z('3',state)/(1-z('3',state))**4);
rhog_pure_define..rhog_pure=e=z_pure('3')/(1-z_pure('3'))**2+d('1')/2*(3*z_pure('2')/(1-z_pure('3'))**2+6*
                  z_pure('2')*z_pure('3')/(1-z_pure('3'))**3)+d('1')*d('1')/4*(4*z_pure('2')**2/(1-z_pure('3'))**3+6*
                  z_pure('2')**2*z_pure('3')/(1-z_pure('3'))**4);

Zhc_define(state)..Zhc(state)=e=mm(state)*Zhs(state)-sum(i,x(i,state)*(m(i)-1)/ghs(i,i,state)*rhogij(i,i,state));
Zhc_pure_define..Zhc_pure=e=m('1')*Zhs_pure-(m('1')-1)/ghs_pure*rhog_pure;

I1z_define(state)..I1z(state)=e=sum(es_j,a(es_j,state)*(es_jp(es_j)+1)*power(z('3',state),es_jp(es_j)));
I1z_pure_define..I1z_pure=e=sum(es_j,a_pure(es_j)*(es_jp(es_j)+1)*power(z_pure('3'),es_jp(es_j)));

I2z_define(state)..I2z(state)=e=sum(es_j,b(es_j,state)*(es_jp(es_j)+1)*power(z('3',state),es_jp(es_j)));
I2z_pure_define..I2z_pure=e=sum(es_j,b_pure(es_j)*(es_jp(es_j)+1)*power(z_pure('3'),es_jp(es_j)));

C2_define(state)..C2(state)=e=-C1(state)*C1(state)*(mm(state)*(-4*z('3',state)**2+20*z('3',state)+8)/(1-z('3',state))**5+(1-mm(state))*(2*z('3',state)**3
                +12*z('3',state)**2-48*z('3',state)+40)/(1-z('3',state))**3/(2-z('3',state))**3);
C2_pure_define..C2_pure=e=-C1_pure*C1_pure*(m('1')*(-4*z_pure('3')**2+20*z_pure('3')+8)/(1-z_pure('3'))**5+(1-m('1'))*(2*z_pure('3')**3
                +12*z_pure('3')**2-48*z_pure('3')+40)/(1-z_pure('3'))**3/(2-z_pure('3'))**3);

Zdisp_define(state)..Zdisp(state)=e=-2*pi*rho(state)*I1z(state)*order1(state)-pi*rho(state)*mm(state)*(C1(state)*I2z(state)
                                  +C2(state)*z('3',state)*I2(state))*order2(state);
Zdisp_pure_define..Zdisp_pure=e=-2*pi*rho_pure*I1z_pure*order1_pure-pi*rho_pure*m('1')*(C1_pure*I2z_pure
                                  +C2_pure*z_pure('3')*I2_pure )*order2_pure;

Ztot_define(state)..Ztot(state)=e=1+Zhc(state)+Zdisp(state)+Zassoc(state);
Ztot_pure_define..Ztot_pure=e=1+Zhc_pure+Zdisp_pure+Zassoc_pure;

PVNRT(state)..0=e=P-Ztot(state)*kb*T*rho(state);
PVNRT_pure..0=e=P-Ztot_pure*kb*T*rho_pure;

dzdx_define(n,i,state)..dzdx(n,i,state)=e=pi/6*rho(state)*m(i)*power(d(i),np(n));

dahsdx_define(i,state)..dahsdx(i,state)=e=-dzdx('0',i,state)/z('0',state)*ahs(state)+1/z('0',state)*(3*(dzdx('1',i,state)*z('2',state)+z('1',state)*dzdx('2',i,state))/(1-z('3',state))
                        +3*z('1',state)*z('2',state)*dzdx('3',i,state)/(1-z('3',state))**2+3*z('2',state)**2*dzdx('2',i,state)/z('3',state)/(1-z('3',state))**2+
                        z('2',state)**3*dzdx('3',i,state)*(3*z('3',state)-1)/z('3',state)**2/(1-z('3',state))**3+((3*z('2',state)**2*z('3',state)*dzdx('2',i,state)-
                        2*z('2',state)**3*dzdx('3',i,state))/z('3',state)**3-dzdx('0',i,state))*log(1-z('3',state))+(z('0',state)-z('2',state)**3/z('3',state)**2)*
                        dzdx('3',i,state)/(1-z('3',state)));

dghsdx_define(i,j,k,state)..dghsdx(i,j,k,state)=e=dzdx('3',k,state)/(1-z('3',state))**2+dij(i,j)*(3*dzdx('2',k,state)/(1-z('3',state))**2+6*z('2',state)*dzdx('3',k,state)/
                        (1-z('3',state))**3)+dij(i,j)*dij(i,j)*(4*z('2',state)*dzdx('2',k,state)/(1-z('3',state))**3+6*z('2',state)**2*dzdx('3',k,state)/(1-z('3',state))**4);

dassocdx_define(i,j,k,state)..dassocdx(i,j,k,state)=e=(d(i)+d(j))**3/8*dghsdx(i,j,k,state)*kappaij(i,j)*(exp(bi_epsij(i,j)/T)-1);

rhoassoc_define(i,j,state)..rhoassoc(i,j,state)=e=(d(i)+d(j))**3/8*rhogij(i,j,state)*kappaij(i,j)*(exp(bi_epsij(i,j)/T)-1);
rhoassoc_pure_define..rhoassoc_pure=e=d('1')**3*rhog_pure*kappa('1')*(exp(bi_eps('1')/T)-1);

dXAdx_define(i,k,state)..dXAdx(i,k,state)=e=-(XA(i,state)*XA(i,state))*rho(state)*(M_site(k)/2*XA(k,state)*assoc(i,k,state)+sum(j,x(j,state)*M_site(j)/2*dXAdx(j,k,state)*assoc(i,j,state)
                     +x(j,state)*M_site(j)/2*XA(j,state)*dassocdx(i,j,k,state)));

rhoXA_define(i,state)..rhoXA(i,state)=e=-(XA(i,state)*XA(i,state))*rho(state)*(sum(j,x(j,state)*M_site(j)/2*XA(j,state)*assoc(i,j,state)+x(j,state)*M_site(j)/2*(rhoXA(j,state)*assoc(i,j,state)+XA(j,state)*rhoassoc(i,j,state))));
rhoXA_pure_define..rhoXA_pure =e=-(XA_pure *XA_pure )*rho_pure*(M_site('1')/2*XA_pure *assoc_pure+M_site('1')/2*(rhoXA_pure *assoc_pure+XA_pure *rhoassoc_pure));

daassocdx_define(k,state)..daassocdx(k,state)=e=M_site(k)*log(XA(k,state))-M_site(k)/2*XA(k,state)+M_site(k)/2+sum(i,x(i,state)*M_site(i)*dXAdx(i,k,state)*(1/XA(i,state)-0.5));

Zassoc_define(state)..Zassoc(state)=e=sum(i,x(i,state)*M_site(i)*rhoXA(i,state)*(1/XA(i,state)-0.5));
Zassoc_pure_define..Zassoc_pure=e=M_site('1')*rhoXA_pure *(1/XA_pure -0.5);

dahcdx_define(i,state)..dahcdx(i,state)=e=m(i)*ahs(state)+mm(state)*dahsdx(i,state)-sum(j,x(j,state)*(m(j)-1)/ghs(j,j,state)*dghsdx(j,j,i,state))-(m(i)-1)*log(ghs(i,i,state));

ord1dx_define(i,state)..ord1dx(i,state)=e=2*m(i)*sum(j,x(j,state)*m(j)*epsij(i,j,state)/T*sigij(i,j)**3);

ord2dx_define(i,state)..ord2dx(i,state)=e=2*m(i)*sum(j,x(j,state)*m(j)*epsij(i,j,state)*epsij(i,j,state)/T/T*sigij(i,j)**3);

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
lnphi_pure_define..lnphi_pure=e=ares_pure+(Ztot_pure-1)-log(Ztot_pure);

acti_coeff_define(state)..lnacti_coeff(state)=e=lnphi('1',state)-lnphi_pure;

model para_est /all/;
para_est.reslim=600;

Option
nlp=conopt;
OPTION DOMLIM =99999;

$if exist matdata.gms $include matdata.gms

solve para_est using nlp minimizing obj;

set stat /modelStat, solveStat/;
parameter returnStat(stat),infes, objval, solver_status;

returnStat('modelStat') = para_est.modelstat;
returnStat('solveStat') = para_est.solvestat;
infes = para_est.suminfes;
objval = para_est.objval;
solver_status = para_est.solvestat;

execute_unload %matout%;

