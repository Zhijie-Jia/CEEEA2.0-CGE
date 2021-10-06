$ Title CEEEA/CGE model 2.0
$Ontext
     China Energy-Environment-Economy Analysis (CEEEA) model 2.0 is based on the
     CGE model used to analyze the impact of energy policies or environmental
     policies. CEEEA/CGE 2.0 model is a recursive dynamic computable general
     equilibrium model with multi-sectors and multi- households. It can simulate
     energy supply and demand, energy-related CO2 emissions as well as the
     fundamental economic changes.

Zhijie Jia
Xi'an Jiaotong University
China

     Any questions please contact: zjjia_cn@163.com
$Offtext

*-------------------------------------------------------------------------------
* 1. Definition of sets for suffix ---------------------------------------------
Set
            u              SAM entry                     /AGR,COL,COLP,O_G,REFO,REFG,OMIN,LGT,CMC,BMTL,STL,MTL_P,MFT,THP,HYP,WDP,NCP,SOP,CST,TSPT,SER,CAP,LAB,IDT,TRF,RUR,URB,GOV,INV,ROW/
            enc(u)         energy consumption ind.       /AGR,COL,COLP,O_G,REFO,REFG,OMIN,LGT,CMC,BMTL,STL,MTL_P,MFT,THP,HYP,WDP,NCP,SOP,CST,TSPT,SER,RUR,URB/
            i(enc)         goods and services            /AGR,COL,COLP,O_G,REFO,REFG,OMIN,LGT,CMC,BMTL,STL,MTL_P,MFT,THP,HYP,WDP,NCP,SOP,CST,TSPT,SER/
            eni(i)         energy input                  /COL,COLP,O_G,REFO,REFG,THP,HYP,WDP,NCP,SOP/
            freni(i)       fossil energy                 /COL,COLP,O_G,REFO,REFG/
            neni(i)        non-energy input              /AGR,OMIN,LGT,CMC,BMTL,STL,MTL_P,MFT,CST,TSPT,SER/
            h(u)           factor input                  /CAP, LAB/
            ei(i)          env. policy-related ind.      /COLP,REFO,REFG,CMC,BMTL,STL,MTL_P,THP,HYP,WDP,NCP,SOP,CST,TSPT/
            nei(i)         non-env. policy-related ind.  /AGR,COL,O_G,OMIN,LGT,MFT,SER/
            prene(eni)     primary energy                /COL,O_G,HYP,WDP,NCP,SOP/
            prfsl(prene)   primary fossil energy         /COL,O_G/
            ele(eni)       electricity input             /THP,HYP,WDP,NCP,SOP/
            prrenew(prene) renewables                    /HYP,WDP,NCP,SOP/
            l(enc)         households                    /RUR,URB/
            ep(enc)        energy processing sector      /COLP,REFO,REFG,THP/
            nep(enc)       non-energy processing ind.    /AGR,COL,O_G,OMIN,LGT,CMC,BMTL,STL,MTL_P,MFT,HYP,WDP,NCP,SOP,CST,TSPT,SER,RUR,URB/

            add            addtional information         /Tariff_rate,depreciation_rate,capital_stock,COL,COLP,O_G,REFO,REFG,THP,HYP,WDP,NCP,SOP,rCAP_INV/
;
Alias (u,v), (i,j,jj), (h,k), (l,al);

*-------------------------------------------------------------------------------
* 2. Data processing -----------------------------------------------------------
* 2.1. Loading data from SAM (excel file) --------------------------------------
execute "gdxxrw 2018SAM(fnlbalanced).xlsx output=2018SAM(fnlbalanced).gdx par=SAM rng=SAM!A1:AE31 cdim=1 rdim=1 par=additional rng=FirmData!A1:X15 cdim=1 rdim=1 par=CO2_factor rng=CO2_factor!A1:E2 cdim=1 rdim=0 par=specificTFP rng=specificTFP!A1:U2 cdim=1 rdim=0  par=frisch rng=frisch!A1:B2 cdim=1 rdim=0  par=LESelas rng=LESelas!A1:U2 cdim=1 rdim=0 ";
Parameter       SAM(u,v)              Social accounting matrix
                additional(add,enc)   Additional firm information
                LESelas(i)            Elasticity in LES utility function
                SpecificTFP(i)        TFP grwoth parameter
                frisch(l)             Frisch parameter
                CO2_factor(freni)     CO2 emission factor from different sources
;
$gdxin   2018SAM(fnlbalanced).gdx
$loaddc  SAM
$loaddc  additional
$loaddc  LESelas
$loaddc  SpecificTFP
$loaddc  frisch
$loaddc  CO2_factor
$gdxin
;

*-------------------------------------------------------------------------------
* 2.2. Loading the initial values ----------------------------------------------
Parameter
        VAE0(j)                Energy-value added composite input
        KE0(j)                 Capital-energy composite input
        ENE0(j)                Energy input
        FOSSIL0(j)             Fossil energy input
        NOS0(j)                Non-fossil energy input
        ELC0(j)                Electricity input
        Renewable0(j)          Renewable energy input
        SOLID0(j)              Solid fossil energy input
        REF0(j)                Refined energy input

        ADEEI0(i)              Additional energy efficiency improvement
        rInvestment(i)         Capital investment rate

        F0(h,j)                Factor input h by firm j
        FS0(h,j)               Factor input (Special)
        FG0(h,j)               Factor input (General)
        X0(i,j)                Intermediate input (with energy input in SAM)
        TX0(j)                 Total intermediate input (without energy input)
        Z0(j)                  Domestic output
        Xp0(i,l)               Household consumption
        Xg0(i)                 Government consumption
        HOHincome0(l)          Household income
        Xv0(i)                 Investment demand
        E0(i)                  Exports
        M0(i)                  Imports
        Q0(i)                  Armington's composite good
        D0(i)                  Domestic consumption from domestic output
        Sp0(l)                 Private saving
        Sg0                    Government saving
        Td0(l)                 Direct tax
        Tz0(j)                 Indirect tax
        Tm0(j)                 Import tariff

        depr(i)                Depreciation rate
        CAPSTK(i)              Capital Stock
        CAPSTKS(i)             Capital
        CAPSTKG
        TOTCAP                 Total Capitial stock

        RFF0(l,h)              Rate of factor input
        FF0(l,h)               Factor endowment of the h-th factor
        FFs0(h)                Special capital endowment
        FFg0(h)                General capital endowment
        r_Sepcial              Rate of special capital
        Sf                     Foreign saving in international price
        pWe(i)                 Export price in international price
        pWm(i)                 Import price in international price
        tauz(i)                Production tax rate
        taum(i)                Import tariff rate

        EM0(enc)               CO2 emissions
        Energy0(enc,enc)       Energy consumption
        coe(eni,enc)           Coefficient between value and physical quantity of energy consumption
        TOT_Energy0(eni)       Total energy consumption in energy source eni
        Emissions0             Total CO2 emissions

        Indentity_Matrix(i,j)  Indentity Matrix
        Cd0(j)                 CO2 emissions per unit final demand
        A0(i,j)                Direct consumption coefficient
        Am0(i,j)               Direct requirement coefficient matrix of the intermediate input from imports
        XLeontiefinverse0(i,j) Leontief Reverse of Direct consumption coefficient Matrix
        Ed0(j)                 Domestic embodied emissions per unit final demand
        FD0(j)                 Final demand
        Eim_intermediate0(j)   Intermediate process of Eim
        Eim0(j)                Emissions of imported intermediate input per unit final demand
        Yim0(j)                Imported directed domestic final consumption
        EEP0                   Total embodied emissions from domestic production
        EEC0                   Total embodied emissions from domestic consumption
        EEE0                   Emissions embodied within exports
        EEI0                   Emissions embodied within imports
        EEB0                   Net embodied emissions of trade balance

        PPI0(j)                Producer price index
        CPI0                   Consumer price index

        Renewable_rate0        Rate of renewable energy in total primary energy
        r_electricity0         Rate of electricity in total energy input
        TOTelectricity0        Total electricity consumption
        GDP0                   Gross Domestic Product

        BAUDummy               Business as Usual Dummy(1 if policy is implemented otherwise 0)
        CFDummy                Counter factual Dummy(1 if policy is implemented otherwise 0)

        CO2_factor             CO2 emission factor of standard coal equivalent
        eff_coal_coke          energy transfer efficiency of coke producer
        eff_coal_thp           efficiency of Grid transmission
        eff_O_G_oil            energy transfer efficiency of product oil producer
        eff_O_G_gas            energy transfer efficiency of product gas producer
;
BAUDummy       =  0;
CFDummy        =  0;

eff_coal_coke  =  1- 0.032593257;
eff_coal_thp   =  1- 0.046771409;
eff_O_G_oil    =  1- 0.04550507;
eff_O_G_gas    =  1- 0;
r_Sepcial      =  0.5;
pWe(i)         =  1;
pWm(i)         =  1;

Td0(l)         =  SAM("GOV",l);
Tz0(j)         =  SAM("IDT",j);
Tm0(j)         =  SAM("TRF",J);

F0(h,j)        =  SAM(h,j);
Fs0(h,j)       =  r_Sepcial*SAM(h,j);
Fg0(h,j)       =  (1-r_Sepcial)*SAM(h,j);
X0(i,j)        =  SAM(i,j);
TX0(j)         =  sum(neni,X0(neni,j));
KE0(j)         =  sum(eni,X0(eni,j))+F0("CAP",j);
ENE0(j)        =  sum(eni,X0(eni,j));
Renewable0(j)  =  X0("HyP",j)+X0("WdP",j)+X0("NcP",j)+X0("SoP",j);
ELC0(j)        =  X0("ThP",j)+X0("HyP",j)+X0("WdP",j)+X0("NcP",j)+X0("SoP",j);
REF0(j)        =  X0("REFO",j)+X0("REFG",j);
NOS0(j)        =  REF0(j)+X0("O_G",j);
SOLID0(j)      =  X0("COL",j)+X0("COLP",j);
FOSSIL0(j)     =  SOLID0(j)+NOS0(j);
VAE0(j)        =  KE0(j)+F0("LAB",j);

Z0(j)          =  VAE0(j) +sum(neni, X0(neni,j));
M0(i)          =  SAM("ROW",i);

tauz(j)        =  Tz0(j)/(Z0(j)+Tz0(j));
taum(j)$(M0(j))=  Tm0(j)/M0(j);

Xp0(i,l)       =  SAM(i,l);
FF0(l,h)       =  SAM(l,h);
FFs0(h)        =  r_Sepcial*sum(l,FF0(l,h));
FFg0(h)        =  (1-r_Sepcial)*sum(l,FF0(l,h));
RFF0(l,h)      =  FF0(l,h)/sum(al,FF0(al,h));
HOHincome0(l)  =  sum((h,j),F0(h,j)*RFF0(l,h));
Xg0(i)         =  SAM(i,"GOV");
Xv0(i)         =  SAM(i,"INV");
E0(i)          =  SAM(i,"ROW");
Q0(i)          =  sum(l,Xp0(i,l))+Xg0(i)+Xv0(i)+sum(j, X0(i,j));
D0(i)          =  Z0(i)/(1-tauz(i))-E0(i);
Sp0(l)         =  SAM("INV",l);
Sg0            =  SAM("INV","GOV");
Sf             =  SAM("INV","ROW");

depr(i)        =  additional("depreciation_rate",i);
CAPSTK(i)      =  F0("CAP",i)/depr(i);
CAPSTKS(i)     =  r_Sepcial*CAPSTK(i);
CAPSTKG        =  (1-r_Sepcial)*sum(i,CAPSTK(i));
TOTCAP         =  sum(i,CAPSTK(i));
rInvestment(i) =  additional("rCAP_INV",i);

Energy0("COL",enc)  =  additional("COL",enc);
Energy0("COLP",enc) =  additional("COLP",enc);
Energy0("O_G",enc)  =  additional("O_G",enc);
Energy0("REFO",enc) =  additional("REFO",enc);
Energy0("REFG",enc) =  additional("REFG",enc);
Energy0("THP",enc)  =  additional("THP",enc);
Energy0("HYP",enc)  =  additional("HYP",enc);
Energy0("WDP",enc)  =  additional("WDP",enc);
Energy0("NCP",enc)  =  additional("NCP",enc);
Energy0("SOP",enc)  =  additional("SOP",enc);
TOT_Energy0(eni)    =  sum(enc,Energy0(eni,enc));
coe(eni,j)$(X0(eni,j))  =  Energy0(eni,j)/X0(eni,j);
coe(eni,l)$(Xp0(eni,l)) =  Energy0(eni,l)/Xp0(eni,l);
EM0(nep)            =  sum(freni,Energy0(freni,nep)*CO2_factor(freni))     +Energy0("COL","THP") *Energy0("THP",nep)   /sum(enc,Energy0("THP",enc))*CO2_factor("COL")  ;
EM0("COLP")         =  sum(freni,Energy0(freni,"COLP")*CO2_factor(freni))  +Energy0("COL","THP") *Energy0("THP","COLP")/sum(enc,Energy0("THP",enc))*CO2_factor("COL") -Energy0("COL","COLP")*eff_coal_coke*CO2_factor("COL");
EM0("THP")          =  sum(freni,Energy0(freni,"THP")*CO2_factor(freni))   +Energy0("COL","THP") *Energy0("THP","THP") /sum(enc,Energy0("THP",enc))*CO2_factor("COL") -Energy0("COL","THP") *eff_coal_thp*CO2_factor("COL");
EM0("REFO")         =  sum(freni,Energy0(freni,"REFO")*CO2_factor(freni))  +Energy0("COL","THP") *Energy0("THP","REFO")/sum(enc,Energy0("THP",enc))*CO2_factor("COL") -Energy0("O_G","REFO")*eff_O_G_oil*CO2_factor("O_G");
EM0("REFG")         =  sum(freni,Energy0(freni,"REFG")*CO2_factor(freni))  +Energy0("COL","THP") *Energy0("THP","REFG")/sum(enc,Energy0("THP",enc))*CO2_factor("COL") -Energy0("O_G","REFG")*eff_O_G_gas*CO2_factor("O_G");

Indentity_Matrix(i,j)   =  0;
Indentity_Matrix(i,i)   =  1;

Cd0(j)              =  sum(freni,Energy0(freni,j)*CO2_factor(freni))/( sum(i,X0(i,j)) +sum(h,F0(h,j)) +Tz0(j) +Tm0(j) );
A0(i,j)             =  X0(i,j)/( sum(jj,X0(jj,j)) +sum(h,F0(h,j)) +Tz0(j) +Tm0(j) );
Am0(i,j)            =  M0(i)/Q0(i)*A0(i,j);

Parameter A0_inter(i,j);
A0_inter(i,j)=Indentity_Matrix(i,j)-A0(i,j);
execute_unload 'gdxforinverse.gdx' i,A0_inter;
execute 'invert gdxforinverse.gdx i A0_inter gdxfrominverse.gdx XLeontiefinverse0';
execute_load 'gdxfrominverse.gdx' , XLeontiefinverse0;

Ed0(j)              =  sum(jj,Cd0(jj)*XLeontiefinverse0(jj,j)) ;
FD0(j)              =  Xg0(j)+Xv0(j)+sum(l,Xp0(j,l))+E0(j)-M0(j);
Eim_intermediate0(j)=  sum(jj, Ed0(jj)*Am0(jj,j));
Eim0(j)             =  sum(jj, Eim_intermediate0(jj) *XLeontiefinverse0(jj,j));
Yim0(j)             =  M0(j)/Q0(j)*FD0(j);

EEP0                =  sum(j,Ed0(j)*FD0(j));
EEC0                =  sum(j,Ed0(j)*(FD0(j)-E0(j))) +sum(j,Eim0(j)*(FD0(j)-E0(j))) +sum(j,Ed0(j)*Yim0(j));
EEE0                =  sum(j,Ed0(j)*E0(j)) +sum(j,Eim0(j)*E0(j));
EEI0                =  sum(j,Ed0(j)*Yim0(j)) +sum(j,Eim0(j)*FD0(j));
EEB0                =  EEP0-EEC0;

PPI0(j)             =  100;
CPI0                =  100;

r_electricity0      =  sum((ele,enc), Energy0(ele,enc))/sum((eni,enc), Energy0(eni,enc));
TOTelectricity0     =  sum((ele,enc),energy0(ele,enc));
Emissions0          =  sum(enc, EM0(enc));
Renewable_rate0     =  sum((prrenew,enc),Energy0(prrenew,enc))/sum((prene,enc),Energy0(prene,enc));
ADEEI0(i)           =  0;
GDP0                =  sum(i,Xv0(i)+Xg0(i)+sum(l,Xp0(i,l))+E0(i)-M0(i));

*-------------------------------------------------------------------------------
* 2.3. Calibration -------------------------------------------------------------
* 2.3.1. Setting elasticity-----------------------------------------------------
Parameter
         sigma(i)              elasticity of substitution in internaional trade
         psi(i)                elasticity of transformation
         eta(i)                substitution elasticity parameter
         phi(i)                transformation elasticity parameter

         rhovae(j)             elasticity of substitution in energy-value added composition
         rhof(j)               elasticity of substitution in factor input
         rhoke(j)              elasticity of substitution in capital-energy input
         rhoene(j)             elasticity of substitution in energy input
         rhofossil(enc)        elasticity of substitution in fossil energy input
         rhosolid(enc)         elasticity of substitution in solid fossil energy input
         rhonos(enc)           elasticity of substitution in non-solid fossil energy input
         rhoref(enc)           elasticity of substitution in refined fossil energy input
         rhoelc(j)             elasticity of substitution in electricity input
         rhorenewable(j)       elasticity of substitution in renewable power input
         rhoz(j)               elasticity of substitution in domestic production function

         sigma_eff(j)          ADEEI elasticity

         Sensitivity           for the use of elasticity sensitivity analysis (default is 1)
;
sigma(i)     =  2;
psi(i)       =  2;
eta(i)       =  (sigma(i)-1)/sigma(i);
phi(i)       =  (psi(i)+1)/psi(i);

sensitivity  =  1.0;

rhovae(j)    =  1-1/1.5*sensitivity;
rhof(j)      =  1-1/0.4*sensitivity;
rhoke(j)     =  1-1/1.3*sensitivity;
rhoene(j)    =  1-1/2.0*sensitivity;
rhofossil(j) =  1-1/1.5*sensitivity;
rhonos(j)    =  1-1/1.5*sensitivity;
rhosolid(j)  =  1-1/1.5*sensitivity;
rhoref(j)    =  1-1/1.5*sensitivity;
rhoelc(j)    =  1-1/10*sensitivity;
rhorenewable(j)  =  1-1/10*sensitivity;
rhoz(j)      =  1-1/1.2*sensitivity;

sigma_eff(j) =  0.1;

parameter        elasticity_EP  Elasticity of substitution in energy processing sectors
;
elasticity_EP    =  0.3;
rhofossil(ep)    =  1-1/elasticity_EP;
rhosolid(ep)     =  1-1/elasticity_EP;
rhonos(ep)       =  1-1/elasticity_EP;
rhoref(ep)       =  1-1/elasticity_EP;

*-------------------------------------------------------------------------------
* 2.3.2. Calibration of the exogenous variables---------------------------------
Parameter
                    deltavae(j)          Share parameter in energy-value added input
                    alphavae(j)          Scale parameter in energy-value added input
                    deltake(j)           Share parameter in capital-energy input
                    alphake(j)           Scale parameter in capital-energy input
                    deltaf(j)            Share parameter in cpital input
                    alphaf(j)            Scale parameter in cpital input
                    deltaene(j)          Share parameter in energy input
                    alphaene(j)          Scale parameter in energy input
                    deltafossil(j)       Share parameter in fossil energy input
                    alphafossil(j)       Scale parameter in fossil energy input
                    deltasolid(j)        Share parameter in solid fossil energy input
                    alphasolid(j)        Scale parameter in solid fossil energy input
                    deltanos(j)          Share parameter in non-solid fossil energy input
                    alphanos(j)          Scale parameter in non-solid fossil energy input
                    deltaref(j)          Share parameter in refined fossil energy input
                    alpharef(j)          Scale parameter in refined fossil energy input
                    deltaelc(j)          Share parameter in electricity input
                    alphaelc(j)          Scale parameter in electricity input
                    deltahyp(j)          Share parameter in hydropower
                    deltawdp(j)          Scale parameter in wind power
                    deltancp(j)          Share parameter in nuclear power
                    deltasop(j)          Share parameter in solar power
                    alpharenewable(j)    Scale parameter in renewable energy
                    deltaz(j)            Share parameter in domestic production function
                    alphaz(j)            Scale parameter in domestic production function

                    ax(i,j)              Intermediate input requirement coefficient
                    mu(i)                Government consumption share
                    lambda(i)            Onvestment demand share
                    deltam(i)            Share parameter in Armington function
                    deltad(i)            Share parameter in Armington function
                    gamma(i)             Scale parameter in Armington function
                    xid(i)               Share parameter in transformation function
                    xie(i)               Share parameter in transformation function
                    theta(i)             Scale parameter in transformation function
                    ssp(l)               Average propensity for private saving
                    ssg                  Average propensity for government saving
                    taud(l)              Direct tax rate

                    bgtshr(i,l)          Budget share in LES utility function
                    LESbeta(i,l)         Marginal consumption share in LES utility function
                    LESsub0(j,l)         Survival consumption of consumption function
                    UU0(l)               Utility in household l
                    SW0                  Social welfare

                    penebau(j)           energy price in BAU scenario
                    pzbau(j)             producer price in BAU scenario
                    pqbau(j)             Armington price in BAU scenario
                    Xpbau(j,l)           Household consumption in BAU scenario
;

deltaf(j)         =  Fg0("CAP",j)**(1-rhof(j))/(Fg0("CAP",j)**(1-rhof(j))+Fs0("CAP",j)**(1-rhof(j)));
alphaf(j)         =  F0("CAP",j)/((deltaf(j)*Fg0("CAP",j)**rhof(j)+(1-deltaf(j))*Fs0("CAP",j)**rhof(j))**(1/rhof(j)));

deltavae(j)       =  F0("LAB",j)**(1-rhovae(j))/(F0("LAB",j)**(1-rhovae(j))+KE0(j)**(1-rhovae(j)));
alphavae(j)       =  VAE0(j)/(deltavae(j)*F0("LAB",j)**rhovae(j)+(1-deltavae(j))*KE0(j)**rhovae(j))**(1/rhovae(j));
deltake(j)        =  F0("CAP",j)**(1-rhoke(j))/(ENE0(j)**(1-rhoke(j))+F0("CAP",j)**(1-rhoke(j)));
alphake(j)        =  KE0(j)/(deltake(j)*F0("CAP",j)**rhoke(j)+(1-deltake(j))*ENE0(j)**rhoke(j))**(1/rhoke(j));

deltaene(j)       =  ELC0(j)**(1-rhoene(j))/(ELC0(j)**(1-rhoene(j))+FOSSIL0(j)**(1-rhoene(j)));
alphaene(j)       =  ENE0(j)/(deltaene(j)*ELC0(j)**rhoene(j)+(1-deltaene(j))*FOSSIL0(j)**rhoene(j))**(1/rhoene(j));

deltafossil(j)    =  SOLID0(j)**(1-rhofossil(j))/(SOLID0(j)**(1-rhofossil(j))+NOS0(j)**(1-rhofossil(j)));
alphafossil(j)    =  FOSSIL0(j)/(deltafossil(j)*SOLID0(j)**rhofossil(j)+(1-deltafossil(j))*NOS0(j)**rhofossil(j))**(1/rhofossil(j));

deltasolid(j)$(X0("COLP",j))  =  X0("COL",j)**(1-rhosolid(j))/(X0("COL",j)**(1-rhosolid(j))+X0("COLP",j)**(1-rhosolid(j))) ;
alphasolid(j)$(X0("COLP",j))  =  SOLID0(j)/(deltasolid(j)*X0("COL",j)**rhosolid(j)+(1-deltasolid(j))*X0("COLP",j)**rhosolid(j))**(1/rhosolid(j));

deltanos(j)$(X0("O_G",j))     =  X0("O_G",j)**(1-rhonos(j))/(X0("O_G",j)**(1-rhonos(j))+REF0(j)**(1-rhonos(j))) ;
alphanos(j)$(X0("O_G",j))     =  NOS0(j)/(deltanos(j)*X0("O_G",j)**rhonos(j)+(1-deltanos(j))*REF0(j)**rhonos(j))**(1/rhonos(j));

deltaref(j)       =  X0("REFO",j)**(1-rhoref(j))/(X0("REFO",j)**(1-rhoref(j))+X0("REFG",j)**(1-rhoref(j)));
alpharef(j)       =  REF0(j)/(deltaref(j)*X0("REFO",j)**rhoref(j)+(1-deltaref(j))*X0("REFG",j)**rhoref(j))**(1/rhoref(j));

deltaelc(j)       =  X0("ThP",j)**(1-rhoelc(j))/(X0("ThP",j)**(1-rhoelc(j))+Renewable0(j)**(1-rhoelc(j)));
alphaelc(j)       =  ELC0(j)/(deltaelc(j)*X0("ThP",j)**rhoelc(j)+(1-deltaelc(j))*Renewable0(j)**rhoelc(j))**(1/rhoelc(j));

deltahyp(j)       =  X0("HyP",j)**(1-rhorenewable(j))/(X0("HyP",j)**(1-rhorenewable(j))+X0("WdP",j)**(1-rhorenewable(j))+X0("NcP",j)**(1-rhorenewable(j))+X0("SoP",j)**(1-rhorenewable(j)) );
deltawdp(j)       =  X0("WdP",j)**(1-rhorenewable(j))/(X0("HyP",j)**(1-rhorenewable(j))+X0("WdP",j)**(1-rhorenewable(j))+X0("NcP",j)**(1-rhorenewable(j))+X0("SoP",j)**(1-rhorenewable(j)) );
deltancp(j)       =  X0("NcP",j)**(1-rhorenewable(j))/(X0("HyP",j)**(1-rhorenewable(j))+X0("WdP",j)**(1-rhorenewable(j))+X0("NcP",j)**(1-rhorenewable(j))+X0("SoP",j)**(1-rhorenewable(j)) );
deltasop(j)       =  X0("SoP",j)**(1-rhorenewable(j))/(X0("HyP",j)**(1-rhorenewable(j))+X0("WdP",j)**(1-rhorenewable(j))+X0("NcP",j)**(1-rhorenewable(j))+X0("SoP",j)**(1-rhorenewable(j)) );
alpharenewable(j) =  Renewable0(j)/(deltahyp(j)*X0("HyP",j)**rhorenewable(j)+deltawdp(j)*X0("WdP",j)**rhorenewable(j)+deltancp(j)*X0("NcP",j)**rhorenewable(j)+deltasop(j)*X0("SoP",j)**rhorenewable(j))**(1/rhorenewable(j));

deltaz(j)         =  VAE0(j)**(1-rhoz(j))/(VAE0(j)**(1-rhoz(j))+TX0(j)**(1-rhoz(j)));
alphaz(j)         =  Z0(j)/(deltaz(j)*VAE0(j)**rhoz(j)+(1-deltaz(j))*TX0(j)**rhoz(j))**(1/rhoz(j)) ;

ax(neni,j)        =  X0(neni,j)/TX0(j);
mu(i)             =  Xg0(i)/sum(j, Xg0(j));
lambda(i)         =  Xv0(i)/(sum(l,Sp0(l))+Sg0+Sf);

deltam(i)$(M0(i)) =  (1+taum(i))*M0(i)**(1-eta(i)) / ((1+taum(i))*M0(i)**(1-eta(i)) +D0(i)**(1-eta(i)));
deltad(i)$(M0(i)) =  D0(i)**(1-eta(i)) / ((1+taum(i))*M0(i)**(1-eta(i)) +D0(i)**(1-eta(i)));
gamma(i)$(M0(i))  =  Q0(i)/(deltam(i)*M0(i)**eta(i)+deltad(i)*D0(i)**eta(i))**(1/eta(i));

xie(i)$(E0(i))    =  E0(i)**(1-phi(i))/(E0(i)**(1-phi(i))+D0(i)**(1-phi(i)));
xid(i)$(E0(i))    =  D0(i)**(1-phi(i))/(E0(i)**(1-phi(i))+D0(i)**(1-phi(i)));
theta(i)$(E0(i))  =  Z0(i) /(xie(i)*E0(i)**phi(i)+xid(i)*D0(i)**phi(i))**(1/phi(i));

ssp(l)            =  Sp0(l)/sum(h, FF0(l,h));
ssg               =  Sg0/(sum(l, Td0(l)) +sum(j, Tz0(j)) +sum(j, Tm0(j)));
taud(l)           =  Td0(l)/sum(h, FF0(l,h));

bgtshr(i,l)       =  Xp0(i,l)/sum(j, Xp0(j,l));
LESbeta(i,l)      =  LESelas(i)*bgtshr(i,l)/sum(j,LESelas(j)*bgtshr(j,l));
LESsub0(i,l)      =  Xp0(i,l)+LESbeta(i,l)/1*((sum((h,j),F0(h,j)*RFF0(l,h))-Sp0(l)-Td0(l))/Frisch(l));
UU0(l)            =  prod(i, (Xp0(i,l) -LESsub0(i,l))**LESbeta(i,l));
SW0               =  SUM(l,UU0(l));

penebau(j)        =  1;
pzbau(j)          =  1;
pqbau(j)          =  1;
Xpbau(j,l)        =  Xp0(j,l);

*-------------------------------------------------------------------------------
* 3. Model declaration----------------------------------------------------------
* 3.1. Endogenous variables-----------------------------------------------------
Variable
                    VAE(j)                Energy-value added composite input
                    F(h,j)                the factor input h by firm j
                    ENE(j)                Energy input
                    KE(j)                 Capital-energy composite input
                    ENE(j)                Energy input
                    FOSSIL(j)             Fossil energy input
                    SOLID(j)              Solid fossil energy input
                    NOS(j)                Non-solid fossil energy input
                    REF(j)                Refined energy input
                    ELC(j)                Electricity input
                    Renewable(j)          Renewable energy input
                    X(i,j)                Intermediate input
                    TX(j)                 Total intermediate input
                    Z(j)                  output of the j-th good

                    ADEEI(i)              Additional energy efficiency improvement

                    FF(l,h)               Factor endowment
                    RFF(l,h)              rate of FF of rural population and citizen
                    Fs(h,j)               Factor input (Special)
                    Fg(h,j)               Factor input (General)
                    FFg(h)                General capital input endowment
                    HOHincome(l)          Household income
                    Xp(i,l)               Household consumption
                    Xg(i)                 Government consumption
                    Xv(i)                 Investment demand
                    E(i)                  Exports
                    M(i)                  Imports
                    Q(i)                  Armington's composite good
                    D(i)                  Domestic consumption from domestic output

                    pke(j)                Price of Capital-energy composite input
                    pene(j)               Price of Energy input
                    pfossil(j)            Price of Fossil energy input
                    psolid(j)             Price of Solid fossil energy input
                    pnos(j)               Price of Non-solid fossil energy input
                    pelc(j)               Price of Electricity input
                    prenewable(j)         Price of Renewable energy input
                    pref(j)               Price of Refined energy input
                    px(eni)               Price of energy input
                    ptx(j)                Price of total intermediate input
                    pf(h,j)               Price of factor input
                    pfg                   Price of general capital input
                    pfs(j)                Price of special capital input
                    pvae(j)               Price of Energy-value added composite input
                    pz(j)                 Price of Domestic output
                    pq(i)                 Price of Armington's composite good
                    pe(i)                 Price of export goods
                    pm(i)                 Price of import goods
                    pd(i)                 Price of domestic consumption from domestic output
                    epsilon               Exchange rate

                    Energy(enc,enc)       Energy consumption (including final consumption and transformation)
                    TOT_Energy(eni)       Total energy consumption
                    EM(enc)               CO2 emissions
                    PLC(i)                Policy cost of enterprises
                    ctr                   Carbon tax rate

                    Sp(l)                 Private saving
                    Sg                    Government saving
                    Td(l)                 Direct tax
                    Tz(j)                 Indirect tax
                    Tm(i)                 Import tariff

                    Cd(j)                 CO2 emissions per unit final demand
                    A(i,j)                Direct consumption coefficient
                    Am(i,j)               Direct requirement coefficient matrix of the intermediate input from imports
                    XLeontiefinverse(i,j) Leontief Reverse of Direct consumption coefficient Matrix
                    Ed(j)                 Domestic embodied emissions per unit final demand
                    FD(j)                 Final demand
                    Eim_intermediate(j)   Intermediate process of Eim
                    Eim(j)                Emissions of imported intermediate input per unit final demand
                    Yim(j)                Imported directed domestic final consumption

                    EEP                   Total embodied emissions from domestic production
                    EEC                   Total embodied emissions from domestic consumption
                    EEE                   Emissions embodied within exports
                    EEI                   Emissions embodied within imports
                    EEB                   Net embodied emissions of trade balance

                    PPI(j)                Producer price index (BAU=100)
                    CPI                   Consumer price index (BAU=100)

                    Renewable_Rate        Rate of renewable energy in total primary energy
                    Emissions             Total CO2 emissions
                    r_electricity         Rate of electricity in total energy input
                    TOTelectricity        Total electricity consumption
                    EV(l)                 Equivalent variation
                    CV(l)                 Compensate variation
                    GDP                   Gross Domestic Product
                    GDPCHK                Model checking variable using different types of GDP accounting (equal or near to 0 if pass)
                    UU(l)                 Utility of household l
                    SW                    Social welfare

                    C_AEEI                Change of Automatic Energy Efficiency Improvement
                    C_TFP                 Change of Total Factor Productivity
                    C_ELE                 Change of Electricity efficiency

                    Walras                Walras Dummy
                    Walras2               Walras Dummy 2

;

* ------------------------------------------------------------------------------
* 3.2. Model Equation-----------------------------------------------------------
Equation
                    eqCESy1(j)            CES production function of energy-value added
                    eqCESy2(j)            First order condition in CES production function of energy-value added
                    eqCESy3(j)            Value balance function of energy-value added
                    eqCESke1(j)           CES production function of capital-energy
                    eqCESke2(j)           First order condition in CES production function of capital-energy
                    eqCESke3(j)           Value balance function of capital-energy
                    eqCESene1(j)          CES production function of energy
                    eqCESene2(j)          First order condition in CES production function of energy
                    eqCESene3(j)          Value balance function of energy
                    eqCESfossil1(j)       CES production function of fossil energy
                    eqCESfossil2(j)       First order condition in CES production function of fossil energy
                    eqCESfossil3(j)       Value balance function of fossil energy
                    eqCESsolid1(j)        CES production function of solid fossil energy
                    eqCESsolid2(j)        First order condition in CES production function of solid fossil energy
                    eqCESsolid3(j)        Value balance function of solid fossil energy
                    eqCESsolid_S1(j)      CES production function of solid fossil energy (for those ind. lacks of coke input)
                    eqCESsolid_S2(j)      First order condition in CES production function of solid fossil energy (for those ind. lacks of coke input)
                    eqCESsolid_S3(j)      Value balance function of solid fossil energy (for those ind. lacks of coke input)
                    eqCESref1(j)          CES production function of refined energy
                    eqCESref2(j)          First order condition in CES production function of refined energy
                    eqCESref3(j)          Value balance function of refined energy
                    eqCESnos1(j)          CES production function of non-solid fossil energy
                    eqCESnos2(j)          First order condition in CES production function of non-solid fossil energy
                    eqCESnos3(j)          Value balance function of non-solid fossil energy
                    eqCESnos_S1(j)        CES production function of non-solid fossil energy (for those ind. lacks of crude oil and gas input)
                    eqCESnos_S2(j)        First order condition in CES production function of non-solid fossil energy (for those ind. lacks of crude oil and gas input)
                    eqCESnos_S3(j)        Value balance function of non-solid fossil energy (for those ind. lacks of crude oil and gas input)
                    eqCESelc1(j)          CES production function of electricity
                    eqCESelc2(j)          First order condition in CES production function of electricity
                    eqCESelc3(j)          Value balance function of  electricity
                    eqCESrenewable1(j)    CES production function of renewable power
                    eqCESrenewable2(j)    First order condition in CES production function of renewable power
                    eqCESrenewable3(j)    First order condition in CES production function of renewable power
                    eqCESrenewable4(j)    First order condition in CES production function of renewable power
                    eqCESrenewable5(j)    Value balance function of renewable power

                    eqpxpq(eni)           Balance of intermediate input price and Armington price
                    eqcoe1(eni,j)         Coefficient between value and physical quantity of energy in sectors
                    eqcoe2(eni,l)         Coefficient between value and physical quantity of energy in households

                    eqADEEI(j)            Addtional Energy Effifiency Improvment cause by energy-environment policy or energy price

                    eqX(i,j)              First-order condition of Leontief production function of intermediate input
                    eqpzs(j)              Value balance in Leontief production function

                    eqZ1(j)               CES production function of domestic output
                    eqZ2(j)               First order condition in CES production function of domestic output
                    eqZ3(j)               Value balance function of domestic output

                    eqem1(enc)            Emission equation by using energy from non-enery processing ind.
                    eqem2                 Emission equation by using energy from coke producer
                    eqem3                 Emission equation by using energy from thermal power producer
                    eqem4                 Emission equation by using energy from refined oil producer
                    eqem5                 Emission equation by using energy from refined gas producer
                    eqPLC1(i)             Policy cost by energy or environment policy
                    eqPLC2(i)             Policy cost by energy or environment policy

                    eqTd(l)               Direct tax revenue function
                    eqTz(j)               Indirect tax revenue function
                    eqTm(i)               Import tariff revenue function
                    eqXg(i)               Government demand function

                    eqXv(i)               Investment demand function
                    eqSp(l)               Private saving function
                    eqSg                  Government saving function

                    eqXp(i,l)             Household demand function
                    eqHOHincome(l)        Household income function

                    eqpe(i)               World export price equation
                    eqpm(i)               World import price equation
                    eqepsilon             Trade deficit function

                    eqpqs(i)              Armington function (CES function)
                    eqM(i)                Import demand function (first-order condition)
                    eqD(i)                Domestic good demand function (first-order condition)
                    eqpqs1(i)             Armington function (CES function)(for goods with no import)
                    eqpqs2(i)             Import demand function (first-order condition)(for goods with no import)
                    eqpqs3(i)             Domestic good demand function (first-order condition)(for goods with no import)

                    eqpzd(i)              Transformation function (CET function)
                    eqDs(i)               Domestic good supply function (first-order condition)
                    eqE(i)                Export supply function (first-order condition)
                    eqpzd1(i)             Transformation function (CET function)(for goods with no export)
                    eqpzd2(i)             Domestic good supply function (first-order condition)(for goods with no export)
                    eqpzd3(i)             Export supply function (first-order condition)(for goods with no export)

                    eqpqd(i)              Market clearing condition for Armington goods
                    eqpf(h)               Factor market clearing condition
                    eqpfLAB               Labor price balance function
                    eqffg                 General capital endowment
                    eqffs(j)              Sepcial capital endowment
                    eqFcap1(j)            CES function of Capital input
                    eqFcap2(j)            First-order condition in CES function of Capital input
                    eqFcap3(j)            Value balance function

                    eqCd(j)               CO2 emissions per unit final demand
                    eqDrtCspCoeff(i,j)    Direct consumption coefficient
                    eqDrtCspCoeff_im(i,j) Direct consumption coefficient matrix of the intermediate input from imports
                    eqLeontiefInverse(i,j) Leontief Reverse of Direct consumption coefficient Matrix
                    eqEd(j)               Domestic embodied emissions per unit final demand
                    eqFD(j)               Final demand
                    eqEim_intermediate(j) Intermediate process of Eim
                    eqEim(j)              Emissions of imported intermediate input per unit final demand
                    eqYim(j)              Imported directed domestic final consumption

                    eqEEP                 Total embodied emissions from domestic production
                    eqEEC                 Total embodied emissions from domestic comsuption
                    eqEEE                 Emissions embodied within exports
                    eqEEI                 Emissions embodied within imports
                    eqEEB                 Net embodied emissions of trade balance

                    eqPPI(j)              Producer price index
                    eqCPI                 Consumer price index

                    eqEmissions           Total CO2 emission
                    eqTOT_energy(eni)     Total energy consumption for energy source eni
                    eqr_electricity       Rate of electricity in total energy input
                    eqTOTelectricity      Total electricity consumption
                    eqrenewable           Rate of renewable energy in total primary energy
                    eqEV(l)               Equivalent Variation of utility
                    eqCV(l)               Compensate Variation of utility
                    eqGDP                 GDP function (calculated by expenditure method)
                    eqGDPCHK              Model check using different GDP calculation
                    eqUU(l)               Utility function
                    eqSW                  Social welfare function
;
*3.2.1. Production block -------------------------------------------------------
eqCESy1(j)..                            VAE(j)  =e=  alphavae(j)*(1+C_TFP*SpecificTFP(j))*(deltavae(j)*F("LAB",j)**rhovae(j)+(1-deltavae(j))*KE(j)**rhovae(j))**(1/rhovae(j));
eqCESy2(j)..                pf("LAB",j)/PKE(J)  =e=  deltavae(j)/(1-deltavae(j))*(KE(j)/F("LAB",j))**(1-rhovae(j));
eqCESy3(j)..                    pvae(j)*VAE(j)  =e=  pf("LAB",j)*F("LAB",j)+PKE(j)*KE(j);

eqCESke1(j)..                            KE(j)  =e=  alphake(j)*(deltake(j)*F("CAP",j)**rhoke(j)+(1-deltake(j))*(ENE(j)/(1-ADEEI(j)))**rhoke(j))**(1/rhoke(j));
eqCESke2(j)..              pf("CAP",j)/PENE(j)  =e=  deltake(j)/(1-deltake(j))*((ENE(j)/(1-ADEEI(j)))/F("CAP",j))**(1-rhoke(j));
eqCESke3(j)..                     PKE(j)*KE(j)  =e=  PENE(j)*ENE(j)+ pf("CAP",j)*F("CAP",j);

eqADEEI(j)..                           ADEEI(j) =e=  (  1-( pene(j)/penebau(j) )**(-sigma_eff(j))  )*(1-BAUDummy);

eqFcap1(j)..                        F("CAP",j)  =e=  alphaf(j)*(deltaf(j)*Fg("CAP",j)**rhof(j)+(1-deltaf(j))*Fs("CAP",j)**rhof(j))**(1/rhof(j));
eqFcap2(j)..                        pfg/pfs(j)  =e=  deltaf(j)/(1-deltaf(j))*(Fs("CAP",j)/Fg("CAP",j))**(1-rhof(j));
eqFcap3(j)..            F("CAP",j)*pf("CAP",j)  =e=  pfg*Fg("CAP",j)+pfs(j)*FS("CAP",j);

eqCESene1(j)..                          ENE(j)  =e=  alphaene(j)*(deltaene(j)*ELC(j)**rhoene(j)+(1-deltaene(j))*FOSSIL(j)**rhoene(j))**(1/rhoene(j));
eqCESene2(j)..              PELC(j)/PFOSSIL(j)  =e=  deltaene(j)/(1-deltaene(j))*(FOSSIL(j)/ELC(j))**(1-rhoene(j));
eqCESene3(j)..                  PENE(J)*ENE(j)  =e=  PELC(j)*ELC(j)+PFOSSIL(j)*FOSSIL(j);

eqCESfossil1(j)..                    FOSSIL(j)  =e=  alphafossil(j)*(1+C_AEEI)*(deltafossil(j)*SOLID(j)**rhofossil(j)+(1-deltafossil(j))*NOS(j)**rhofossil(j))**(1/rhofossil(j));
eqCESfossil2(j)..            PSOLID(j)/PNOS(j)  =e=  deltafossil(j)/(1-deltafossil(j))*(NOS(j)/SOLID(j))**(1-rhofossil(j));
eqCESfossil3(j)..         PFOSSIL(j)*FOSSIL(j)  =e=  PSOLID(j)*SOLID(j)+PNOS(j)*NOS(j)+PLC(j);

eqCESsolid1(j)$(X0("COLP",j)<>0)..                 SOLID(j)  =e=  alphasolid(j)*(deltasolid(j)*X("COL",j)**rhosolid(j)+(1-deltasolid(j))*X("COLP",j)**rhosolid(j))**(1/rhosolid(j));
eqCESsolid2(j)$(X0("COLP",j)<>0)..     PX("COL")/PX("COLP")  =e=  deltasolid(j)/(1-deltasolid(j))*(X("COLP",j)/X("COL",j))**(1-rhosolid(j));
eqCESsolid3(j)$(X0("COLP",j)<>0)..       PSOLID(j)*SOLID(j)  =e=  PX("COL")*X("COL",j)+PX("COLP")*X("COLP",j);
eqCESsolid_S1(j)$(X0("COLP",j)=0)..                SOLID(j)  =e=  X("COL",j);
eqCESsolid_S2(j)$(X0("COLP",j)=0)..               PSOLID(j)  =e=  PX("COL");
eqCESsolid_S3(j)$(X0("COLP",j)=0)..                       0  =e=  X("COLP",j);

eqCESnos1(j)$(X0("O_G",j)<>0)..                      NOS(j)  =e=  alphanos(j)*(deltanos(j)*X("O_G",j)**rhonos(j)+(1-deltanos(j))*REF(j)**rhonos(j))**(1/rhonos(j));
eqCESnos2(j)$(X0("O_G",j)<>0)..           PX("O_G")/PREF(j)  =e=  deltanos(j)/(1-deltanos(j))*(REF(j)/X("O_G",j))**(1-rhonos(j));
eqCESnos3(j)$(X0("O_G",j)<>0)..              PNOS(j)*NOS(j)  =e=  PX("O_G")*X("O_G",j)+PREF(j)*REF(j);
eqCESnos_S1(j)$(X0("O_G",j)=0)..                    PNOS(j)  =e=  PREF(j);
eqCESnos_S2(j)$(X0("O_G",j)=0)..                     NOS(j)  =e=  REF(j);
eqCESnos_S3(j)$(X0("O_G",j)=0)..                          0  =e=  X("O_G",j);

eqCESref1(j)..                          REF(j)  =e=  alpharef(j)*(deltaref(j)*X("REFO",j)**rhoref(j)+(1-deltaref(j))*X("REFG",j)**rhoref(j))**(1/rhoref(j));
eqCESref2(j)..           PX("REFO")/PX("REFG")  =e=  deltaref(j)/(1-deltaref(j))*(X("REFG",j)/X("REFO",j))**(1-rhoref(j));
eqCESref3(j)..                  PREF(j)*REF(j)  =e=  PX("REFO")*X("REFO",j)+PX("REFG")*X("REFG",j);

eqCESelc1(j)..                          ELC(j)  =e=  alphaelc(j)*(1+C_ELE)*(deltaelc(j)*X("ThP",j)**rhoelc (j)+(1-deltaelc(j))*Renewable(j)**rhoelc(j))**(1/rhoelc(j));
eqCESelc2(j)..         PX("ThP")/prenewable(j)  =e=  deltaelc(j)/(1-deltaelc(j))*(Renewable(j)/X("ThP",j))**(1-rhoelc(j));
eqCESelc3(j)..                  ELC(j)*PELC(j)  =e=  X("ThP",j)*PX("ThP")+PRenewable(j)*Renewable(j) ;

eqCESrenewable1(j)..              Renewable(j)  =e=  alpharenewable(j)*(deltahyp(j)*X("HyP",j)**rhorenewable(j)+deltawdp(j)*X("WdP",j)**rhorenewable(j)+deltancp(j)*X("NcP",j)**rhorenewable(j)+deltasop(j)*X("SoP",j)**rhorenewable(j))**(1/rhorenewable(j));
eqCESrenewable2(j)..       PX("SoP")/PX("HyP")  =e=  deltasop(j)/deltahyp(j)*(X("HyP",j)/X("SoP",j))**(1-rhorenewable(j));
eqCESrenewable3(j)..       PX("SoP")/PX("WdP")  =e=  deltasop(j)/deltawdp(j)*(X("WdP",j)/X("SoP",j))**(1-rhorenewable(j));
eqCESrenewable4(j)..       PX("SoP")/PX("NcP")  =e=  deltasop(j)/deltancp(j)*(X("NcP",j)/X("SoP",j))**(1-rhorenewable(j));
eqCESrenewable5(j).. PRenewable(j)*Renewable(j) =e= PX("HyP")*X("HyP",j) +PX("WdP")*X("WdP",j)+PX("NcP")*X("NcP",j)+PX("SoP")*X("SoP",j);

eqpxpq(eni)..                          px(eni)  =e=  pq(eni);

eqX(neni,j)..                        X(neni,j)  =e=  ax(neni,j)*TX(j);
eqpzs(j)..                              ptx(j)  =e=  sum(neni, ax(neni,j)*pq(neni));

eqZ1(j)..                                 Z(j)  =e=  alphaz(j)*(deltaz(j)*VAE(j)**rhoz(j)+(1-deltaz(j))*TX(j)**rhoz(j))**(1/rhoz(j));
eqZ2(j)..                       pvae(j)/ptx(j)  =e=  deltaz(j)/(1-deltaz(j))*(TX(j)/VAE(j))**(1-rhoz(j));
eqZ3(j)..                           pz(j)*Z(j)  =e=  pvae(j)*VAE(j)+ptx(j)*TX(j);

* 3.2.2. Income & Expenditure block --------------------------------------------
eqTd(l)..                                Td(l)  =e=  taud(l)*sum((h,j), pf(h,j)*F(h,j)*RFF(l,h));
eqTz(i)..                                Tz(i)  =e=  tauz(i)/(1-tauz(i))*pz(i)*Z(i);
eqTm(i)..                                Tm(i)  =e=  taum(i)*pm(i)*M(i);
eqXg(i)..                                Xg(i)  =e=  mu(i)*(sum(l,Td(l)) +sum(j, Tz(j)) +sum(j, Tm(j))  -Sg +sum(j,PLC(j)) )/pq(i);

eqXv(i)..                                Xv(i)  =e=  lambda(i)*(sum(l,Sp(l)) +Sg +epsilon*Sf)/pq(i);

eqSp(l)..                                Sp(l)  =e=  ssp(l)*sum((h,j), pf(h,j)*F(h,j)*RFF(l,h));
eqSg..                                      Sg  =e=  ssg*(sum(l,Td(l)) +sum(j, Tz(j)) +sum(j, Tm(j)) +sum(j,PLC(j)) );

eqXp(i,l)..                            Xp(i,l)  =e=  LESsub0(i,l)+LESbeta(i,l)/pq(i)*( sum((h,j),pf(h,j)*F(h,j)*RFF(l,h))-Sp(l)-Td(l) -sum(j,LESsub0(j,l)*pq(j)) );
eqHOHincome(l)..                  HOHincome(l)  =e=  sum((h,j),pf(h,j)*F(h,j)*RFF(l,h));

* 3.2.3. Trade block -----------------------------------------------------------
eqpe(i)..                                pe(i)  =e=  epsilon*pWe(i);
eqpm(i)..                                pm(i)  =e=  epsilon*pWm(i);
eqepsilon..            sum(i, pWe(i)*E(i)) +Sf  =e=  sum(i, pWm(i)*M(i));

eqpqs(i)$(M0(i)<>0)..                     Q(i)  =e=  gamma(i)*(deltam(i)*M(i)**eta(i)+deltad(i) *D(i)**eta(i))**(1/eta(i));
eqM(i)$(M0(i)<>0)..                       M(i)  =e=  (gamma(i)**eta(i)*deltam(i)*pq(i) /((1+taum(i))*pm(i)))**(1/(1-eta(i)))*Q(i);
eqD(i)$(M0(i)<>0)..                       D(i)  =e=  (gamma(i)**eta(i)*deltad(i)*pq(i)/pd(i))**(1/(1-eta(i)))*Q(i);
eqpqs1(i)$(M0(i)=0)..                     Q(i)  =e=  D(i);
eqpqs2(i)$(M0(i)=0)..                    pq(i)  =e=  pd(i);
eqpqs3(i)$(M0(i)=0)..                     M(i)  =e=  0;

eqpzd(i)$(E0(i)<>0)..                     Z(i)  =e=  theta(i)*(xie(i)*E(i)**phi(i)+xid(i)*D(i)**phi(i))**(1/phi(i));
eqE(i)$(E0(i)<>0)..                       E(i)  =e=  (theta(i)**phi(i)*xie(i)*1/(1-tauz(i))*pz(i)/pe(i))**(1/(1-phi(i)))*Z(i);
eqDs(i)$(E0(i)<>0)..                      D(i)  =e=  (theta(i)**phi(i)*xid(i)*1/(1-tauz(i))*pz(i)/pd(i))**(1/(1-phi(i)))*Z(i);
eqpzd1(i)$(E0(i)=0)..         Z(i)/(1-tauz(i))  =e=  D(i);
eqpzd2(i)$(E0(i)=0)..                    pz(i)  =e=  pd(i);
eqpzd3(i)$(E0(i)=0)..                     E(i)  =e=  0;

* 3.2.4. Energy-Environment block ----------------------------------------------
eqcoe1(eni,j)..            coe(eni,j)*X(eni,j)  =e=  Energy(eni,j);
eqcoe2(eni,l)..           coe(eni,l)*Xp(eni,l)  =e=  Energy(eni,l);

eqem1(nep)..                           EM(nep)  =e=  sum(freni,Energy(freni,nep)*CO2_factor(freni))     +Energy("COL","THP") *Energy("THP",nep)   /sum(enc,Energy("THP",enc))*CO2_factor("COL")  ;
eqem2("COLP")..                     EM("COLP")  =e=  sum(freni,Energy(freni,"COLP")*CO2_factor(freni))  +Energy("COL","THP") *Energy("THP","COLP")/sum(enc,Energy("THP",enc))*CO2_factor("COL") -Energy("COL","COLP")*eff_coal_coke*CO2_factor("COL");
eqem3("THP")..                       EM("THP")  =e=  sum(freni,Energy(freni,"THP")*CO2_factor(freni))   +Energy("COL","THP") *Energy("THP","THP") /sum(enc,Energy("THP",enc))*CO2_factor("COL") -Energy("COL","THP") *eff_coal_thp*CO2_factor("COL");
eqem4("REFO")..                     EM("REFO")  =e=  sum(freni,Energy(freni,"REFO")*CO2_factor(freni))  +Energy("COL","THP") *Energy("THP","REFO")/sum(enc,Energy("THP",enc))*CO2_factor("COL") -Energy("O_G","REFO")*eff_O_G_oil*CO2_factor("O_G");
eqem5("REFG")..                     EM("REFG")  =e=  sum(freni,Energy(freni,"REFG")*CO2_factor(freni))  +Energy("COL","THP") *Energy("THP","REFG")/sum(enc,Energy("THP",enc))*CO2_factor("COL") -Energy("O_G","REFG")*eff_O_G_gas*CO2_factor("O_G");

eqPLC1(ei)..                           PLC(ei)  =e=  ctr*EM(ei)*CFDummy ;
eqPLC2(nei)..                         PLC(nei)  =e=  0;

eqEmissions..                        Emissions  =e=  sum(enc, EM(enc));
eqTOT_energy(eni)..            TOT_ENERGY(eni)  =e=  sum(enc, Energy(eni,enc));

* 3.2.5. market-clearing & Macroscopic closure block ---------------------------
eqpqd(i)..                                Q(i)  =e=  sum(l, Xp(i,l)) +Xg(i) +Xv(i) +sum(j, X(i,j))+Walras;
eqpf("LAB")..               sum(j, F("LAB",j))  =e=  sum(l, FF(l,"LAB"));
eqpfLAB(j)..                       pf("LAB",j)  =e=  sum(i,pf("LAB",i))/card(j)+Walras2;
eqffg..                     sum(j,FG("CAP",j))  =e=  FFG("CAP");
eqffs(j)..                         FS("CAP",j)  =e=  CAPSTKS(j)*depr(j);

* 3.2.6 Macro-indicators function ----------------------------------------------
eqCd(j)..                                Cd(j)  =e=  sum(freni,Energy(freni,j)*CO2_factor(freni))/( sum(jj,pq(jj)*X(jj,j)) +sum(h,pf(h,j)*F(h,j)) +Tz(j) +Tm(j) +PLC(j));
eqDrtCspCoeff(i,j)..                    A(i,j)  =e=  pq(i)*X(i,j)/( sum(jj,pq(jj)*X(jj,j)) +sum(h,pf(h,j)*F(h,j)) +Tz(j) +Tm(j) +PLC(j));
eqDrtCspCoeff_im(i,j)..                Am(i,j)  =e=  pm(i)*M(i)/(pq(i)*Q(i))*A(i,j);
eqLeontiefInverse(i,j).. Indentity_Matrix(i,j)  =e=  sum(jj, (Indentity_Matrix(i,jj)-A(i,jj))*XLeontiefinverse(jj,j) );
eqEd(j)..                                Ed(j)  =e=  sum(jj,Cd(jj)*XLeontiefinverse(jj,j)) ;
eqFD(j)..                                FD(j)  =e=  pq(j) *(Xg(j)+Xv(j)+sum(l,Xp(j,l))) +pe(j)*E(j)-pm(j)*M(j);
eqEim_intermediate(j)..    Eim_intermediate(j)  =e=  sum(jj, Ed(jj)*Am(jj,j));
eqEim(j)..                              Eim(j)  =e=  sum(jj, Eim_intermediate(jj) *XLeontiefinverse(jj,j));
eqYim(j)..                              Yim(j)  =e=  pm(j)*M(j)/(pq(j)*Q(j))*FD(j);

eqEEP..                                    EEP  =e=  sum(j,Ed(j)*FD(j));
eqEEC..                                    EEC  =e=  sum(j,Ed(j)*(FD(j)-pe(j)*E(j))) +sum(j,Eim(j)*(FD(j)-pe(j)*E(j))) +sum(j,Ed(j)*Yim(j));
eqEEE..                                    EEE  =e=  sum(j,Ed(j)*pe(j)*E(j)) +sum(j,Eim(j)*pe(j)*E(j));
eqEEI..                                    EEI  =e=  sum(j,Ed(j)*Yim(j)) +sum(j,Eim(j)*FD(j));
eqEEB..                                    EEB  =e=  EEP-EEC;

eqPPI(j)..                              PPI(j)  =e=  pz(j)/pzbau(j)*100;
eqCPI..                                    CPI  =e=  sum(j, pq(j)*sum(l,Xpbau(j,l))/sum((i,l),Xpbau(i,l))) / sum(j, pqbau(j)*sum(l,Xpbau(j,l))/sum((i,l),Xpbau(i,l)))*100;

eqrenewable..                   Renewable_rate  =e=  sum((prrenew,enc),Energy(prrenew,enc))/sum((prene,enc),Energy(prene,enc));
eqr_electricity..                r_electricity  =e=  sum((ele,enc), Energy(ele,enc))/sum((eni,enc), Energy(eni,enc));
eqTOTelectricity..              TOTelectricity  =e=  sum((ele,enc),energy(ele,enc));

eqEV(l)..                                EV(l)  =e=  (UU.l(l)-UU0(l))*prod(i$(LESbeta(i,l)), (1/LESbeta(i,l))**LESbeta(i,l) ) ;
eqCV(l)..                                CV(l)  =e=  (UU.l(l)-UU0(l))*prod(i$(LESbeta(i,l)),(pq(i)/LESbeta(i,l))**LESbeta(i,l)) ;
eqGDP..                                    GDP  =e=  sum(i,Xv(i)+Xg(i)+sum(l,Xp(i,l))+E(i)-M(i));
eqGDPCHK..                              GDPCHK  =e=  sum((h,j), pf(h,j)*F(h,j)) +sum(j, Tz(j)+Tm(j) +PLC(j))  -sum(i, pq(i)*(Xv(i)+Xg(i)+sum(l,Xp(i,l))) +pe(i)*E(i) -pm(i)*M(i) ) ;

eqUU(l)..                                UU(l)  =e=  prod(i, (Xp(i,l) -LESsub0(i,l))**LESbeta(i,l));
eqSW..                                      SW  =e=  sum(l, UU(l));

* ------------------------------------------------------------------------------
* 3.3. Initializing endogenous variables ---------------------------------------
VAE.l(j)                =  VAE0(j);
F.l(h,j)                =  F0(h,j);
FS.l(h,j)               =  FS0(h,j);
FG.l(h,j)               =  FG0(h,j);
KE.l(j)                 =  KE0(j);
ENE.l(j)                =  ENE0(j);
FOSSIL.l(j)             =  FOSSIL0(j);
SOLID.l(j)              =  SOLID0(j);
NOS.l(j)                =  NOS0(j);
ELC.l(j)                =  ELC0(j);
Renewable.l(j)          =  Renewable0(j);
REF.l(j)                =  REF0(j);
X.l(i,j)                =  X0(i,j);
TX.l(j)                 =  TX0(j);
Z.l(j)                  =  Z0(j);
HOHincome.l(l)          =  HOHincome0(l);
Xp.l(i,l)               =  Xp0(i,l);
Xg.l(i)                 =  Xg0(i);
Xv.l(i)                 =  Xv0(i);
E.l(i)                  =  E0(i);
M.l(i)                  =  M0(i);
Q.l(i)                  =  Q0(i);
D.l(i)                  =  D0(i);
pf.l(h,j)               =  1;
pfg.l                   =  1;
pfs.l(j)                =  1;
pvae.l(j)               =  1;
pz.l(j)                 =  1;
pq.l(i)                 =  1;
pe.l(i)                 =  1;
pke.l(j)                =  1;
pene.l(i)               =  1;
pfossil.l(i)            =  1;
psolid.l(i)             =  1;
pnos.l(i)               =  1;
pelc.l(j)               =  1;
prenewable.l(j)         =  1;
pref.l(j)               =  1;
px.l(eni)               =  1;
ptx.l(j)                =  1;
FFg.l(h)                =  FFg0(h);
RFF.l(l,h)              =  RFF0(l,h);
pz.l(j)                 =  1;
pq.l(i)                 =  1;
pe.l(i)                 =  1;
pm.l(i)                 =  1;
pd.l(i)                 =  1;
epsilon.l               =  1;
Sp.l(l)                 =  Sp0(l);
Sg.l                    =  Sg0;
Td.l(l)                 =  Td0(l);
Tz.l(j)                 =  Tz0(j);
Tm.l(i)                 =  Tm0(i);
Energy.l(eni,enc)       =  Energy0(eni,enc);
TOT_Energy.l(eni)       =  TOT_Energy0(eni);
EM.l(enc)               =  EM0(enc);
PLC.l(i)                =  0 ;
Cd.l(j)                 =  Cd0(j);
A.l(i,j)                =  A0(i,j);
Am.l(i,j)               =  Am0(i,j);
XLeontiefinverse.l(i,j) =  XLeontiefinverse0(i,j);
Ed.l(j)                 =  Ed0(j);
FD.l(j)                 =  FD0(j);
Eim_intermediate.l(j)   =  Eim_intermediate0(j);
Eim.l(j)                =  Eim0(j);
Yim.l(j)                =  Yim0(j);
EEP.l                   =  EEP0;
EEC.l                   =  EEC0;
EEE.l                   =  EEE0;
EEI.l                   =  EEI0;
EEB.l                   =  EEB0;
PPI.l(j)                =  100;
CPI.l                   =  100;
Renewable_rate.l        =  Renewable_rate0;
Emissions.l             =  Emissions0;
r_electricity.l         =  r_electricity0;
TOTelectricity.l        =  TOTelectricity0;
EV.l(l)                 =  0;
CV.l(l)                 =  0;
GDP.l                   =  GDP0;
GDPCHK.l                =  0;
UU.l(l)                 =  UU0(l);
SW.l                    =  SW0;
ADEEI.l(i)              =  ADEEI0(i);
Walras.l                =  0;
Walras2.l               =  0;

* ------------------------------------------------------------------------------
* 3.4. Model setting -----------------------------------------------------------
* Instead of pf.fx("LAB","AGR"), You can also set other price as benchmark price, such as pq.fx("AGR").
pf.fx("LAB","AGR")      =  1;
*pq.fx("AGR")           =  1;

* FF and FFG are exogenous variables but decleared as endogenous ones, so here we fix the value.
FFG.fx("CAP")           =  FFG0("CAP");
FF.fx(l,h)              =  FF0(l,h);

* ctr is the carbon tax rate, which is endogenous in this example. Here we assume that government will set a appropriate rate to meet the emission traget that is exogenous set.
ctr.fx                  =  0;

* The dyanmic parameter simulated by DyanimcParameter.gms file based on the growth path needed.
C_TFP.fx                =  0;
C_AEEI.fx               =  0;
C_ELE.fx                =  0;

* RFF is the share of factor income between rural and urban households. It is also can be dynamic set.
RFF.fx(l,h)             =  RFF0(l,h);

* BAUDummy is a Dummy variable only exist in BAU scenario and CF scenarios before treating the policy.
BAUDummy                =  1;

* ------------------------------------------------------------------------------
* 3.5. Defining and solving the model ------------------------------------------
Model CEEEA /all/;
Solve CEEEA using mcp;

* ------------------------------------------------------------------------------
* 4. Dynamics ------------------------------------------------------------------
* 4.1. setting dynamic variables -----------------------------------------------
set                  t                           Time priod      /2018 * 2060/ ;
scalar               gr                          Labor endowment growth         /0.0061/      ;
Parameter            Energy1(t,enc,enc)
                     TOT_Energy1(t,eni)
                     EM1(t,enc)
                     PLC1(t,i)
                     ctr1(t)
                     VAE1(t,j)
                     F1(t,h,j)
                     ENE1(t,j)
                     KE1(t,j)
                     ENE1(t,j)
                     FOSSIL1(t,j)
                     SOLID1(t,j)
                     NOS1(t,j)
                     REF1(t,j)
                     ELC1(t,j)
                     X1(t,i,j)
                     TX1(t,j)
                     Z1(t,j)

                     ADEEI1(t,i)

                     HOHincome1(t,l)
                     Xp1(t,i,l)
                     Xg1(t,i)
                     Xv1(t,i)
                     E1(t,i)
                     M1(t,i)
                     Q1(t,i)
                     D1(t,i)

                     pke1(t,j)
                     pene1(t,j)
                     pfossil1(t,j)
                     psolid1(t,j)
                     pnos1(t,j)
                     pelc1(t,j)
                     pref1(t,j)
                     px1(t,eni)
                     pf1(t,h,j)
                     pvae1(t,j)
                     ptx1(t,j)
                     pz1(t,j)
                     pq1(t,i)
                     pe1(t,i)
                     pm1(t,i)
                     pd1(t,i)
                     epsilon1(t)

                     Sp1(t,l)
                     Sg1(t)
                     Td1(t,l)
                     Tz1(t,j)
                     Tm1(t,i)

                     Cd1(t,j)
                     A1(t,i,j)
                     Am1(t,i,j)
                     XLeontiefinverse1(t,i,j)
                     Ed1(t,j)
                     FD1(t,j)
                     Eim_intermediate1(t,j)
                     Eim1(t,j)
                     Yim1(t,j)
                     EEP1(t)
                     EEC1(t)
                     EEE1(t)
                     EEI1(t)
                     EEB1(t)
                     PPI1(t,j)
                     CPI1(t)
                     Renewable_rate1(t)
                     Emissions1(t)
                     r_electricity1(t)
                     TOTelectricity1(t)
                     EV1(t,l)
                     CV1(t,l)
                     GDP1(t)
                     GDPCHK1(t)
                     UU1(t,l)
                     SW1(t)

                     C_AEEI1(t)
                     C_TFP1(t)
                     C_ELE1(t)

                     Walras1(t)
                     Walras21(t)
                     FF1(t,l,h)
                     RFF1(t,l,h)

                     penebau1(t,j)
                     pzbau1(t,j)
                     pqbau1(t,j)
                     Xpbau1(t,j,l)
;

* ------------------------------------------------------------------------------
* 4.2. Extract data of growth path set in need ---------------------------------

$ontext
In BAU scenario, we only consider changes of C_TFP, C_AEEI and C_ELE, which is
simulated by DynamicParameter.gms file.
$offtext

*parameter
*                    CFemissions(t)
*                    BAUemissions(t)
*                    BAUGDP(t)
*                    BAUTOTelectricity(t)
*;
*execute "gdxxrw Growth_Path.xlsx output=Growth_Path.gdx par=CFemissions rng=CFemissions!A1:AQ2 cdim=1 rdim=0  par=BAUemissions rng=BAUemissions!A1:AQ2 cdim=1 rdim=0  par=BAUGDP rng=BAUGDP!A1:AQ2 cdim=1 rdim=0 par=BAUTOTelectricity rng=BAUTOTelectricity!A1:AQ2 cdim=1 rdim=0";
*$gdxin Growth_Path.gdx
*$loaddc CFemissions
*$loaddc BAUemissions
*$loaddc BAUGDP
*$loaddc BAUTOTelectricity
*;
$gdxin output-DynamicParameter.gdx
$loaddc C_TFP1
$loaddc C_AEEI1
*$loaddc penebau1
$loaddc C_ELE1
$loaddc pzbau1
$loaddc pqbau1
$loaddc Xpbau1

*-------------------------------------------------------------------------------
* 4.3. Running the dynamic model -----------------------------------------------
* 4.3.1 Assigning output values ------------------------------------------------
Loop(t,
Energy1(t,enc,enc)       =        Energy.l(enc,enc);
TOT_Energy1(t,eni)       =        TOT_Energy.l(eni);
EM1(t,enc)               =        EM.l(enc);
PLC1(t,i)                =        PLC.l(i);
ctr1(t)                  =        ctr.l;
VAE1(t,j)                =        VAE.l(j);
F1(t,h,j)                =        F.l(h,j);
ENE1(t,j)                =        ENE.l(j);
KE1(t,j)                 =        KE.l(j);
ENE1(t,j)                =        ENE.l(j);
FOSSIL1(t,j)             =        FOSSIL.l(j);
SOLID1(t,j)              =        SOLID.l(j);
NOS1(t,j)                =        NOS.l(j);
REF1(t,j)                =        REF.l(j);
ELC1(t,j)                =        ELC.l(j);
X1(t,i,j)                =        X.l(i,j);
TX1(t,j)                 =        TX.l(j);
Z1(t,j)                  =        Z.l(j);

ADEEI1(t,i)              =        ADEEI.l(i);

HOHincome1(t,l)          =        HOHincome.l(l);
Xp1(t,i,l)               =        Xp.l(i,l);
Xg1(t,i)                 =        Xg.l(i);
Xv1(t,i)                 =        Xv.l(i);
E1(t,i)                  =        E.l(i);
M1(t,i)                  =        M.l(i);
Q1(t,i)                  =        Q.l(i);
D1(t,i)                  =        D.l(i);

pke1(t,j)                =        pke.l(j);
pene1(t,j)               =        pene.l(j);
pfossil1(t,j)            =        pfossil.l(j);
psolid1(t,j)             =        psolid.l(j);
pnos1(t,j)               =        pnos.l(j);
pelc1(t,j)               =        pelc.l(j);
pref1(t,j)               =        pref.l(j);
px1(t,eni)               =        px.l(eni);
pf1(t,h,j)               =        pf.l(h,j);
pvae1(t,j)               =        pvae.l(j);
ptx1(t,j)                =        ptx.l(j);
pz1(t,j)                 =        pz.l(j);
pq1(t,i)                 =        pq.l(i);
pe1(t,i)                 =        pe.l(i);
pm1(t,i)                 =        pm.l(i);
pd1(t,i)                 =        pd.l(i);
epsilon1(t)              =        epsilon.l;

Sp1(t,l)                 =        Sp.l(l);
Sg1(t)                   =        Sg.l;
Td1(t,l)                 =        Td.l(l);
Tz1(t,j)                 =        Tz.l(j);
Tm1(t,i)                 =        Tm.l(i);

Cd1(t,j)                 =        Cd.l(j);
A1(t,i,j)                =        A.l(i,j);
Am1(t,i,j)               =        Am.l(i,j);
XLeontiefinverse1(t,i,j) =        XLeontiefinverse.l(i,j);
Ed1(t,j)                 =        Ed.l(j);
FD1(t,j)                 =        FD.l(j);
Eim_intermediate1(t,j)   =        Eim_intermediate.l(j);
Eim1(t,j)                =        Eim.l(j);
Yim1(t,j)                =        Yim.l(j);
EEP1(t)                  =        EEP.l;
EEC1(t)                  =        EEC.l;
EEE1(t)                  =        EEE.l;
EEI1(t)                  =        EEI.l;
EEB1(t)                  =        EEB.l;
PPI1(t,j)                =        PPI.l(j);
CPI1(t)                  =        CPI.l;
Renewable_rate1(t)       =        Renewable_rate.l;
Emissions1(t)            =        Emissions.l;
r_electricity1(t)        =        r_electricity.l;
TOTelectricity1(t)       =        TOTelectricity.l;
EV1(t,l)                 =        EV.l(l);
CV1(t,l)                 =        CV.l(l);
GDP1(t)                  =        GDP.l;
GDPCHK1(t)               =        GDPCHK.l;
UU1(t,l)                 =        UU.l(l);
SW1(t)                   =        SW.l;

C_AEEI1(t)               =        C_AEEI.l;
C_TFP1(t)                =        C_TFP.l;
C_ELE1(t)                =        C_ELE.l;

Walras1(t)               =        Walras.l;
Walras21(t)              =        Walras2.l;
FF1(t,l,h)               =        FF.l(l,h);
RFF1(t,l,h)              =        RFF.l(l,h);

* 4.3.2 Solve the dynamic model  -----------------------------------------------
* The following "if" command is only used for all CF scenario. If it is BAU scenario, or pre-BAU scenario for simuating dynamic parameter, we should not use them.
*if(        ord(t)  >= 3,
*         BAUDummy  =  0;
*          CFDummy  =  1;
*           ctr.lo  =  -inf;
*           ctr.up  =  +inf;
*       );

* You can set benchmark price changes to simulate the most real inflation
*pf.fx("LAB")      =  pf.l("LAB")*1.0;

FF.fx(l,"LAB")     =  FF.l(l,"LAB")*(1+gr);
CAPSTKg            =  CAPSTKg+sum(j,Xv.l(j)*(1-r_Sepcial)-Fg.l("CAP",j));
FFg.fx("CAP")      =  CAPSTKg*sum(j,depr(j)/card(j));
CAPSTKs(j)         =  CAPSTKs(j) +rInvestment(j)*sum(i,Xv.l(i)*r_Sepcial-Fs.l("CAP",i) ) ;

* The following equations is only used for DyanimcParameter.gms file
*C_TFP.lo=-inf;
*C_AEEI.lo=-inf;
*C_ELE.lo=-inf;
*C_TFP.up=+inf;
*C_AEEI.up=+inf;
*C_ELE.up=+inf;

* The following formula dynamically gives the exogenous dynamic parameter value of each recursive period of the model.
if(ord(t)=        1        ,C_TFP.fx=C_TFP1("2019");        C_AEEI.fx=C_AEEI1("2019");        C_ELE.fx=C_ELE1("2019");        pzbau(j)=pzbau1("2019",j);        pqbau(j)=pqbau1("2019",j);        xpbau(j,l)=xpbau1("2019",j,l);        );
if(ord(t)=        2        ,C_TFP.fx=C_TFP1("2020");        C_AEEI.fx=C_AEEI1("2020");        C_ELE.fx=C_ELE1("2020");        pzbau(j)=pzbau1("2020",j);        pqbau(j)=pqbau1("2020",j);        xpbau(j,l)=xpbau1("2020",j,l);        );
if(ord(t)=        3        ,C_TFP.fx=C_TFP1("2021");        C_AEEI.fx=C_AEEI1("2021");        C_ELE.fx=C_ELE1("2021");        pzbau(j)=pzbau1("2021",j);        pqbau(j)=pqbau1("2021",j);        xpbau(j,l)=xpbau1("2021",j,l);        );
if(ord(t)=        4        ,C_TFP.fx=C_TFP1("2022");        C_AEEI.fx=C_AEEI1("2022");        C_ELE.fx=C_ELE1("2022");        pzbau(j)=pzbau1("2022",j);        pqbau(j)=pqbau1("2022",j);        xpbau(j,l)=xpbau1("2022",j,l);        );
if(ord(t)=        5        ,C_TFP.fx=C_TFP1("2023");        C_AEEI.fx=C_AEEI1("2023");        C_ELE.fx=C_ELE1("2023");        pzbau(j)=pzbau1("2023",j);        pqbau(j)=pqbau1("2023",j);        xpbau(j,l)=xpbau1("2023",j,l);        );
if(ord(t)=        6        ,C_TFP.fx=C_TFP1("2024");        C_AEEI.fx=C_AEEI1("2024");        C_ELE.fx=C_ELE1("2024");        pzbau(j)=pzbau1("2024",j);        pqbau(j)=pqbau1("2024",j);        xpbau(j,l)=xpbau1("2024",j,l);        );
if(ord(t)=        7        ,C_TFP.fx=C_TFP1("2025");        C_AEEI.fx=C_AEEI1("2025");        C_ELE.fx=C_ELE1("2025");        pzbau(j)=pzbau1("2025",j);        pqbau(j)=pqbau1("2025",j);        xpbau(j,l)=xpbau1("2025",j,l);        );
if(ord(t)=        8        ,C_TFP.fx=C_TFP1("2026");        C_AEEI.fx=C_AEEI1("2026");        C_ELE.fx=C_ELE1("2026");        pzbau(j)=pzbau1("2026",j);        pqbau(j)=pqbau1("2026",j);        xpbau(j,l)=xpbau1("2026",j,l);        );
if(ord(t)=        9        ,C_TFP.fx=C_TFP1("2027");        C_AEEI.fx=C_AEEI1("2027");        C_ELE.fx=C_ELE1("2027");        pzbau(j)=pzbau1("2027",j);        pqbau(j)=pqbau1("2027",j);        xpbau(j,l)=xpbau1("2027",j,l);        );
if(ord(t)=        10        ,C_TFP.fx=C_TFP1("2028");        C_AEEI.fx=C_AEEI1("2028");        C_ELE.fx=C_ELE1("2028");        pzbau(j)=pzbau1("2028",j);        pqbau(j)=pqbau1("2028",j);        xpbau(j,l)=xpbau1("2028",j,l);        );
if(ord(t)=        11        ,C_TFP.fx=C_TFP1("2029");        C_AEEI.fx=C_AEEI1("2029");        C_ELE.fx=C_ELE1("2029");        pzbau(j)=pzbau1("2029",j);        pqbau(j)=pqbau1("2029",j);        xpbau(j,l)=xpbau1("2029",j,l);        );
if(ord(t)=        12        ,C_TFP.fx=C_TFP1("2030");        C_AEEI.fx=C_AEEI1("2030");        C_ELE.fx=C_ELE1("2030");        pzbau(j)=pzbau1("2030",j);        pqbau(j)=pqbau1("2030",j);        xpbau(j,l)=xpbau1("2030",j,l);        );
if(ord(t)=        13        ,C_TFP.fx=C_TFP1("2031");        C_AEEI.fx=C_AEEI1("2031");        C_ELE.fx=C_ELE1("2031");        pzbau(j)=pzbau1("2031",j);        pqbau(j)=pqbau1("2031",j);        xpbau(j,l)=xpbau1("2031",j,l);        );
if(ord(t)=        14        ,C_TFP.fx=C_TFP1("2032");        C_AEEI.fx=C_AEEI1("2032");        C_ELE.fx=C_ELE1("2032");        pzbau(j)=pzbau1("2032",j);        pqbau(j)=pqbau1("2032",j);        xpbau(j,l)=xpbau1("2032",j,l);        );
if(ord(t)=        15        ,C_TFP.fx=C_TFP1("2033");        C_AEEI.fx=C_AEEI1("2033");        C_ELE.fx=C_ELE1("2033");        pzbau(j)=pzbau1("2033",j);        pqbau(j)=pqbau1("2033",j);        xpbau(j,l)=xpbau1("2033",j,l);        );
if(ord(t)=        16        ,C_TFP.fx=C_TFP1("2034");        C_AEEI.fx=C_AEEI1("2034");        C_ELE.fx=C_ELE1("2034");        pzbau(j)=pzbau1("2034",j);        pqbau(j)=pqbau1("2034",j);        xpbau(j,l)=xpbau1("2034",j,l);        );
if(ord(t)=        17        ,C_TFP.fx=C_TFP1("2035");        C_AEEI.fx=C_AEEI1("2035");        C_ELE.fx=C_ELE1("2035");        pzbau(j)=pzbau1("2035",j);        pqbau(j)=pqbau1("2035",j);        xpbau(j,l)=xpbau1("2035",j,l);        );
if(ord(t)=        18        ,C_TFP.fx=C_TFP1("2036");        C_AEEI.fx=C_AEEI1("2036");        C_ELE.fx=C_ELE1("2036");        pzbau(j)=pzbau1("2036",j);        pqbau(j)=pqbau1("2036",j);        xpbau(j,l)=xpbau1("2036",j,l);        );
if(ord(t)=        19        ,C_TFP.fx=C_TFP1("2037");        C_AEEI.fx=C_AEEI1("2037");        C_ELE.fx=C_ELE1("2037");        pzbau(j)=pzbau1("2037",j);        pqbau(j)=pqbau1("2037",j);        xpbau(j,l)=xpbau1("2037",j,l);        );
if(ord(t)=        20        ,C_TFP.fx=C_TFP1("2038");        C_AEEI.fx=C_AEEI1("2038");        C_ELE.fx=C_ELE1("2038");        pzbau(j)=pzbau1("2038",j);        pqbau(j)=pqbau1("2038",j);        xpbau(j,l)=xpbau1("2038",j,l);        );
if(ord(t)=        21        ,C_TFP.fx=C_TFP1("2039");        C_AEEI.fx=C_AEEI1("2039");        C_ELE.fx=C_ELE1("2039");        pzbau(j)=pzbau1("2039",j);        pqbau(j)=pqbau1("2039",j);        xpbau(j,l)=xpbau1("2039",j,l);        );
if(ord(t)=        22        ,C_TFP.fx=C_TFP1("2040");        C_AEEI.fx=C_AEEI1("2040");        C_ELE.fx=C_ELE1("2040");        pzbau(j)=pzbau1("2040",j);        pqbau(j)=pqbau1("2040",j);        xpbau(j,l)=xpbau1("2040",j,l);        );
if(ord(t)=        23        ,C_TFP.fx=C_TFP1("2041");        C_AEEI.fx=C_AEEI1("2041");        C_ELE.fx=C_ELE1("2041");        pzbau(j)=pzbau1("2041",j);        pqbau(j)=pqbau1("2041",j);        xpbau(j,l)=xpbau1("2041",j,l);        );
if(ord(t)=        24        ,C_TFP.fx=C_TFP1("2042");        C_AEEI.fx=C_AEEI1("2042");        C_ELE.fx=C_ELE1("2042");        pzbau(j)=pzbau1("2042",j);        pqbau(j)=pqbau1("2042",j);        xpbau(j,l)=xpbau1("2042",j,l);        );
if(ord(t)=        25        ,C_TFP.fx=C_TFP1("2043");        C_AEEI.fx=C_AEEI1("2043");        C_ELE.fx=C_ELE1("2043");        pzbau(j)=pzbau1("2043",j);        pqbau(j)=pqbau1("2043",j);        xpbau(j,l)=xpbau1("2043",j,l);        );
if(ord(t)=        26        ,C_TFP.fx=C_TFP1("2044");        C_AEEI.fx=C_AEEI1("2044");        C_ELE.fx=C_ELE1("2044");        pzbau(j)=pzbau1("2044",j);        pqbau(j)=pqbau1("2044",j);        xpbau(j,l)=xpbau1("2044",j,l);        );
if(ord(t)=        27        ,C_TFP.fx=C_TFP1("2045");        C_AEEI.fx=C_AEEI1("2045");        C_ELE.fx=C_ELE1("2045");        pzbau(j)=pzbau1("2045",j);        pqbau(j)=pqbau1("2045",j);        xpbau(j,l)=xpbau1("2045",j,l);        );
if(ord(t)=        28        ,C_TFP.fx=C_TFP1("2046");        C_AEEI.fx=C_AEEI1("2046");        C_ELE.fx=C_ELE1("2046");        pzbau(j)=pzbau1("2046",j);        pqbau(j)=pqbau1("2046",j);        xpbau(j,l)=xpbau1("2046",j,l);        );
if(ord(t)=        29        ,C_TFP.fx=C_TFP1("2047");        C_AEEI.fx=C_AEEI1("2047");        C_ELE.fx=C_ELE1("2047");        pzbau(j)=pzbau1("2047",j);        pqbau(j)=pqbau1("2047",j);        xpbau(j,l)=xpbau1("2047",j,l);        );
if(ord(t)=        30        ,C_TFP.fx=C_TFP1("2048");        C_AEEI.fx=C_AEEI1("2048");        C_ELE.fx=C_ELE1("2048");        pzbau(j)=pzbau1("2048",j);        pqbau(j)=pqbau1("2048",j);        xpbau(j,l)=xpbau1("2048",j,l);        );
if(ord(t)=        31        ,C_TFP.fx=C_TFP1("2049");        C_AEEI.fx=C_AEEI1("2049");        C_ELE.fx=C_ELE1("2049");        pzbau(j)=pzbau1("2049",j);        pqbau(j)=pqbau1("2049",j);        xpbau(j,l)=xpbau1("2049",j,l);        );
if(ord(t)=        32        ,C_TFP.fx=C_TFP1("2050");        C_AEEI.fx=C_AEEI1("2050");        C_ELE.fx=C_ELE1("2050");        pzbau(j)=pzbau1("2050",j);        pqbau(j)=pqbau1("2050",j);        xpbau(j,l)=xpbau1("2050",j,l);        );
if(ord(t)=        33        ,C_TFP.fx=C_TFP1("2051");        C_AEEI.fx=C_AEEI1("2051");        C_ELE.fx=C_ELE1("2051");        pzbau(j)=pzbau1("2051",j);        pqbau(j)=pqbau1("2051",j);        xpbau(j,l)=xpbau1("2051",j,l);        );
if(ord(t)=        34        ,C_TFP.fx=C_TFP1("2052");        C_AEEI.fx=C_AEEI1("2052");        C_ELE.fx=C_ELE1("2052");        pzbau(j)=pzbau1("2052",j);        pqbau(j)=pqbau1("2052",j);        xpbau(j,l)=xpbau1("2052",j,l);        );
if(ord(t)=        35        ,C_TFP.fx=C_TFP1("2053");        C_AEEI.fx=C_AEEI1("2053");        C_ELE.fx=C_ELE1("2053");        pzbau(j)=pzbau1("2053",j);        pqbau(j)=pqbau1("2053",j);        xpbau(j,l)=xpbau1("2053",j,l);        );
if(ord(t)=        36        ,C_TFP.fx=C_TFP1("2054");        C_AEEI.fx=C_AEEI1("2054");        C_ELE.fx=C_ELE1("2054");        pzbau(j)=pzbau1("2054",j);        pqbau(j)=pqbau1("2054",j);        xpbau(j,l)=xpbau1("2054",j,l);        );
if(ord(t)=        37        ,C_TFP.fx=C_TFP1("2055");        C_AEEI.fx=C_AEEI1("2055");        C_ELE.fx=C_ELE1("2055");        pzbau(j)=pzbau1("2055",j);        pqbau(j)=pqbau1("2055",j);        xpbau(j,l)=xpbau1("2055",j,l);        );
if(ord(t)=        38        ,C_TFP.fx=C_TFP1("2056");        C_AEEI.fx=C_AEEI1("2056");        C_ELE.fx=C_ELE1("2056");        pzbau(j)=pzbau1("2056",j);        pqbau(j)=pqbau1("2056",j);        xpbau(j,l)=xpbau1("2056",j,l);        );
if(ord(t)=        39        ,C_TFP.fx=C_TFP1("2057");        C_AEEI.fx=C_AEEI1("2057");        C_ELE.fx=C_ELE1("2057");        pzbau(j)=pzbau1("2057",j);        pqbau(j)=pqbau1("2057",j);        xpbau(j,l)=xpbau1("2057",j,l);        );
if(ord(t)=        40        ,C_TFP.fx=C_TFP1("2058");        C_AEEI.fx=C_AEEI1("2058");        C_ELE.fx=C_ELE1("2058");        pzbau(j)=pzbau1("2058",j);        pqbau(j)=pqbau1("2058",j);        xpbau(j,l)=xpbau1("2058",j,l);        );
if(ord(t)=        41        ,C_TFP.fx=C_TFP1("2059");        C_AEEI.fx=C_AEEI1("2059");        C_ELE.fx=C_ELE1("2059");        pzbau(j)=pzbau1("2059",j);        pqbau(j)=pqbau1("2059",j);        xpbau(j,l)=xpbau1("2059",j,l);        );
if(ord(t)=        42        ,C_TFP.fx=C_TFP1("2060");        C_AEEI.fx=C_AEEI1("2060");        C_ELE.fx=C_ELE1("2060");        pzbau(j)=pzbau1("2060",j);        pqbau(j)=pqbau1("2060",j);        xpbau(j,l)=xpbau1("2060",j,l);        );

* solve the recursive dynamic model
Solve CEEEA using mcp;
);

*-------------------------------------------------------------------------------
* 5. Results output ------------------------------------------------------------
execute_unload 'output.gdx';
execute 'gdxxrw output.gdx output=C:\Users\jasperjia\Documents\gamsdir\P6BAU.xlsx Par=PPI1 RNG=PPI!A1 Par=CPI1 RNG=CPI!A1 Par=Energy1 RNG=Energy!A1 Par=TOT_Energy1 RNG=TOT_Energy!A1 Par=EM1 RNG=EM!A1 Par=PLC1 RNG=PLC!A1 Par=ctr1 RNG=ctr!A1 Par=VAE1 RNG=VAE!A1 Par=F1 RNG=F!A1 Par=ENE1 RNG=ENE!A1 Par=KE1 RNG=KE!A1 Par=ENE1 RNG=ENE!A1 Par=FOSSIL1 RNG=FOSSIL!A1 Par=SOLID1 RNG=SOLID!A1 Par=NOS1 RNG=NOS!A1 Par=REF1 RNG=REF!A1 Par=ELC1 RNG=ELC!A1 Par=X1 RNG=X!A1 Par=Z1 RNG=Z!A1 Par=ADEEI1 RNG=ADEEI!A1 Par=HOHincome1 RNG=HOHincome!A1 Par=Xp1 RNG=Xp!A1 Par=Xg1 RNG=Xg!A1 Par=Xv1 RNG=Xv!A1 Par=E1 RNG=E!A1 Par=M1 RNG=M!A1 Par=Q1 RNG=Q!A1 Par=D1 RNG=D!A1 Par=pke1 RNG=pke!A1 Par=pene1 RNG=pene!A1 Par=pfossil1 RNG=pfossil!A1 Par=psolid1 RNG=psolid!A1 Par=pnos1 RNG=pnos!A1 Par=pelc1 RNG=pelc!A1 Par=pref1 RNG=pref!A1 Par=px1 RNG=px!A1 Par=pf1 RNG=pf!A1 Par=pvae1 RNG=pvae!A1 Par=pz1 RNG=pz!A1 Par=pq1 RNG=pq!A1 Par=pe1 RNG=pe!A1 Par=pm1 RNG=pm!A1 Par=pd1 RNG=pd!A1 Par=epsilon1 RNG=epsilon!A1 Par=Sp1 RNG=Sp!A1 Par=Sg1 RNG=Sg!A1 Par=Td1 RNG=Td!A1 Par=Tz1 RNG=Tz!A1 Par=Tm1 RNG=Tm!A1 Par=Cd1 RNG=Cd!A1 Par=A1 RNG=A!A1 Par=Am1 RNG=Am!A1 Par=XLeontiefinverse1 RNG=XLeontiefinverse!A1 Par=Ed1 RNG=Ed!A1 Par=FD1 RNG=FD!A1 Par=Eim_intermediate1 RNG=Eim_intermediate!A1 Par=Eim1 RNG=Eim!A1 Par=Yim1 RNG=Yim!A1 Par=EEP1 RNG=EEP!A1 Par=EEC1 RNG=EEC!A1 Par=EEE1 RNG=EEE!A1 Par=EEI1 RNG=EEI!A1 Par=EEB1 RNG=EEB!A1 Par=Renewable_rate1 RNG=Renewable_rate!A1 Par=Emissions1 RNG=Emissions!A1 Par=r_electricity1 RNG=r_electricity!A1 Par=TOTelectricity1 RNG=TOTelectricity!A1 Par=EV1 RNG=EV!A1 Par=CV1 RNG=CV!A1 Par=GDP1 RNG=GDP!A1 Par=GDPCHK1 RNG=GDPCHK!A1 Par=UU1 RNG=UU!A1 Par=SW1 RNG=SW!A1 Par=C_AEEI1 RNG=C_AEEI!A1 Par=C_TFP1 RNG=C_TFP!A1 Par=C_ELE1 RNG=C_ELE!A1 Par=Walras1 RNG=Walras!A1 Par=Walras21 RNG=Walras2!A1 Par=FF1 RNG=FF!A1 Par=RFF1 RNG=RFF!A1';

* Only used in DynamicParameter.gms file
*execute_unload 'output-DynamicParameter.gdx';
