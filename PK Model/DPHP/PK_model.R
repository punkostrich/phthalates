#install.packages("mrgsolve")
library(mrgsolve) 
############## DPHP #################
rm(list=ls())
################ DPHP ################
HumanPK.code <- '
$PARAM @annotated
Vp                  :   3.2    : Volume of distribution parent compound L/kg bw
Vm                  :   3.2    : Volume of distribution metabolite compound 
ka                  :   0.3    : stomach absorption rate parent 1/h
BW                  :   82.3   : kg, Bodyweight (EPA Factors Handbook, 2011)

lambda_p            :   0.12     : lambda for parent 
lambda_zm           :   0.17     : lambda for metabolite 
lambda_u            :   1e-4    : lambda for metabolite urine 
 
lambda_z_OH         :   0.27       : h-1
lambda_z_oxo        :   0.16       : h-1

lambda_u_OH         :   4.6e-3      : h-1
lambda_u_oxo        :   4.4e-3      : h-1

Fr_OH               :   0.5     : h-1
Fr_oxo              :   0.5    : h-1
V_OH                :   3.2    : Volume of distribution parent compound L/kg bw
V_oxo               :   3.2    : Volume of distribution metabolite compound 



$MAIN
double Vp_all    = Vp*BW;                                  // DPHP
double Vm_all    = Vm*BW;                                    // MPHP
double V_OH_all  = V_OH*BW;                                 // OH
double V_oxo_all = V_oxo*BW;                             // oxo
double lambda_m_MPHP  = lambda_zm    - lambda_u;         // metabolisam rate of metabolite
double lambda_m_OH    = lambda_z_OH  - lambda_u_OH;         // metabolisam rate of metabolite
double lambda_m_oxo   = lambda_z_oxo - lambda_u_oxo;         // metabolisam rate of metabolite


$CMT  AC_DPHP AST AC_MPHP AC_OH AC_oxo Aurine_MPHP Aurine_OH Aurine_oxo AUCCA_DPHP AUCCA_MPHP

$ODE
// Concentrations in the tissues and in the capillary blood of the tissues
double CA_DPHP = AC_DPHP/Vp_all;            // ug/L, Concentration of total DEHP in central compartment
double CA_MPHP = AC_MPHP/Vm_all;         // ug/L, Concentration of total MEHP in central compartment
double CA_OH   = AC_OH/V_OH_all;            // ug/L, Concentration of total DEHP in central compartment
double CA_oxo  = AC_oxo/V_oxo_all;         // ug/L, Concentration of total MEHP in central compartment


// {Stomach}
double RST = -ka*AST;                    // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AST = RST;                                  // mg, Amount in Stomach


// {DPHP distribution in central compartment}
double RC_DPHP = ka*AST - lambda_p*AC_DPHP;           // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_DPHP = RC_DPHP;                                  // ug, Amount in Stomach
dxdt_AUCCA_DPHP = CA_DPHP; // ug*h/L, Area under curve of DEHP in plasma compartment

// {MPHP distribution in central compartment}
double RC_MPHP = lambda_p*AC_DPHP- lambda_m_MPHP*AC_MPHP -lambda_u*AC_MPHP ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_MPHP = RC_MPHP;                                  // mg, Amount in central
dxdt_AUCCA_MPHP = CA_MPHP;               // ug*h/L, Area under curve of MEHP in plasma compartment

// {OH distribution in central compartment}
double RC_OH = Fr_OH * lambda_m_MPHP * AC_MPHP- lambda_m_OH * AC_OH -lambda_u_OH * AC_OH ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_OH = RC_OH;                                  // mg, Amount in central

// {oxo distribution in central compartment}
double RC_oxo = Fr_oxo * lambda_m_OH * AC_OH- lambda_m_oxo* AC_oxo -lambda_u_oxo * AC_oxo ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_oxo = RC_oxo;                                             // mg, Amount in central
   
                                  // mg, Amount in central
// Urine elimination
double RAurine_MPHP = lambda_u*AC_MPHP;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_MPHP = RAurine_MPHP;                //  ug, Amount in urine

double RAurine_OH = lambda_u_OH*AC_OH;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_OH = RAurine_OH;                //  ug, Amount in urine

double RAurine_oxo = lambda_u_oxo*AC_oxo;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_oxo = RAurine_oxo;                //  ug, Amount in urine



$TABLE
capture Urine_MPHP = Aurine_MPHP;
capture Urine_OH   = Aurine_OH;
capture Urine_oxo  = Aurine_oxo;
capture Plasma_DPHP = CA_DPHP;
capture Plasma_MPHP = CA_MPHP;
capture Plasma_OH   = CA_OH;
capture Plasma_oxo  = CA_oxo;
capture AUC_CA  = AUCCA_DPHP;
capture AUC_CAM = AUCCA_MPHP;

'


saveRDS(HumanPK.code,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DPHP/humanPK.RDS')

