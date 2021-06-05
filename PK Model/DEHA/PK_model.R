#install.packages("mrgsolve")
library(mrgsolve) 

########### DEHA #############
rm(list=ls())

HumanPK.code <- '
$PARAM @annotated
Vp                  :   5.3     : Volume of distribution parent compound L/kg bw
Vm                  :   0.5    : Volume of distribution metabolite compound 

ka                  :   0.2     : stomach absorption rate parent 1/h

BW                  :   82.3    : kg, Bodyweight (EPA Factors Handbook, 2011)

lambda_p            :   0.08    : lambda for parent (h-1)


lambda_zm           :   0.2       : h-1
lambda_z_OH         :   0.4       : h-1
lambda_z_cx         :   0.4       : h-1
lambda_z_oxo        :   0.4       : h-1


lambda_u_OH         :   2.8E-4     : h-1
lambda_u_cx         :   8E-4     : h-1
lambda_u_oxo        :   2E-4     : h-1


Fr_OH               :   0.5     : h-1
Fr_oxo              :   0.5     : h-1

$MAIN
double Vp_all = Vp*BW;
double Vm_all = Vm*BW;
double lambda_m_OH    = lambda_z_OH - lambda_u_OH;         // metabolisam rate of metabolite
double lambda_m_cx    = lambda_z_cx - lambda_u_cx;         // metabolisam rate of metabolite
double lambda_m_oxo   = lambda_z_oxo - lambda_u_oxo;         // metabolisam rate of metabolite


$CMT  AC_DEHA AST AC_MEHA AC_OH AC_cx AC_oxo Aurine_OH Aurine_cx Aurine_oxo AUCCA AUCCA_M

$ODE
double CA_DEHA = AC_DEHA/Vp_all;            // ug/L, Concentration of total DEHP in central compartment
double CA_M = AC_MEHA/Vm_all;         // ug/L, Concentration of total MEHP in central compartment
double CA_cx = AC_cx/Vm_all;         // ug/L, Concentration of total MEHP in central compartment

// {Stomach}
double RST = -ka*AST;                    // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AST = RST;                                  // mg, Amount in Stomach


// {DEHA distribution in central compartment}
double RC_DEHA = ka*AST - lambda_p*AC_DEHA;           // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_DEHA = RC_DEHA;                                  // ug, Amount in Stomach
dxdt_AUCCA = CA_DEHA; // ug*h/L, Area under curve of DEHP in plasma compartment


// {MEHA distribution in central compartment}
double RC_MEHA = lambda_p*AC_DEHA - lambda_zm * AC_MEHA;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_MEHA = RC_MEHA;                                  // mg, Amount in central
dxdt_AUCCA_M = CA_M; // ug*h/L, Area under curve of DEHP in plasma compartment


// {OH distribution in central compartment}
double RC_OH = Fr_OH * lambda_zm * AC_MEHA- lambda_m_OH * AC_OH -lambda_u_OH * AC_OH ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_OH = RC_OH;                                  // mg, Amount in central

// {cx distribution in central compartment}
double RC_cx = (1 - Fr_OH) * lambda_zm * AC_MEHA- lambda_m_cx* AC_cx -lambda_u_cx * AC_cx ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_cx = RC_cx;                                             // mg, Amount in central
//dxdt_AUCCA_cx = CA_cx;               // ug*h/L, Area under curve of MEHP in plasma compartment

// {oxo distribution in central compartment}
double RC_oxo = Fr_oxo * lambda_m_OH * AC_OH- lambda_m_oxo* AC_oxo -lambda_u_oxo * AC_oxo ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_oxo = RC_oxo;                                             // mg, Amount in central
                                  // mg, Amount in central


// Urine elimination
double RAurine_OH = lambda_u_OH * AC_OH;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_OH = RAurine_OH;                //  ug, Amount in urine

double RAurine_cx = lambda_u_cx * AC_cx;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_cx = RAurine_cx;                //  ug, Amount in urine

double RAurine_oxo = lambda_u_oxo * AC_oxo;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_oxo = RAurine_oxo;                //  ug, Amount in urine



$TABLE
capture Urine_OH = Aurine_OH;
capture Urine_cx = Aurine_cx;
capture Urine_oxo = Aurine_oxo;
capture Plasma_DEHA = CA_DEHA;
capture Plasma_cx = CA_cx;
capture AUC_CA = AUCCA;
capture AUC_CAM = AUCCA_M;

'


saveRDS(HumanPK.code,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHA/humanPK.RDS')
