#install.packages("mrgsolve")
library(mrgsolve) 

######### DINCH #############
rm(list=ls())

HumanPK.code <- '
$PARAM @annotated
Vp                  :   5.3     : Volume of distribution parent compound L/kg bw
Vm                  :   0.05    : Volume of distribution metabolite compound 

ka                  :   0.2     : stomach absorption rate parent 1/h

BW                  :   82.3    : kg, Bodyweight (EPA Factors Handbook, 2011)

lambda_p            :   0.08    : lambda for parent (h-1)


lambda_zm           :   0.2     : h-1
lambda_u            :   0.1     : h-1

lambda_z_MHNCH      :   1       : h-1
lambda_z_cx         :   1       : h-1
lambda_z_oxo        :   1       : h-1
lambda_z_CHDA       :   1       : h-1


lambda_u_MHNCH      :   0.1     : h-1
lambda_u_cx         :   0.2     : h-1
lambda_u_oxo        :   0.2     : h-1
lambda_u_CHDA       :   1       : h-1

Fr_MINCH            :   0.5     : h-1
Fr_MHNCH            :   0.5     : h-1
Fr_cx               :   0.17    : h-1
Fr_oxo              :   0.13    : h-1


$MAIN
double Vp_all = Vp*BW;
double Vm_all = Vm*BW;
double lambda_m_MINCH = lambda_zm - lambda_u;         // metabolisam rate of metabolite
double lambda_m_MHNCH = lambda_z_MHNCH - lambda_u_MHNCH;         // metabolisam rate of metabolite
double lambda_m_cx    = lambda_z_cx - lambda_u_cx;         // metabolisam rate of metabolite
double lambda_m_oxo   = lambda_z_oxo - lambda_u_oxo;         // metabolisam rate of metabolite
double lambda_m_CHDA  = lambda_z_CHDA - lambda_u_CHDA;         // metabolisam rate of metabolite


$CMT  AC_DINCH AST AC_MINCH AC_MHNCH AC_cx AC_oxo AC_CHDA Aurine_MINCH Aurine_MHNCH Aurine_cx Aurine_oxo Aurine_CHDA AUCCA_DINCH AUCCA_MINCH 

$ODE
double CA_DINCH = AC_DINCH/Vp_all;            // ug/L, Concentration of total DEHP in central compartment
double CA_MINCH = AC_MINCH/Vm_all;         // ug/L, Concentration of total MEHP in central compartment

// {Stomach}
double RST = -ka*AST;                    // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AST = RST;                                  // mg, Amount in Stomach


// {DINCH distribution in central compartment}
double RC_DINCH = ka*AST - lambda_p*AC_DINCH;           // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_DINCH = RC_DINCH;                                  // ug, Amount in Stomach
dxdt_AUCCA_DINCH = CA_DINCH; // ug*h/L, Area under curve of DEHP in plasma compartment


// {MINCH distribution in central compartment}
double RC_MINCH = lambda_p*AC_DINCH- lambda_m_MINCH*AC_MINCH -lambda_u*AC_MINCH ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_MINCH = RC_MINCH;                                  // mg, Amount in central
dxdt_AUCCA_MINCH = CA_MINCH;               // ug*h/L, Area under curve of MEHP in plasma compartment

// {MHNCH distribution in central compartment}
double RC_MHNCH = Fr_MHNCH * lambda_m_MINCH * AC_MINCH- lambda_m_MHNCH * AC_MHNCH -lambda_u_MHNCH * AC_MHNCH ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_MHNCH = RC_MHNCH;                                  // mg, Amount in central

// {cx distribution in central compartment}
double RC_cx = Fr_cx * lambda_m_MINCH * AC_MINCH- lambda_m_cx* AC_cx -lambda_u_cx * AC_cx ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_cx = RC_cx;                                             // mg, Amount in central

// {oxo distribution in central compartment}
double RC_oxo = Fr_oxo * lambda_m_MHNCH * AC_MHNCH- lambda_m_oxo* AC_oxo -lambda_u_oxo * AC_oxo ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_oxo = RC_oxo;                                             // mg, Amount in central
                                  // mg, Amount in central

// {CHDA distribution in central compartment}
double RC_CHDA = (1-Fr_oxo) * lambda_m_MHNCH * AC_MHNCH + lambda_m_oxo* AC_oxo + lambda_m_cx* AC_cx + (1 - Fr_MHNCH - Fr_cx) * lambda_m_MINCH * AC_MINCH - lambda_u_CHDA *AC_CHDA ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC_CHDA = RC_CHDA;                                  // mg, Amount in central


// Urine elimination
double RAurine_MINCH = lambda_u*AC_MINCH;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_MINCH = RAurine_MINCH;                //  ug, Amount in urine

double RAurine_MHNCH = lambda_u_MHNCH*AC_MHNCH;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_MHNCH = RAurine_MHNCH;                //  ug, Amount in urine

double RAurine_cx = lambda_u_cx*AC_cx;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_cx = RAurine_cx;                //  ug, Amount in urine

double RAurine_oxo = lambda_u_oxo*AC_oxo;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_oxo = RAurine_oxo;                //  ug, Amount in urine

double RAurine_CHDA = lambda_u_CHDA*AC_CHDA;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_CHDA = RAurine_CHDA;                //  ug, Amount in urine


$TABLE
capture Urine_MINCH = Aurine_MINCH;
capture Urine_MHNCH = Aurine_MHNCH;
capture Urine_cx = Aurine_cx;
capture Urine_oxo = Aurine_oxo;
capture Urine_CHDA = Aurine_CHDA;
capture Plasma_DINCH = CA_DINCH;
capture Plasma_MINCH = CA_MINCH;
capture AUC_CA = AUCCA_DINCH;
capture AUC_CAM = AUCCA_MINCH;

'


saveRDS(HumanPK.code,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DINCH/humanPK.RDS')

