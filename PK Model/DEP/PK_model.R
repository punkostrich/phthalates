#install.packages("mrgsolve")
library(mrgsolve) 

rm(list=ls())

HumanPK.code <- '
$PARAM @annotated
Vp                  :   17     : Volume of distribution parent compound L/kg bw
Vm                  :   0.375  : Volume of distribution metabolite compound  L/kg bw

ka                  :   2       : stomach absorption rate parent 1/h

BW                  :   82.3    : kg, Bodyweight (EPA Factors Handbook, 2011)

lambda_p            :   0.7     : lambda for parent (h-1)
lambda_zm           :   0.74    : lambda for metabolite (h-1)
lambda_u            :   0.6     : lambda for metabolite urine (h-1)
 
MW                  :   222     : g/mol, DEP molecular mass 
Free                :   0.0001  : Free fraction of DEP in plasma 
FreeM               :   0.007   : Free fraction of MEP in plasma 


$MAIN
double Vp_all = Vp*BW;
double Vm_all = Vm*BW;
double lambda_m = lambda_zm - lambda_u;         // metabolisam rate of metabolite


$CMT  AC AST ACM AUCCA AUCCA_M  Aurine_M 

$ODE
// Concentrations in the tissues and in the capillary blood of the tissues
double CA_free = AC/Vp_all*Free;      // ug/L, Free DEHP oncentration in the central compartment
double CA = CA_free/Free;            // ug/L, Concentration of total DEHP in central compartment
double CA_freeM = ACM/Vm_all*FreeM;    // ug/L, Free MEHP oncentration in the central compartment
double CA_M = CA_freeM/FreeM;         // ug/L, Concentration of total MEHP in central compartment

// {Stomach}
double RST = -ka*AST;                    // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AST = RST;                                  // mg, Amount in Stomach


// {DEHP distribution in central compartment}
double RC = ka*AST - lambda_p*AC;           // ug/h,   Rate of chagne in stomach caomprtment
dxdt_AC = RC;                                  // ug, Amount in Stomach
dxdt_AUCCA = CA; // ug*h/L, Area under curve of DEHP in plasma compartment


// {MEHP distribution in central compartment}
double RCM = lambda_p*AC- lambda_m*ACM -lambda_u*ACM ;          // ug/h,   Rate of chagne in stomach caomprtment
dxdt_ACM = RCM;                                  // mg, Amount in Stomach
dxdt_AUCCA_M = CA_M;                // ug*h/L, Area under curve of MEHP in plasma compartment


// Urine elimination
double RAurine = lambda_u*ACM;          // ug/h, Rate of urine clearance rate of MEHP
dxdt_Aurine_M = RAurine;                //  ug, Amount in urine


$TABLE
capture Plasma = CA;
capture Plasma_MEP = CA_M;
capture AUC_CA = AUCCA;
capture AUC_CAM = AUCCA_M;
capture Urine = Aurine_M;

'


saveRDS(HumanPK.code,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEP/humanPK.RDS')

