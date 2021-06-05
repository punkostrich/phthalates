#install.packages("mrgsolve")
library(mrgsolve) 

rm(list=ls())

HumanPK.code <- '
$PARAM @annotated
Vp                  :   23.5    : Volume of distribution parent compound L/kg bw
Vm                  :   23.5    : Volume of distribution metabolite compound 

ka                  :   1.6     : stomach absorption rate parent 1/h
ks   		            :   5     : skin Absorption rate to plasma (1/h)

BW                  :   82.3   : kg, Bodyweight (EPA Factors Handbook, 2011)

lambda_p            :   0.2     : lambda for parent 
lambda_zm           :   0.2    : lambda for metabolite 
lambda_u            :   0.004   : lambda for metabolite urine 
 
MW                  :   391     : g/mol, DEHP molecular mass 
Free                :   0.0001  : Free fraction of DEHP in plasma 
FreeM               :   0.007   : Free fraction of MEHP in plasma 

Kvoid               :   0.06974 : (L/hr), Daily urine volume rate (L/hr)
SA		              :   2.1     : skin area exposred m2
Cair                :   0.013   : ug/m3

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
capture Plasma_MEHP = CA_M;
capture AUC_CA = AUCCA;
capture AUC_CAM = AUCCA_M;
capture Urine = Aurine_M;

'


saveRDS(HumanPK.code,file='C:/Users/Punkostrich/Dropbox/NM/Class/DEHP/PK/DEHP/humanPK.RDS')

