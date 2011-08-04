#include<stdio.h>
#include<conio.h>
#include<math.h>
#include<iostream.h>
#include<fstream.h>

int choice;
int i,n,j,count,count2;
int aa,bb,cc,dd,ee,ff,gg;
int eoftest;
int eoftest2;
int display;
int istate;
int state[20];
int tempi;
int tempj;
int tempk;
int templ;
int points[12];
int savetype;
int ii,nn,jj,kk,ll;
int absorber;
int boiler;
int refrig;
int source;
int turbine;
int flag;
int index;
int single;
int eofsort;
int sorttempi;
int oldsys;
int uncert;
int cycletype;	//1 designates cycle with partial boiling (desired), 2 all vapor
				//and 3 no vapor formed

long double sum1=0,sum2=0;
long double limit;
double factor;	//factor to include refrigeration in 2nd law analysis
double GC;		//GC area ratio found in calibration
double incr;
double incr2;
double incr3;
double hst;		//heat source temperature
double bp;		//boiler pressure
double bt;		//boiler exit temperature
double at;		//ambient temperature
double abt;		//absorber temperature
double ap;		//absorber pressure
double bsmf;	//basic solution mass fraction
double hsfr;	//heat source flow ratio...m(hs)/m(str)
double upperhst;	//hst used in range calculations
double upperbp;		//bp used in range calculations
double upperbt;		//bt used in range calculations
double upperat;		//at used in range calculations
double upperabt;	//abt used in range calculations
double upperap;		//ap used in range calculations
double upperbsmf;	//bsmf used in range calculations
double upperhsfr;	//hsfr used in range calculations
double uppertie;	//tie used in range calculations
double upperpie;	//pie used in range calculations
double upperrecHE;	//recHE used in range calculations
double upperboilerHE;	//boilerHE used in range calculations
double bufhst;
double bufbp;
double bufbt;
double bufat;
double bufabt;
double bufap;
double bufbsmf;
double bufhsfr;
double buftie;
double bufpie;
double bufrecHE;
double bufboilerHE;
double step[12];
double T[20];
double P[20];
double x[20];
double m[20];
double xw,xv;		//weak and vapor mfs expected based on measured boiler T,P,x
double calcT[20];
double calcP[20];
double calcx[20];
double calcm[20];
double thermo[25];	//used when pulling raw thermocouple values from data files
double trans[8];	//used when pulling raw transducer values from data files
double designT[20];
double designP[20];
double designx[20];
double designm[20];
double lossT[10];	//temperature approach limits
double lossP[10];	//pressure losses (%)
double tie;			//turbine isentropic efficiency
double pie;			//pump isentropic efficiency
double xNH3l[20],xNH3v[20];		//mass fraction of ammonia as liquid and vapor
double xH2Ol[20],xH2Ov[20];		//mass fraction of water as liquid and vapor
								//note, the above 4 mass fractions add to 1
double guess;
double h[20];		//enthalpy matrix
double s[20];		//entropy matrix
double v[20];		//specific volume matrix
double av[20];		//availability matrix
double Q[20];		//quality of the mixture
double temp1,temp2,temp3,temp4,temp5;	//all purpose junk variables
double temp,tempx,tempT,tempP,temps,temph,tempv;
double tempxb,tempTb,tempPb,tempsb,temphb,tempvb;
double tempxh,tempTh,tempPh,tempsh,temphh,tempvh;
double tempxl,tempTl,tempPl,tempsl,temphl,tempvl;
double tempxv,tempTv,tempPv,tempsv,temphv,tempvv;
double qboiler;		//boiler heat input per unit strong
double Qboiler;		//boiler heat input
double wturbine;	//turbine work output per unit strong
double Wturbine;	//turbine work output
double qabsorber;	//heat removal from absorber per unit strong
double Qabsorber;	//heat removal from absorber
double Qrecmax;		//heat recovery for 100% effectiveness
double Qboilermax;	//boiler heat input for 100% effectiveness
double Qrec;		//recovery unit heat transfer
double qrec;		//recovery per unit bs
double Wpump;		//pump work
double wpump;		//pump work per unit strong solution
double Qrefrig;		//refrigeration
double qrefrig;		//refrigeration per unit strong solution
double Irecovery;
double Irefrig;
double Iboiler;
double Iturbine;
double Ipump;
double Iabsorber;
double Isep;
double Ithrottle;
double firstlaw;	//1st law efficiency
double secondlaw;	//2nd law efficiency
double error2;
double sortvalue;
double sortlast;
double sortbt;
double sortbp;
double sortabt;
double sortap;
double sortbsmf;
double sorthst;
double sorthsfr;
double buffdate;
double upperfdate;
double sortfdate;
double fdate;
double day;
double month;
double year;
double Time;
double recHE;		//recovery heat exchanger effectiveness
double boilerHE;	//boiler heat exchanger effectiveness
double recHEdes;	//used in comaparison of exp to theor in results page
double boilerHEdes;
//uncertainties
double UrecHE;
double UboilerHE;
double UQboiler;
double UQabsorber;
double UQrec;
double UQ;			//uncertainty in the quality calculated by properties from measured T,P,x
double Uxw;			//uncertainty in the weak mf calculated by props from measured boiler T,P,x
double Uxv;			//uncertainty in the vapor mf calculated by props from measured boiler T,P,x
double UWpump;
double UWturbine;
double UQrefrig;
double Uqboiler;
double Uqabsorber;
double Uqrec;
double Uwpump;
double Uwturbine;
double Uqrefrig;
double Ufirstlaw;
double Usecondlaw;
double Uvaporfrac;
double UP[20];		//uncertainty of pressures
double UT[20];		//uncertainty of temperatures
double Um[20];		//uncertainty of mass flow rates
double Ux[20];		//uncertainty of mass fractions
double UGC;			//uncertainty of the calibration constant
double Tlowcal[20]; //calibration values at low T - with Tlowcal[0]=0 C actual
double Thighcal[20];//calibration values at high T - with Thighcal[0]=100 C actual
double UTreading[20];	//uncertainty in the thermocouple readings
double Plowcal[10];		//calibration value at atmospheric pressure
double UPreading[10];	//uncertainty in the pressure transducer readings
double Umreading[10];	//uncertainty in the flow meter readings (all 2% full scale)
double Uxvapor;		//uncertainty in GC vapor readings
double Uxliquid;	//uncertainty in GC liquid readings
double NPreading[10];	//number of pressure readings (note 1 gauge reading = 1 sample)
						//and 1 transducer reading = 250 samples
double NTreading[20];	//number of 250 sample thermocouple readings
double Nmreading[10];	//number of mass flow meter readings
double Nxstrong;		//number of strong samples used to find x
double Nxweak;			//number of weak samples used to find x
double Nxvapor;			//number of vapor samples used to find x
double calcNPreading[10];	//number of pressure readings (note 1 gauge reading = 1 sample)
						//and 1 transducer reading = 250 samples
double calcNTreading[20];	//number of 250 sample thermocouple readings
double calcNmreading[10];	//number of mass flow meter readings
double calcNxstrong;		//number of strong samples used to find x
double calcNxweak;			//number of weak samples used to find x
double calcNxvapor;			//number of vapor samples used to find x
double UPFS[10];			//full scale reading of gauge/transducer
double UmFS[10];			//full scale reading of flow meter

//partial derivatives
double dQdT4,dQdP4,dQdx1;
double dxwdT4,dxwdP4,dxwdx1;
double dxvdT4,dxvdP4,dxvdx1;
double dQhdT15,dQhdT16,dQhdP15,dQhdm15;
double dqhdT15,dqhdT16,dqhdP15,dqhdm15,dqhdm1;
double dQrecdT8,dQrecdT10,dQrecdP4,dQrecdP10,dQrecdx8,dQrecdm8;
double dqrecdT8,dqrecdT10,dqrecdP4,dqrecdP10,dqrecdx8,dqrecdm8,dqrecdm1;
double dQcdT5,dQcdT17,dQcdP4,dQcdP7,dQcdx5,dQcdm5;
double dqcdT5,dqcdT17,dqcdP4,dqcdP7,dqcdx5,dqcdm5,dqcdm1;
double dQadT1,dQadT5,dQadT10,dQadP1,dQadP4,dQadP7,dQadP10;
double dqadT1,dqadT5,dqadT10,dqadP1,dqadP4,dqadP7,dqadP10;
double dQadx1,dQadx5,dQadx8,dQadm1,dQadm5,dQadm8;
double dqadx1,dqadx5,dqadx8,dqadm1,dqadm5,dqadm8;
double dWpdT1,dWpdP1,dWpdP4,dWpdx1,dWpdm1;
double dwpdT1,dwpdP1,dwpdP4,dwpdx1,dwpdm1;
double dWtdT5,dWtdP4,dWtdP7,dWtdx5,dWtdm5;
double dwtdT5,dwtdP4,dwtdP7,dwtdx5,dwtdm5,dwtdm1;
double dvfdm1,dvfdm5;
double dHEbdT3,dHEbdT15,dHEbdT16,dHEbdP4,dHEbdP15,dHEbdx1,dHEbdm15;
double dHErdT1,dHErdT8,dHErdT10,dHErdP1,dHErdP4,dHErdP10,dHErdx1,dHErdx8,dHErdm8;
double d1stdT1,d1stdT5,d1stdT15,d1stdT16,d1stdT17;
double d1stdP1,d1stdP4,d1stdP7,d1stdP15;
double d1stdx1,d1stdx5,d1stdm1,d1stdm5,d1stdm15;
double d2nddT1,d2nddT5,d2nddT15,d2nddT16,d2nddT17;
double d2nddP1,d2nddP4,d2nddP7,d2nddP15;
double d2nddx1,d2nddx5,d2nddm1,d2nddm5,d2nddm15;


char filename[80];	//filename root with .txt extension
char fname[20];		//filename root
char rfilename[80];	//filename root with r.txt extension
char rangefile[80];	//range values filename
char resultsfile[80];	//results filename to send values to
char rawfile[80];	//data file with all experimental data
char propfile[80];	//full address properties filename
char Tunits[10];	//temperature units
char Punits[10];	//pressure units
char xunits[10];	//mass fraction units
char munits[10];	//mass flow units
char hunits[10];	//enthalpy units
char sunits[10];	//entropy units
char vunits[10];	//specific volume units
char tempTunits[10];
char tempPunits[10];
char tempmunits[10];
char Tunitstemp[10];
char Punitstemp[10];
char schoice[80];	//string choice
char svalue[20];	//temp file string values
char buf[256];		//temp string
char buf2[80];		//temp string
char buf3[80];		//temp string
char buf4[80];		//temp string
char comments[71];	//file list comments
char error[80];
char fluid[30];
char function[30];
char activeh[3];
char activeP[3];
char actives[3];
char activeT[3];
char activev[3];
char activex[3];
char sortby[10];
char sortfile[80];
char date[10];
char lastfile[10];
char rlastfile[10];
char histfile[20];

//variables used in property evaluation
double TT,PP,xx;	//temp (K), pressure (bar), mass fraction of ammonia
double qm;			//quality of the mixture in the saturation region
double q;			//quality of a pure component
double y;			//mole fraction of ammonia
double Tc,Pc;		//critical temp and pressure of mixture
double Tb,Td;		//bubble and dew point temps of mixture
double Tb1,Td1;		//bubble and dew point temps of mixture
double Pb,Pd;		//bubble and dew point pressures of mixture
double sb,sd;		//bubble and dew point entropies
double hb,hd;		//bubble and dew point enthalpies
double vb,vd;		//bubble and dew point spec volumes
double xb,xd;		//bubble and dew point amm mass fractions of mixture
double yb,yd;		//bubble and dew point amm mole fractions of mixture
double M;			//molecular weight of mixture
double Mb,Md;		//molecular weight of mixture at bubble and dew points
double Tsatw;		//saturation temperature of water at P (K)
double Tsata;		//saturation temperature of ammonia at P (K)
double Tcw=1165.1;	//critical temperature of water (R)
double Pcw=3208;	//critical pressure of water (psi)
double Tca=729.9;	//critical temperature of ammonia (R)
double Pca=1646;	//critical pressure of ammonia (psi)
double Ptpa=.04;      	//triple point of ammonia (bar)
double Ttpa=197;	//triple point of ammonia (K)
double Ptpw=.006113;   //triple point of water (bar)
double Ttpw=273.16;	//triple point of water (K)
double Tr,Pr;		//reduced temperature and pressure
double TB=100;		//reference temperature (K)
double PB=10;		//reference pressure (bar)
double R=8.31451;	//gas constant (kJ/kmolK)
double hLw,hLa;		//enthalpy of liquid water and ammonia
double hgw,hga;		//enthalpy of gaseous water and ammonia
double sLw,sLa;		//entropy of liquid water and ammonia
double sgw,sga;		//entropy of gaseous water and ammonia
double vLw,vLa;		//specific volume of liquid water and ammonia
double vgw,vga;		//specific volume of gaseous water and ammonia
double hE;			//excess enthalpy of mixture
double sE;			//excess entropy of mixture
double vE;			//excess specific volume of mixture
double smix;		//entropy of mixing
double hLm;			//liq enthalpy of mixture
double sLm;			//liq entropy of mixture
double vLm;			//liq spec volume of mixture
double hgm;			//gas enthalpy of mixture
double sgm;			//gas entropy of mixture
double vgm;			//gas spec volume of mixture
double hw,ha,hm;	//enthalpy of water,ammonia,mixture
double sw,sa,sm;	//entropy of water,ammonia,mixture
double vw,va,vm;	//specific volume of water,ammonia,mixture
double xxNH3l,xxNH3v;
double xxH2Ol,xxH2Ov;


//matrix constants, with values defined in main program
long double a[5];
long double b[9];
long double Ci[8];
long double ai[7];
long double Aij[7][5];
long double Cij[8][11];

//functions
double bubbleT(double P,double x);						//K
double dewT(double P,double x);							//K
double bubbleP(double T,double x);						//bar
double dewP(double T,double x);							//bar
double bubblex(double T,double P);						//mass NH3/mass total
double dewx(double T,double P);							//mass NH3/mass total
double crittemp(double x);								//F
double critpress(double x);								//psi
double enthalpylw(double Tr,double Pr);					//kJ/kmol
double enthalpyla(double Tr,double Pr);					//kJ/kmol
double entropylw(double Tr,double Pr);					//kJ/kmolK
double entropyla(double Tr,double Pr);					//kJ/kmolK
double specvollw(double Tr,double Pr);					//m3/kmol
double specvolla(double Tr,double Pr);					//m3/kmol
double enthalpygw(double Tr,double Pr);					//kJ/kmol
double enthalpyga(double Tr,double Pr);					//kJ/kmol
double entropygw(double Tr,double Pr);					//kJ/kmolK
double entropyga(double Tr,double Pr);					//kJ/kmolK
double specvolgw(double Tr,double Pr);					//m3/kmol
double specvolga(double Tr,double Pr);					//m3/kmol
double enthalpyE(double Tr,double Pr,double y);			//kJ/kmol
double entropyE(double Tr,double Pr,double y);			//kJ/kmolK
double specvolE(double Tr,double Pr,double y);			//m3/kmol
double entropymixing(double y);							//kJ/kmolK
double enthalpyL(double hLa,double hLw,double hE,double y);	//kJ/kmol
double entropyL(double sLa,double sLw,double sE,double smix,double y);	//kJ/kmolK
double specvolL(double vLa,double vLw,double vE,double y);		//m3/kmol
double enthalpyg(double hga,double hgw,double y);				//kJ/kmol
double entropyg(double sga,double sgw,double smix,double y);	//kJ/kmolK
double specvolg(double vga,double vgw,double y);				//m3/kmol

//some constants for ammmonia and water
long double A1a=3.971423e-2,A1w=2.748796e-2;
long double A2a=-1.790557e-5,A2w=-1.016665e-5;
long double A3a=-1.308905e-2,A3w=-4.452025e-3;
long double A4a=3.752836e-3,A4w=8.389246e-4;
long double B1a=1.634519e+1,B1w=1.214557e+1;
long double B2a=-6.508119,B2w=-1.898065;
long double B3a=1.448937,B3w=2.911966e-1;
long double C1a=-1.049377e-2,C1w=2.136131e-2;
long double C2a=-8.288224,C2w=-3.169291e+1;
long double C3a=-6.647257e+2,C3w=-4.634611e+4;
long double C4a=-3.045352e+3,C4w=0.0;
long double D1a=3.673647,D1w=4.019170;
long double D2a=9.989629e-2,D2w=-5.175550e-2;
long double D3a=3.617622e-2,D3w=1.951939e-2;
long double hLroa=4.878573,hLrow=21.821141;
long double hgroa=26.468873,hgrow=60.965058;
long double sLroa=1.644773,sLrow=5.733498;
long double sgroa=8.339026,sgrow=13.453430;

double Troa=3.2252,Trow=5.0705;
double Proa=2.0000,Prow=3.0000;

long double E1 = -41.733398;
long double E2 = 0.02414;
long double E3 = 6.702285;
long double E4 = -0.011475;
long double E5 = 63.608967;
long double E6 = -62.490768;
long double E7 = 1.761064;
long double E8 = 0.008626;
long double E9 = 0.387983;
long double E10 = -0.004772;
long double E11 = -4.648107;
long double E12 = 0.836376;
long double E13 = -3.553627;
long double E14 = 0.000904;
long double E15 = 24.361723;
long double E16 = -20.736547;

//alternate values - these were found somewhere
//long double E1 = -46.26129;
//long double E2 = 0.02060225;
//long double E3 = 7.292369;
//long double E4 = -0.01032613;
//long double E5 = 80.74824;
//long double E6 = -84.61214;
//long double E7 = 24.52882;
//long double E8 = 0.009598767;
//long double E9 = -1.475383;
//long double E10 = -0.005038107;
//long double E11 = -96.40398;
//long double E12 = 122.6973;
//long double E13 = -7.582637;
//long double E14 = 0.0006012445;
//long double E15 = 54.87018;
//long double E16 = -76.67596;

//enthalpy of pure water in liquid phase	OK
double enthalpylw(double Tr,double Pr) {
	hLw=-R*TB*(-hLrow+B1w*(Trow-Tr)+B2w/2*(Trow*Trow-Tr*Tr)+B3w/3*
		(pow(Trow,3)-pow(Tr,3))-(A1w+A4w*Tr*Tr)*(Pr-Prow)-
		A2w/2*(Pr*Pr-Prow*Prow));
	return hLw;		// kJ/kmol(water)
}

//enthalpy of pure ammonia in liquid phase	OK
double enthalpyla(double Tr,double Pr) {
	hLa=-R*TB*(-hLroa+B1a*(Troa-Tr)+B2a/2*(Troa*Troa-Tr*Tr)+B3a/3*
		(pow(Troa,3)-pow(Tr,3))-(A1a+A4a*Tr*Tr)*(Pr-Proa)-
		A2a/2*(Pr*Pr-Proa*Proa));
	return hLa;		// kJ/kmol(ammonia)
}

//entropy of pure water in liquid phase		OK
double entropylw(double Tr,double Pr) {
	sLw=-R*(-sLrow-B1w*log(Tr/Trow)+B2w*(Trow-Tr)+
		B3w/2*(Trow*Trow-Tr*Tr)+(Pr-Prow)*(A3w+2*A4w*Tr));
	return sLw;		//kJ/kmol K (water)
}

//entropy of pure ammonia in liquid phase	OK
double entropyla(double Tr,double Pr) {
	sLa=-R*(-sLroa-B1a*log(Tr/Troa)+B2a*(Troa-Tr)+
		B3a/2*(Troa*Troa-Tr*Tr)+(Pr-Proa)*(A3a+2*A4a*Tr));
	return sLa;
}

//specific volume of water in liquid phase	OK
double specvollw(double Tr,double Pr) {
	vLw=R*TB/PB*(A1w+A3w*Tr+A4w*Tr*Tr+A2w*Pr);
	vLw=vLw/100;		//m3/kmol (water)
	return vLw;
}

//specific volume of ammonia in liquid phase	OK
double specvolla(double Tr,double Pr) {
	vLa=R*TB/PB*(A1a+A3a*Tr+A4a*Tr*Tr+A2a*Pr);
	vLa=vLa/100;
	return vLa;
}

//enthalpy of pure water in gaseous phase	OK
double enthalpygw(double Tr,double Pr) {
	hgw=-R*TB*(-hgrow+D1w*(Trow-Tr)+D2w/2*(Trow*Trow-Tr*Tr)+
		D3w/3*(pow(Trow,3)-pow(Tr,3))-C1w*(Pr-Prow)-
		4*C2w*(Pr*pow(Tr,-3)-Prow*pow(Trow,-3))-
		12*C3w*(Pr*pow(Tr,-11)-Prow*pow(Trow,-11))-
		4*C4w*(pow(Pr,3)*pow(Tr,-11)-pow(Prow,3)*pow(Trow,-11)));
	return hgw;
}

//enthalpy of pure ammonia in gaseous phase		OK
double enthalpyga(double Tr,double Pr) {
	hga=-R*TB*(-hgroa+D1a*(Troa-Tr)+D2a/2*(Troa*Troa-Tr*Tr)+
		D3a/3*(pow(Troa,3)-pow(Tr,3))-C1a*(Pr-Proa)-
		4*C2a*(Pr*pow(Tr,-3)-Proa*pow(Troa,-3))-
		12*C3a*(Pr*pow(Tr,-11)-Proa*pow(Troa,-11))-
		4*C4a*(pow(Pr,3)*pow(Tr,-11)-pow(Proa,3)*pow(Troa,-11)));
	return hga;
}

//entropy of pure water in gaseous phase	OK
double entropygw(double Tr,double Pr) {
	if(Pr==0)
		sgw=0;
	else
		sgw=-R*(-sgrow-D1w*log(Tr/Trow)+D2w*(Trow-Tr)+D3w/2*(Trow*Trow-Tr*Tr)+
			log(Pr/Prow)-3*C2w*(Pr*pow(Tr,-4)-Prow*pow(Trow,-4))-11*C3w*
			(Pr*pow(Tr,-12)-Prow*pow(Trow,-12))-11/3*C4w*(pow(Pr,3)*
			pow(Tr,-12)-pow(Prow,3)*pow(Trow,-12)));
	return sgw;
}

//entropy of pure ammonia in gaseous phase	OK
double entropyga(double Tr,double Pr) {
	if(Pr==0)
		sga=0;
	else
		sga=-R*(-sgroa-D1a*log(Tr/Troa)+D2a*(Troa-Tr)+D3a/2*(Troa*Troa-Tr*Tr)+
			log(Pr/Proa)-3*C2a*(Pr*pow(Tr,-4)-Proa*pow(Troa,-4))-11*C3a*
			(Pr*pow(Tr,-12)-Proa*pow(Troa,-12))-11/3*C4a*(pow(Pr,3)*
			pow(Tr,-12)-pow(Proa,3)*pow(Troa,-12)));
	return sga;
}

//specific volume of pure water in gaseous phase	OK
double specvolgw(double Tr,double Pr) {
	if(Pr==0)
		vgw=0;
	else
		vgw=R*TB/PB*(Tr/Pr+C1w+C2w*pow(Tr,-3)+C3w*pow(Tr,-11)+
			C4w*pow(Pr,2)*pow(Tr,-11));
	vgw=vgw/100;
	return vgw;
}

//specific volume of pure ammonia in gaseous phase	OK
double specvolga(double Tr,double Pr) {
	if(Pr==0)
		vga=0;
	else
		vga=R*TB/PB*(Tr/Pr+C1a+C2a*pow(Tr,-3)+C3a*pow(Tr,-11)+
			C4a*pow(Pr,2)*pow(Tr,-11));
	vga=vga/100;
	return vga;
}

//excess enthalpy of the mixture (in liquid phase)	OK
double enthalpyE(double Tr,double Pr,double y) {
	hE=-R*TB*(1-y)*y*(-E1-E2*Pr-2*E5/Tr-3*E6*pow(Tr,-2)+
		(2*y-1)*(-E7-E8*Pr-2*E11/Tr-3*E12*pow(Tr,-2))+
		pow((2*y-1),2)*(-E13-E14*Pr-2*E15/Tr-3*E16*pow(Tr,-2)));
	return hE;		// kJ/kmol(mix)
}

//excess entropy of the mixture (in liquid phase)	OK
double entropyE(double Tr,double Pr,double y) {
	sE=-R*(1-y)*y*(E3+E4*Pr-E5*pow(Tr,-2)-2*E6*pow(Tr,-3)+
		(2*y-1)*(E9+E10*Pr-E11*pow(Tr,-2)-2*E12*pow(Tr,-3))+
		pow((2*y-1),2)*(-E15*pow(Tr,-2)-2*E16*pow(Tr,-3)));
	return sE;
}

//excess specific volume of the mixture (in liquid phase)	OK
double specvolE(double Tr,double Pr,double y) {
	vE=R*TB/PB*(1-y)*y*(E2+E4*Tr+(2*y-1)*(E8+E10*Tr)+pow((2*y-1),2)*E14);
	vE=vE/100;
	return vE;
}

//entropy of mixing		OK
double entropymixing(double y) {
	if(y==0||y==1)
		smix=0;
	else
		smix=-R*(y*log(y)+(1-y)*log(1-y));
	return smix;
}

//enthalpy of the mixture in liquid phase	OK
double enthalpyL(double hLa,double hLw,double hE,double y) {
	hLm=y*hLa+(1-y)*hLw+hE;		// kJ/kmol(mix)
	return hLm;
}

//entropy of the mixture in liquid phase	OK
double entropyL(double sLa,double sLw,double sE,double smix,double y) {
	sLm=y*sLa+(1-y)*sLw+sE+smix;
	return sLm;
}

//specfic volume of the mixture in liquid phase		OK
double specvolL(double vLa,double vLw,double vE,double y) {
	vLm=y*vLa+(1-y)*vLw+vE;
	return vLm;
}

//enthalpy of mixture in gaseous phase
double enthalpyg(double hga,double hgw,double y) {
	hgm=y*hga+(1-y)*hgw;
	return hgm;
}

//entropy of mixture in gaseous phase
double entropyg(double sga,double sgw,double smix,double y) {
	sgm=y*sga+(1-y)*sgw+smix;
	return sgm;
}

//specific volume of mixture in gaseous phase
double specvolg(double vga,double vgw,double y) {
	vgm=y*vga+(1-y)*vgw;
	return vgm;
}

//bubble point temperature of the mixture		OK
double bubbleT(double PP,double xx) {			//P in bar
	Ci[1] = 153.634521459;
	Ci[2] = -13.0305543892;
	Ci[3] = -1.14845282991;
	Ci[4] = .550358094447;
	Ci[5] = -.0753450148427;
	Ci[6] = .0048111666267;
	Ci[7] = -.000120433757177;
	Cij[1][1] = -462.460321366;
	Cij[1][2] = 23739.9986309;
	Cij[1][3] = -194504.35292;
	Cij[1][4] = 639383.528867;
	Cij[1][5] = -523748.057636;
	Cij[1][6] = -2328271.47551;
	Cij[1][7] = 7562418.53499;
	Cij[1][8] = -9668295.89504;
	Cij[1][9] = 5922081.87086;
	Cij[1][10] = -1432405.52125;
	Cij[2][1] = 421.443122208;
	Cij[2][2] = -14560.354925;
	Cij[2][3] = 53051.4495633;
	Cij[2][4] = 382763.793582;
	Cij[2][5] = -3583589.86875;
	Cij[2][6] = 12243265.3815;
	Cij[2][7] = -22307970.0156;
	Cij[2][8] = 22896656.8499;
	Cij[2][9] = -12483324.8091;
	Cij[2][10] = 2813311.71633;
	Cij[3][1] = -248.783804168;
	Cij[3][2] = 4807.07241098;
	Cij[3][3] = 13565.1003309;
	Cij[3][4] = -466407.780832;
	Cij[3][5] = 2827083.44764;
	Cij[3][6] = -8469715.15799;
	Cij[3][7] = 14459588.8962;
	Cij[3][8] = -14281087.5331;
	Cij[3][9] = 7596403.59678;
	Cij[3][10] = -1684002.64482;
	Cij[4][1] = 126.965580728;
	Cij[4][2] = -2090.45270574;
	Cij[4][3] = 1993.17101166;
	Cij[4][4] = 100706.510396;
	Cij[4][5] = -687388.808612;
	Cij[4][6] = 2132412.46959;
	Cij[4][7] = -3699199.65914;
	Cij[4][8] = 3688365.22546;
	Cij[4][9] = -1975122.39296;
	Cij[4][10] = 440201.446068;
	Cij[5][1] = -33.5343446156;
	Cij[5][2] = 601.878586689;
	Cij[5][3] = -3064.82070658;
	Cij[5][4] = 71.7954752052;
	Cij[5][5] = 51780.666659;
	Cij[5][6] = -209714.899856;
	Cij[5][7] = 405011.985355;
	Cij[5][8] = -428310.461566;
	Cij[5][9] = 238153.698326;
	Cij[5][10] = -54497.0973336;
	Cij[6][1] = 3.97454953787;
	Cij[6][2] = -77.026846469;
	Cij[6][3] = 541.19105807;
	Cij[6][4] = -1696.60270972;
	Cij[6][5] = 1713.45942707;
	Cij[6][6] = 4019.01019872;
	Cij[6][7] = -14844.7928004;
	Cij[6][8] = 19481.0094551;
	Cij[6][9] = -12107.0794501;
	Cij[6][10] = 2966.92804386;
	Cij[7][1] = -.170806170177;
	Cij[7][2] = 3.48182859299;
	Cij[7][3] = -27.7957587743;
	Cij[7][4] = 113.762064546;
	Cij[7][5] = -258.750496922;
	Cij[7][6] = 311.002585218;
	Cij[7][7] = -123.917993454;
	Cij[7][8] = -123.480627492;
	Cij[7][9] = 154.375042114;
	Cij[7][10] = -48.5083828701;
	Tc=crittemp(xx);
	Pc=critpress(xx);
	i=1;
	sum2=0;
	while (i<=7) {
		j=1;
		sum1=0;
		while (j<=10) {
			sum1=sum1+Cij[i][j]*pow(xx,j);
			j++;
		}
		sum2=sum2+(Ci[i]+sum1)*pow(log(Pc/PP),i);
		i++;
	}
	Tb=Tc-sum2/1.8;							//gives temp in K
	return Tb;
}

//dew point temperature of the mixture		OK
double dewT(double PP,double xx) {		//P in bar
	ai[1] = 153.17055346;
	ai[2] = -11.7705687461;
	ai[3] = -1.78126355957;
	ai[4] = .647385455059;
	ai[5] = -.0719950751898;
	ai[6] = .00285423950786;
	Aij[1][1] = 194.793913463;
	Aij[1][2] = 74.236124188;
	Aij[1][3] = 9.84103819552;
	Aij[1][4] = .436843852745;
	Aij[2][1] = -74.3508283362;
	Aij[2][2] = -33.2941879809;
	Aij[2][3] = -4.78866918581;
	Aij[2][4] = -.225416733476;
	Aij[3][1] = 13.0175447367;
	Aij[3][2] = 6.1586564117;
	Aij[3][3] = .789740337141;
	Aij[3][4] = .0321510834958;
	Aij[4][1] = -.90857587517;
	Aij[4][2] = -.356752691147;
	Aij[4][3] = .0238067275502;
	Aij[4][4] = .00495593933952;
	Aij[5][1] = -.00071863574153;
	Aij[5][2] = -.0251026383533;
	Aij[5][3] = -.0191664613304;
	Aij[5][4] = -.0017014253867;
	Aij[6][1] = .00195441702983;
	Aij[6][2] = .00280533348937;
	Aij[6][3] = .0013899436563;
	Aij[6][4] = .000116422611616;
	Tc=crittemp(xx);
	Pc=critpress(xx);
	i=1;
	sum2=0;
	while (i<=6) {
		j=1;
		sum1=0;
		while (j<=4) {
			sum1=sum1+Aij[i][j]*pow(log(1.0001-xx),j);
			j++;
		}
		sum2=sum2+(ai[i]+sum1)*pow(log(Pc/PP),i);
		i++;
	}
	Td=Tc-sum2/1.8;					//gives temp in K
	return Td;
}

//bubble point pressure of the mixture		OK
double bubbleP(double TT,double xx) {
	Pb=critpress(xx);
	Tb=1000;
	incr=10;
	while(fabs(TT-Tb)>.0001) {
		while(TT-Tb<0) {
			Pb=Pb-incr;
			while(Pb<=incr) {
				Pb=Pb+incr;
				incr=incr/10;
			}
			Tb=bubbleT(Pb,xx);
		}
		Pb=Pb+incr;
		Tb=bubbleT(Pb,xx);
		incr=incr/10;
	}
	return Pb;
}

//dew point pressure of the mixture		OK
double dewP(double TT,double xx) {
	Pd=critpress(xx);
	Td=1000;
	incr=10;
	while(fabs(TT-Td)>.0001) {
		while(TT-Td<0) {
			Pd=Pd-incr;
			while(Pd<=incr) {
				Pd=Pd+incr;
				incr=incr/10;
			}
			Td=dewT(Pd,xx);
		}
		Pd=Pd+incr;
		Td=dewT(Pd,xx);
		incr=incr/10;
	}
	return Pd;
}

//bubble point ammonia mass fraction of the mixture		OK
double bubblex(double TT,double PP) {
	xb=0;
	Tc=crittemp(xb);
	Pc=critpress(xb);
	Tb=10000;
	incr=.1;
	while(fabs(TT-Tb)>.0001) {
		while(TT-Tb<0) {
			xb=xb+incr;
			if(xb>1) {
				xb=2;
				goto done;
			}
			Tc=crittemp(xb);
			Pc=critpress(xb);
			Tb=bubbleT(PP,xb);
		}
		xb=xb-incr;
		if(xb<0) {
			xb=-1;
			goto done;
		}
		Tc=crittemp(xb);
		Pc=critpress(xb);
		Tb=bubbleT(PP,xb);
		incr=incr/10;
	}
done:
	return xb;
}

//dew point ammonia mass fraction of the mixture	OK
double dewx(double TT,double PP) {
	xd=0;
	Tc=crittemp(xd);
	Pc=critpress(xd);
	Td=10000;
	incr=.1;
	while(fabs(TT-Td)>.0001) {
		while(TT-Td<0) {
			xd=xd+incr;
			if(xd>1) {
				xd=2;
				goto done;
			}
			Tc=crittemp(xd);
			Pc=critpress(xd);
			Td=dewT(PP,xd);
		}
		xd=xd-incr;
		if(xd<0) {
			xd=-1;
			goto done;
		}
		Tc=crittemp(xd);
		Pc=critpress(xd);
		Td=dewT(PP,xd);
		incr=incr/10;
	}
done:
	return xd;
}

//critical temperature of the mixture		OK
double crittemp(double xx) {
	a[1] = 205.8889;
	a[2] = 280.930556;
	a[3] = -317.0138889;
	a[4] = 263.194444;
	Tc=0;
	i=1;
	sum1=0;
	while (i<=4) {
		sum1=sum1+(a[i]*pow(xx,i));
		i++;
	}
	Tc=(Tcw-sum1)/1.8;			//gives temperature in K
	return Tc;
}

//critical pressure of the mixture		OK
double critpress(double xx) {
	b[1] = .368105523897;
	b[2] = -3.6679548875;
	b[3] = 46.6000470809;
	b[4] = -262.921061996;
	b[5] = 732.99536936;
	b[6] = -1076.0613489;
	b[7] = 797.948078048;
	b[8] = -235.903904222;
	Pc=0;
	i=1;
	sum1=0;
	while (i<=8) {
		sum1=sum1+(b[i]*pow(xx,i));
		i++;
	}
	Pc=Pcw*exp(sum1)/14.51135;	//gives pressure in bar
	return Pc;
}

void PTx (double PP,double TT,double xx) {
	Tr=TT/TB;
	Pr=PP/PB;
	M=18.015*17.031/((1-xx)*17.031+xx*18.015);
	y=xx*18.015/(xx*18.015+(1-xx)*17.031);
	Tb=bubbleT(PP,xx);
	Td=dewT(PP,xx);
	if(TT<Tb) {
		istate=1;
		hm=enthalpyL(enthalpyla(Tr,Pr),enthalpylw(Tr,Pr),
			enthalpyE(Tr,Pr,y),y)/M;
		sm=entropyL(entropyla(Tr,Pr),entropylw(Tr,Pr),
			entropyE(Tr,Pr,y),entropymixing(y),y)/M;
		vm=specvolL(specvolla(Tr,Pr),specvollw(Tr,Pr),
			specvolE(Tr,Pr,y),y)/M;
		xxNH3v=xxH2Ov=0;
		xxNH3l=xx;
		xxH2Ol=1-xx;
	}
	else if(TT>Td) {
		istate=3;
		hm=enthalpyg(enthalpyga(Tr,Pr),enthalpygw(Tr,Pr),y)/M;
		sm=entropyg(entropyga(Tr,Pr),entropygw(Tr,Pr),
			entropymixing(y),y)/M;
		vm=specvolg(specvolga(Tr,Pr),specvolgw(Tr,Pr),y)/M;
		xxNH3l=xxH2Ol=0;
		xxNH3v=xx;
		xxH2Ov=1-xx;
	}
	else {
		istate=2;
		xb=bubblex(TT,PP);
		xd=dewx(TT,PP);
		qm=(xx-xb)/(xd-xb);			/* quality of mixture */
		yb=xb*18.015/(xb*18.015+(1-xb)*17.031);
		yd=xd*18.015/(xd*18.015+(1-xd)*17.031);
		Mb=18.015*17.031/((1-xb)*17.031+xb*18.015);
		Md=18.015*17.031/((1-xd)*17.031+xd*18.015);
		hm=(1-qm)/Mb*enthalpyL(enthalpyla(Tr,Pr),enthalpylw(Tr,Pr),
			enthalpyE(Tr,Pr,yb),yb)+qm/Md*enthalpyg(enthalpyga
			(Tr,Pr),enthalpygw(Tr,Pr),yd);  //kJ/kg(mix at current state)
		sm=(1-qm)/Mb*entropyL(entropyla(Tr,Pr),entropylw(Tr,Pr),
			entropyE(Tr,Pr,yb),entropymixing(yb),yb)+qm/Md*
			entropyg(entropyga(Tr,Pr),entropygw(Tr,Pr),
			entropymixing(yd),yd);
		vm=(1-qm)/Mb*specvolL(specvolla(Tr,Pr),specvollw(Tr,Pr),
			specvolE(Tr,Pr,yb),yb)+qm/Md*specvolg(specvolga(Tr,Pr),
			specvolgw(Tr,Pr),yd);
		xxNH3v=(xx-xb)/(xd-xb)*xd;  //mass NH3 vap per mass total
		xxH2Ov=(xx-xb)/(xd-xb)*(1-xd);
		xxNH3l=(1-(xx-xb)/(xd-xb))*xb;
		xxH2Ol=(1-(xx-xb)/(xd-xb))*(1-xb);
	}
}

//gives pure component properties
void purePTx (double PP,double TT,double xx) {
	Tr=TT/TB;
	Pr=PP/PB;
	M=18.015*17.031/((1-xx)*17.031+xx*18.015);
	y=xx*18.015/(xx*18.015+(1-xx)*17.031);
	Tb=bubbleT(PP,xx);
	Td=dewT(PP,xx);
	Tb=.5*(Tb+Td);
	if(TT<Tb) {
		istate=1;
		hm=enthalpyL(enthalpyla(Tr,Pr),enthalpylw(Tr,Pr),
			enthalpyE(Tr,Pr,y),y)/M;
		sm=entropyL(entropyla(Tr,Pr),entropylw(Tr,Pr),
			entropyE(Tr,Pr,y),entropymixing(y),y)/M;
		vm=specvolL(specvolla(Tr,Pr),specvollw(Tr,Pr),
			specvolE(Tr,Pr,y),y)/M;
		xxNH3v=xxH2Ov=0;
		xxNH3l=xx;
		xxH2Ol=1-xx;
	}
	else if(TT>Tb){
		istate=3;
		hm=enthalpyg(enthalpyga(Tr,Pr),enthalpygw(Tr,Pr),y)/M;
		sm=entropyg(entropyga(Tr,Pr),entropygw(Tr,Pr),
			entropymixing(y),y)/M;
		vm=specvolg(specvolga(Tr,Pr),specvolgw(Tr,Pr),y)/M;
		xxNH3l=xxH2Ol=0;
		xxNH3v=xx;
		xxH2Ov=1-xx;
	}
	else {
		istate=2;
		printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
		printf("Enter quality of mixture --> ");
		scanf("%lf",Q[1]);
		y=xx*18.015/(xx*18.015+(1-xx)*17.031);
		M=18.015*17.031/((1-xx)*17.031+xx*18.015);
		hm=(1-Q[1])/M*enthalpyL(enthalpyla(Tr,Pr),enthalpylw(Tr,Pr),
			enthalpyE(Tr,Pr,y),y)+Q[1]/M*enthalpyg(enthalpyga
			(Tr,Pr),enthalpygw(Tr,Pr),y);
		sm=(1-Q[1])/M*entropyL(entropyla(Tr,Pr),entropylw(Tr,Pr),
			entropyE(Tr,Pr,y),entropymixing(y),y)+Q[1]/M*
			entropyg(entropyga(Tr,Pr),entropygw(Tr,Pr),
			entropymixing(y),y);
		vm=(1-Q[1])/M*specvolL(specvolla(Tr,Pr),specvollw(Tr,Pr),
			specvolE(Tr,Pr,y),y)+Q[1]/M*specvolg(specvolga(Tr,Pr),
			specvolgw(Tr,Pr),y);
		if(xx<.0001) {
			xxH2Ov=Q[1];xxH2Ol=1-Q[1];
		}
		else {
			xxNH3v=Q[1];xxNH3l=1-Q[1];
		}
	}
}

//gives sat properties for T,x
void Tx (double TT,double xx) {
	Tr=TT/TB;
	y=xx*18.015/(xx*18.015+(1-xx)*17.031);
	Pb=bubbleP(TT,xx);
	Pd=dewP(TT,xx);
	M=18.015*17.031/((1-xx)*17.031+xx*18.015);
//	if(TT>dewT(Pb,xx)&&xx<.1)
//		xd=xx;
//	else if(TT<dewT(Pb,1))
//		xd=1;
//	else					//xx is the sat liq mass fraction
//		xd=dewx(TT,Pb);		//xd is the sat vap mass fraction
//	Md=18.015*17.031/((1-xd)*17.031+xd*18.015);
//	yd=xd*18.015/(xd*18.015+(1-xd)*17.031);
	PP=Pb;
// sat liq h
	hm=hb=enthalpyL(enthalpyla(Tr,Pb/PB),enthalpylw(Tr,Pb/PB),enthalpyE(Tr,Pb/PB,y),y)/M;
// sat vap h
	hd=enthalpyg(enthalpyga(Tr,Pd/PB),enthalpygw(Tr,Pd/PB),y)/M;
// sat liq s 
	sm=sb=entropyL(entropyla(Tr,Pb/PB),entropylw(Tr,Pb/PB),entropyE(Tr,Pb/PB,y),entropymixing(y),y)/M;
// sat vap s
	sd=entropyg(entropyga(Tr,Pd/PB),entropygw(Tr,Pd/PB),entropymixing(y),y)/M;
// sat liq v
	vm=vb=specvolL(specvolla(Tr,Pb/PB),specvollw(Tr,Pb/PB),specvolE(Tr,Pb/PB,y),y)/M;
// sat vap v
	vd=specvolg(specvolga(Tr,Pd/PB),specvolgw(Tr,Pd/PB),y)/M;
}

//gives sat properties for P,x
void Px (double PP,double xx) {
	Pr=PP/PB;
	y=xx*18.015/(xx*18.015+(1-xx)*17.031);
	Tb=bubbleT(PP,xx);
	Td=dewT(PP,xx);
	M=18.015*17.031/((1-xx)*17.031+xx*18.015);
//	if(PP<dewP(Tb,xx)&&xx<.1)
//		xd=xx;
//	else if(PP>dewP(Tb,1))
//		xd=1;
//	else					//xx is the sat liq mass fraction
//		xd=dewx(Tb,PP);		//xd is the sat vap mass fraction
//	Md=18.015*17.031/((1-xd)*17.031+xd*18.015);
//	yd=xd*18.015/(xd*18.015+(1-xd)*17.031);
	TT=Tb;
// sat liq h
	hm=hb=enthalpyL(enthalpyla(Tb/TB,Pr),enthalpylw(Tb/TB,Pr),enthalpyE(Tb/TB,Pr,y),y)/M;
// sat vap h
	hd=enthalpyg(enthalpyga(Td/TB,Pr),enthalpygw(Td/TB,Pr),y)/M;
// sat liq s 
	sm=sb=entropyL(entropyla(Tb/TB,Pr),entropylw(Tb/TB,Pr),entropyE(Tb/TB,Pr,y),entropymixing(y),y)/M;
// sat vap s
	sd=entropyg(entropyga(Td/TB,Pr),entropygw(Td/TB,Pr),entropymixing(y),y)/M;
// sat liq v
	vm=vb=specvolL(specvolla(Tb/TB,Pr),specvollw(Tb/TB,Pr),specvolE(Tb/TB,Pr,y),y)/M;
// sat vap v
	vd=specvolg(specvolga(Td/TB,Pr),specvolgw(Td/TB,Pr),y)/M;
}


//gives sat liq properties for P,T
void PT (double PP,double TT) {
	Tr=TT/TB;
	Pr=PP/PB;
	xb=bubblex(TT,PP);
	if(xb<0) {
		hm=hb=hd=sm=sb=sd=vm=vb=vd=xb=0;
		goto done;
	}
	xx=xb;
	y=xx*18.015/(xx*18.015+(1-xx)*17.031);
	M=18.015*17.031/((1-xx)*17.031+xx*18.015);
	xd=dewx(TT,PP);		//xd is the sat vap mass fraction
	if(xd>1) {
		hm=hb=hd=sm=sb=sd=vm=vb=vd=0;
		xd=1;
		goto done;
	}
	Md=18.015*17.031/((1-xd)*17.031+xd*18.015);
	yd=xd*18.015/(xd*18.015+(1-xd)*17.031);
// sat liq h
	hm=hb=enthalpyL(enthalpyla(Tr,Pr),enthalpylw(Tr,Pr),enthalpyE(Tr,Pr,y),y)/M;
// sat vap h
	hd=enthalpyg(enthalpyga(Tr,Pr),enthalpygw(Tr,Pr),yd)/Md;
// sat liq s 
	sm=sb=entropyL(entropyla(Tr,Pr),entropylw(Tr,Pr),entropyE(Tr,Pr,y),entropymixing(y),y)/M;
// sat vap s
	sd=entropyg(entropyga(Tr,Pr),entropygw(Tr,Pr),entropymixing(yd),yd)/Md;
// sat liq v
	vm=vb=specvolL(specvolla(Tr,Pr),specvollw(Tr,Pr),specvolE(Tr,Pr,y),y)/M;
// sat vap v 
	vd=specvolg(specvolga(Tr,Pr),specvolgw(Tr,Pr),yd)/Md;
done:
	;
}

//gives mixture properties, given P,s,x
void Psx (double PP,double temps,double xx) {
	TT=400;
	incr2=10;
	PTx(PP,TT,xx);
	while (fabs(temps-sm)>.001) {
		while(sm>temps) {
			TT=TT-incr2;
			PTx(PP,TT,xx);
		}
		incr2=incr2/10;
		while(sm<temps) {
			TT=TT+incr2;
			PTx(PP,TT,xx);
		}
		incr2=incr2/10;
	}
}

//gives mixture properties, given T,v,x
void Tvx (double TT,double tempv,double xx) {
	PP=1;
	incr2=.01;
	limit=.000001;
	PTx(PP,TT,xx);
	while (fabs(tempv-vm)>limit) {
		while(vm<tempv) {
			PP=PP-incr2;
			PTx(PP,TT,xx);
		}
		incr2=incr2/10;
		n=1;
		while(vm>tempv) {
			if(n>=11) {
				incr2=incr2*10;
				n=1;
			}
			PP=PP+incr2;
			PTx(PP,TT,xx);
			n=n+1;
		}
		incr2=incr2/10;
	}
}

//gives mixture properties, given P,v,x
void Pvx (double PP,double tempv,double xx) {
	TT=300;
	incr2=1;
	limit=.000001;
	PTx(PP,TT,xx);
	while (fabs(tempv-vm)>limit) {
		n=1;
		while(vm>tempv) {
			if(n>=11) {
				incr2=incr2*10;
				n=1;
			}
			TT=TT-incr2;
			PTx(PP,TT,xx);
			n=n+1;
		}
		incr2=incr2/10;
		n=1;
		while(vm<tempv) {
			if(n>=11) {
				incr2=incr2*10;
				n=1;
			}
			TT=TT+incr2;
			PTx(PP,TT,xx);
			n=n+1;
		}
		incr2=incr2/10;
	}
}

//gives mixture properties, given P,T,v
void PTv (double PP,double TT,double tempv) {
	xx=.5;
	incr2=.01;
	limit=.000001;
	PTx(PP,TT,xx);
	while (fabs(tempv-vm)>limit) {
		n=1;
		while(vm>tempv) {
			if(n>=11) {
				incr2=incr2*10;
				n=1;
			}
			xx=xx-incr2;
			if(xx<0) {
				printf("mixture not possible\n");
				getch();
				goto done;
			}
			PTx(PP,TT,xx);
			n=n+1;
		}
		incr2=incr2/10;
		n=1;
		while(vm<tempv) {
			if(n>=11) {
				incr2=incr2*10;
				n=1;
			}
			xx=xx+incr2;
			if(xx>1) {
				printf("mixture not possible\n");
				getch();
				goto done;
			}
			PTx(PP,TT,xx);
			n=n+1;
		}
		incr2=incr2/10;
	}
done:
	;
}


//gives pure properties, given T,v,x
void pureTvx (double TT,double tempv,double xx) {
	Tr=TT/TB;
	M=18.015*17.031/((1-xx)*17.031+xx*18.015);
	y=xx*18.015/(xx*18.015+(1-xx)*17.031);
	Tx(TT,xx);
	if(tempv<vb) {
		limit=.000000001;
		PP=1;
		incr2=.1;
		while (fabs(tempv-vm)>limit) {
			while(vm<tempv&&PP>=0) {
				PP=PP-incr2;
				purePTx(PP,TT,xx);
			}
			incr2=incr2/10;
			n=1;
			while(vm>tempv) {
				if(n>=11) {
					incr2=incr2*10;
					n=1;
				}
				PP=PP+incr2;
				purePTx(PP,TT,xx);
				n=n+1;
			}
			incr2=incr2/10;
		}
		purePTx(PP,TT,xx);
	}
	else if(tempv>vd){
		limit=.0000001;
		PP=1;
		incr2=.1;
		while (fabs(tempv-vm)>limit) {
			while(vm<tempv&&PP>=0) {
				PP=PP-incr2;
				purePTx(PP,TT,xx);
			}
			incr2=incr2/10;
			n=1;
			while(vm>tempv) {
				if(n>=11) {
					incr2=incr2*10;
					n=1;
				}
				PP=PP+incr2;
				purePTx(PP,TT,xx);
				n=n+1;
			}
			incr2=incr2/10;
		}
		purePTx(PP,TT,xx);
	}
	else {
		Tx(TT,xx);
		PP=.5*(Pb+Pd);
		istate=2;
		Q[1]=(tempv-vb)/(vd-vb);	//quality of pure fluid
		hm=(1-Q[1])/M*enthalpyL(enthalpyla(Tr,Pr),enthalpylw(Tr,Pr),
			enthalpyE(Tr,Pr,y),y)+Q[1]/M*enthalpyg(enthalpyga
			(Tr,Pr),enthalpygw(Tr,Pr),y);
		sm=(1-Q[1])/M*entropyL(entropyla(Tr,Pr),entropylw(Tr,Pr),
			entropyE(Tr,Pr,y),entropymixing(y),y)+Q[1]/M*
			entropyg(entropyga(Tr,Pr),entropygw(Tr,Pr),
			entropymixing(y),y);
		vm=(1-Q[1])/M*specvolL(specvolla(Tr,Pr),specvollw(Tr,Pr),
			specvolE(Tr,Pr,y),y)+Q[1]/M*specvolg(specvolga(Tr,Pr),
			specvolgw(Tr,Pr),y);
		if(xx<0.0001) {
			xxH2Ov=Q[1];xxH2Ol=1-Q[1];
		}
		else {
			xxNH3v=Q[1];xxNH3l=1-Q[1];
		}
	}
}

void hPx (double temph,double PP,double xx) {
	TT=400;
	incr2=10;
	PTx(PP,TT,xx);
	while (fabs(temph-hm)>.01) {
		while(hm>temph) {
			TT=TT-incr2;
			PTx(PP,TT,xx);
		}
		incr2=incr2/10;
		while(hm<temph) {
			TT=TT+incr2;
			PTx(PP,TT,xx);
		}
		incr2=incr2/10;
	}
}
