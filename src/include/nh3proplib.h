#include<cstdio>

//remove old MS-DOS consol io header
//#include<conio.h>

#include<cmath>
#include<iostream>
#include<fstream>

using namespace std;

class Props {

private:
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

long double sum1,sum2;
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
double Tcw;//critical temperature of water (R)
double Pcw;	//critical pressure of water (psi)
double Tca;	//critical temperature of ammonia (R)
double Pca;	//critical pressure of ammonia (psi)
double Ptpa;      	//triple point of ammonia (bar)
double Ttpa;	//triple point of ammonia (K)
double Ptpw;   //triple point of water (bar)
double Ttpw;	//triple point of water (K)
double Tr,Pr;		//reduced temperature and pressure
double TB;		//reference temperature (K)
double PB;		//reference pressure (bar)
double R;	//gas constant (kJ/kmolK)
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

//some constants for ammmonia and water
long double A1a,A1w;
long double A2a,A2w;
long double A3a,A3w;
long double A4a,A4w;
long double B1a,B1w;
long double B2a,B2w;
long double B3a,B3w;
long double C1a,C1w;
long double C2a,C2w;
long double C3a,C3w;
long double C4a,C4w;
long double D1a,D1w;
long double D2a,D2w;
long double D3a,D3w;
long double hLroa,hLrow;
long double hgroa,hgrow;
long double sLroa,sLrow;
long double sgroa,sgrow;

double Troa,Trow;
double Proa,Prow;

long double E1;
long double E2;
long double E3;
long double E4;
long double E5;
long double E6;
long double E7;
long double E8;
long double E9;
long double E10;
long double E11;
long double E12;
long double E13;
long double E14;
long double E15;
long double E16;


public:

//functions

Props();    //Constructor
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

void PTx (double PP,double TT,double xx);
void purePTx (double PP,double TT,double xx);
void Tx (double TT,double xx);
void Px (double PP,double xx);
void PT (double PP,double TT);
void Psx (double PP,double temps,double xx);
void Tvx (double TT,double tempv,double xx);
void Pvx (double PP,double tempv,double xx);
void PTv (double PP,double TT,double tempv);
void pureTvx (double TT,double tempv,double xx);
void hPx (double temph,double PP,double xx);

double geth();
double getv();



};


