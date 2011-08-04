#include "nh3proplib.h"

//constructor
Props::Props() {
        sum1=0;
        sum2=0;
         Tcw=1165.1;	//critical temperature of water (R)
         Pcw=3208;	//critical pressure of water (psi)
         Tca=729.9;	//critical temperature of ammonia (R)
         Pca=1646;	//critical pressure of ammonia (psi)
         Ptpa=.04;      	//triple point of ammonia (bar)
         Ttpa=197;	//triple point of ammonia (K)
         Ptpw=.006113;   //triple point of water (bar)
         Ttpw=273.16;	//triple point of water (K)
         TB=100;		//reference temperature (K)
         PB=10;		//reference pressure (bar)
         R=8.31451;	//gas constant (kJ/kmolK)
          A1a=3.971423e-2,A1w=2.748796e-2;
          A2a=-1.790557e-5,A2w=-1.016665e-5;
          A3a=-1.308905e-2,A3w=-4.452025e-3;
          A4a=3.752836e-3,A4w=8.389246e-4;
          B1a=1.634519e+1,B1w=1.214557e+1;
          B2a=-6.508119,B2w=-1.898065;
          B3a=1.448937,B3w=2.911966e-1;
          C1a=-1.049377e-2,C1w=2.136131e-2;
          C2a=-8.288224,C2w=-3.169291e+1;
          C3a=-6.647257e+2,C3w=-4.634611e+4;
          C4a=-3.045352e+3,C4w=0.0;
          D1a=3.673647,D1w=4.019170;
          D2a=9.989629e-2,D2w=-5.175550e-2;
          D3a=3.617622e-2,D3w=1.951939e-2;
          hLroa=4.878573,hLrow=21.821141;
          hgroa=26.468873,hgrow=60.965058;
          sLroa=1.644773,sLrow=5.733498;
          sgroa=8.339026,sgrow=13.453430;

         Troa=3.2252,Trow=5.0705;
         Proa=2.0000,Prow=3.0000;

          E1 = -41.733398;
          E2 = 0.02414;
          E3 = 6.702285;
          E4 = -0.011475;
          E5 = 63.608967;
          E6 = -62.490768;
          E7 = 1.761064;
          E8 = 0.008626;
          E9 = 0.387983;
          E10 = -0.004772;
          E11 = -4.648107;
          E12 = 0.836376;
          E13 = -3.553627;
          E14 = 0.000904;
          E15 = 24.361723;
          E16 = -20.736547;



        }



//enthalpy of pure water in liquid phase	OK
double Props::enthalpylw(double Tr,double Pr) {
	hLw=-R*TB*(-hLrow+B1w*(Trow-Tr)+B2w/2*(Trow*Trow-Tr*Tr)+B3w/3*
		(pow(Trow,3)-pow(Tr,3))-(A1w+A4w*Tr*Tr)*(Pr-Prow)-
		A2w/2*(Pr*Pr-Prow*Prow));
	return hLw;		// kJ/kmol(water)
}

//enthalpy of pure ammonia in liquid phase	OK
double Props::enthalpyla(double Tr,double Pr) {
	hLa=-R*TB*(-hLroa+B1a*(Troa-Tr)+B2a/2*(Troa*Troa-Tr*Tr)+B3a/3*
		(pow(Troa,3)-pow(Tr,3))-(A1a+A4a*Tr*Tr)*(Pr-Proa)-
		A2a/2*(Pr*Pr-Proa*Proa));
	return hLa;		// kJ/kmol(ammonia)
}

//entropy of pure water in liquid phase		OK
double Props::entropylw(double Tr,double Pr) {
	sLw=-R*(-sLrow-B1w*log(Tr/Trow)+B2w*(Trow-Tr)+
		B3w/2*(Trow*Trow-Tr*Tr)+(Pr-Prow)*(A3w+2*A4w*Tr));
	return sLw;		//kJ/kmol K (water)
}

//entropy of pure ammonia in liquid phase	OK
double Props::entropyla(double Tr,double Pr) {
	sLa=-R*(-sLroa-B1a*log(Tr/Troa)+B2a*(Troa-Tr)+
		B3a/2*(Troa*Troa-Tr*Tr)+(Pr-Proa)*(A3a+2*A4a*Tr));
	return sLa;
}

//specific volume of water in liquid phase	OK
double Props::specvollw(double Tr,double Pr) {
	vLw=R*TB/PB*(A1w+A3w*Tr+A4w*Tr*Tr+A2w*Pr);
	vLw=vLw/100;		//m3/kmol (water)
	return vLw;
}

//specific volume of ammonia in liquid phase	OK
double Props::specvolla(double Tr,double Pr) {
	vLa=R*TB/PB*(A1a+A3a*Tr+A4a*Tr*Tr+A2a*Pr);
	vLa=vLa/100;
	return vLa;
}

//enthalpy of pure water in gaseous phase	OK
double Props::enthalpygw(double Tr,double Pr) {
	hgw=-R*TB*(-hgrow+D1w*(Trow-Tr)+D2w/2*(Trow*Trow-Tr*Tr)+
		D3w/3*(pow(Trow,3)-pow(Tr,3))-C1w*(Pr-Prow)-
		4*C2w*(Pr*pow(Tr,-3)-Prow*pow(Trow,-3))-
		12*C3w*(Pr*pow(Tr,-11)-Prow*pow(Trow,-11))-
		4*C4w*(pow(Pr,3)*pow(Tr,-11)-pow(Prow,3)*pow(Trow,-11)));
	return hgw;
}

//enthalpy of pure ammonia in gaseous phase		OK
double Props::enthalpyga(double Tr,double Pr) {
	hga=-R*TB*(-hgroa+D1a*(Troa-Tr)+D2a/2*(Troa*Troa-Tr*Tr)+
		D3a/3*(pow(Troa,3)-pow(Tr,3))-C1a*(Pr-Proa)-
		4*C2a*(Pr*pow(Tr,-3)-Proa*pow(Troa,-3))-
		12*C3a*(Pr*pow(Tr,-11)-Proa*pow(Troa,-11))-
		4*C4a*(pow(Pr,3)*pow(Tr,-11)-pow(Proa,3)*pow(Troa,-11)));
	return hga;
}

//entropy of pure water in gaseous phase	OK
double Props::entropygw(double Tr,double Pr) {
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
double Props::entropyga(double Tr,double Pr) {
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
double Props::specvolgw(double Tr,double Pr) {
	if(Pr==0)
		vgw=0;
	else
		vgw=R*TB/PB*(Tr/Pr+C1w+C2w*pow(Tr,-3)+C3w*pow(Tr,-11)+
			C4w*pow(Pr,2)*pow(Tr,-11));
	vgw=vgw/100;
	return vgw;
}

//specific volume of pure ammonia in gaseous phase	OK
double Props::specvolga(double Tr,double Pr) {
	if(Pr==0)
		vga=0;
	else
		vga=R*TB/PB*(Tr/Pr+C1a+C2a*pow(Tr,-3)+C3a*pow(Tr,-11)+
			C4a*pow(Pr,2)*pow(Tr,-11));
	vga=vga/100;
	return vga;
}

//excess enthalpy of the mixture (in liquid phase)	OK
double Props::enthalpyE(double Tr,double Pr,double y) {
	hE=-R*TB*(1-y)*y*(-E1-E2*Pr-2*E5/Tr-3*E6*pow(Tr,-2)+
		(2*y-1)*(-E7-E8*Pr-2*E11/Tr-3*E12*pow(Tr,-2))+
		pow((2*y-1),2)*(-E13-E14*Pr-2*E15/Tr-3*E16*pow(Tr,-2)));
	return hE;		// kJ/kmol(mix)
}

//excess entropy of the mixture (in liquid phase)	OK
double Props::entropyE(double Tr,double Pr,double y) {
	sE=-R*(1-y)*y*(E3+E4*Pr-E5*pow(Tr,-2)-2*E6*pow(Tr,-3)+
		(2*y-1)*(E9+E10*Pr-E11*pow(Tr,-2)-2*E12*pow(Tr,-3))+
		pow((2*y-1),2)*(-E15*pow(Tr,-2)-2*E16*pow(Tr,-3)));
	return sE;
}

//excess specific volume of the mixture (in liquid phase)	OK
double Props::specvolE(double Tr,double Pr,double y) {
	vE=R*TB/PB*(1-y)*y*(E2+E4*Tr+(2*y-1)*(E8+E10*Tr)+pow((2*y-1),2)*E14);
	vE=vE/100;
	return vE;
}

//entropy of mixing		OK
double Props::entropymixing(double y) {
	if(y==0||y==1)
		smix=0;
	else
		smix=-R*(y*log(y)+(1-y)*log(1-y));
	return smix;
}

//enthalpy of the mixture in liquid phase	OK
double Props::enthalpyL(double hLa,double hLw,double hE,double y) {
	hLm=y*hLa+(1-y)*hLw+hE;		// kJ/kmol(mix)
	return hLm;
}

//entropy of the mixture in liquid phase	OK
double Props::entropyL(double sLa,double sLw,double sE,double smix,double y) {
	sLm=y*sLa+(1-y)*sLw+sE+smix;
	return sLm;
}

//specfic volume of the mixture in liquid phase		OK
double Props::specvolL(double vLa,double vLw,double vE,double y) {
	vLm=y*vLa+(1-y)*vLw+vE;
	return vLm;
}

//enthalpy of mixture in gaseous phase
double Props::enthalpyg(double hga,double hgw,double y) {
	hgm=y*hga+(1-y)*hgw;
	return hgm;
}

//entropy of mixture in gaseous phase
double Props::entropyg(double sga,double sgw,double smix,double y) {
	sgm=y*sga+(1-y)*sgw+smix;
	return sgm;
}

//specific volume of mixture in gaseous phase
double Props::specvolg(double vga,double vgw,double y) {
	vgm=y*vga+(1-y)*vgw;
	return vgm;
}

//bubble point temperature of the mixture		OK
double Props::bubbleT(double PP,double xx) {			//P in bar
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
double Props::dewT(double PP,double xx) {		//P in bar
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
double Props::bubbleP(double TT,double xx) {
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
double Props::dewP(double TT,double xx) {
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
double Props::bubblex(double TT,double PP) {
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
double Props::dewx(double TT,double PP) {
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
double Props::crittemp(double xx) {
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
double Props::critpress(double xx) {
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

void Props::PTx (double PP,double TT,double xx) {
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
		qm=(xx-xb)/(xd-xb);			// quality of mixture
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
void Props::purePTx (double PP,double TT,double xx) {
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
		//printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
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
void Props::Tx (double TT,double xx) {
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
void Props::Px (double PP,double xx) {
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
void Props::PT (double PP,double TT) {
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
void Props::Psx (double PP,double temps,double xx) {
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
void Props::Tvx (double TT,double tempv,double xx) {
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
void Props::Pvx (double PP,double tempv,double xx) {
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
void Props::PTv (double PP,double TT,double tempv) {
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
void Props::pureTvx (double TT,double tempv,double xx) {
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

void Props::hPx (double temph,double PP,double xx) {
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


double Props::geth() {
return(hm);
}



double Props::getv() {
return(vm);
}









