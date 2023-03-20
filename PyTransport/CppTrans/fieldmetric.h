//#This file is part of PyTransport.

//#PyTransport is free software: you can redistribute it and/or modify
//#it under the terms of the GNU General Public License as published by
//#the Free Software Foundation, either version 3 of the License, or
//#(at your option) any later version.

//#PyTransport is distributed in the hope that it will be useful,
//#but WITHOUT ANY WARRANTY; without even the implied warranty of
//#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//#GNU General Public License for more details.

//#You should have received a copy of the GNU General Public License
//#along with PyTransport.  If not, see <http://www.gnu.org/licenses/>.

// This file contains a prototype of the potential.h file of PyTransport -- it is edited by the PyTransScripts module

#ifndef FIELDMETRIC_H  // Prevents the class being re-defined
#define FIELDMETRIC_H


#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;


class fieldmetric
{
private:
	int nF; // field number
    int nP; // params number which definFs potential

public:
    fieldmetric()
   {
// #FP
nF=2;
nP=3;
	
   }
	
	
	//calculates fieldmetic()
	vector<double> fmetric(vector<double> f, vector<double> p)
	{
		vector<double> FM((2*nF)*(2*nF),0.0) ;
        
// metric
  double x0 = std::pow(-std::pow(f[0], 2) - std::pow(f[1], 2) + 1, 2);
  double x1 = (1.0/6.0)*x0/p[0];
  double x2 = 6*p[0]/x0;

 FM[0]=x1;

 FM[1]=0;

 FM[2]=1;

 FM[3]=0;

 FM[4]=0;

 FM[5]=x1;

 FM[6]=0;

 FM[7]=1;

 FM[8]=1;

 FM[9]=0;

 FM[10]=x2;

 FM[11]=0;

 FM[12]=0;

 FM[13]=1;

 FM[14]=0;

 FM[15]=x2;

         return FM;
	}
	
	
	
	//calculates ChristoffelSymbole()
	vector<double> Chroff(vector<double> f, vector<double> p)
	{
		vector<double> CS((2*nF)*(2*nF)*(2*nF),0.0);
	
// Christoffel
  double x0 = std::pow(f[0], 2) + std::pow(f[1], 2) - 1;
  double x1 = 2/x0;
  double x2 = f[0]*x1;
  double x3 = -x2;
  double x4 = f[1]*x1;
  double x5 = -x4;

 CS[0]=0;

 CS[1]=0;

 CS[2]=0;

 CS[3]=0;

 CS[4]=0;

 CS[5]=0;

 CS[6]=0;

 CS[7]=0;

 CS[8]=0;

 CS[9]=0;

 CS[10]=x3;

 CS[11]=x5;

 CS[12]=0;

 CS[13]=0;

 CS[14]=x5;

 CS[15]=x2;

 CS[16]=0;

 CS[17]=0;

 CS[18]=0;

 CS[19]=0;

 CS[20]=0;

 CS[21]=0;

 CS[22]=0;

 CS[23]=0;

 CS[24]=0;

 CS[25]=0;

 CS[26]=x4;

 CS[27]=x3;

 CS[28]=0;

 CS[29]=0;

 CS[30]=x3;

 CS[31]=x5;

 CS[32]=0;

 CS[33]=0;

 CS[34]=0;

 CS[35]=0;

 CS[36]=0;

 CS[37]=0;

 CS[38]=0;

 CS[39]=0;

 CS[40]=0;

 CS[41]=0;

 CS[42]=0;

 CS[43]=0;

 CS[44]=0;

 CS[45]=0;

 CS[46]=0;

 CS[47]=0;

 CS[48]=0;

 CS[49]=0;

 CS[50]=0;

 CS[51]=0;

 CS[52]=0;

 CS[53]=0;

 CS[54]=0;

 CS[55]=0;

 CS[56]=0;

 CS[57]=0;

 CS[58]=0;

 CS[59]=0;

 CS[60]=0;

 CS[61]=0;

 CS[62]=0;

 CS[63]=0;
        
		return CS;
	}
    

	
	// calculates RiemannTensor()
	vector<double> Riemn(vector<double> f, vector<double> p)
	{
		vector<double> RM((nF)*(nF)*(nF)*(nF),0.0);
		
// Riemann
  double x0 = std::pow(f[0], 2);
  double x1 = std::pow(f[1], 2);
  double x2 = x0 + x1 - 1;
  double x3 = -x2;
  double x4 = 24*p[0]*(x0*x3 + x1*x3 + std::pow(x2, 2))/std::pow(x2, 5);
  double x5 = -x4;

 RM[0]=0;

 RM[1]=0;

 RM[2]=0;

 RM[3]=0;

 RM[4]=0;

 RM[5]=x4;

 RM[6]=x5;

 RM[7]=0;

 RM[8]=0;

 RM[9]=x5;

 RM[10]=x4;

 RM[11]=0;

 RM[12]=0;

 RM[13]=0;

 RM[14]=0;

 RM[15]=0;
     
        return RM;
	}

	// calculates RiemannTensor() covariant derivatives
	vector<double> Riemncd(vector<double> f, vector<double> p)
	{
		vector<double> RMcd((nF)*(nF)*(nF)*(nF)*(nF),0.0);
		
// Riemanncd
  double x0 = std::pow(f[0], 2);
  double x1 = std::pow(f[1], 2);
  double x2 = x0 + x1 - 1;
  double x3 = -x2;
  double x4 = std::pow(x2, 2);
  double x5 = x0*x3 + x1*x3 + x4;
  double x6 = std::pow(x2, 3) - 4*std::pow(x3, 2)*x5 + 5*x4*x5;
  double x7 = 48*p[0]/std::pow(x2, 8);
  double x8 = -x6*x7;
  double x9 = f[0]*x8;
  double x10 = f[1]*x8;
  double x11 = x6*x7;
  double x12 = f[0]*x11;
  double x13 = f[1]*x11;

 RMcd[0]=0;

 RMcd[1]=0;

 RMcd[2]=0;

 RMcd[3]=0;

 RMcd[4]=0;

 RMcd[5]=0;

 RMcd[6]=0;

 RMcd[7]=0;

 RMcd[8]=0;

 RMcd[9]=0;

 RMcd[10]=x9;

 RMcd[11]=x10;

 RMcd[12]=x12;

 RMcd[13]=x13;

 RMcd[14]=0;

 RMcd[15]=0;

 RMcd[16]=0;

 RMcd[17]=0;

 RMcd[18]=x12;

 RMcd[19]=x13;

 RMcd[20]=x9;

 RMcd[21]=x10;

 RMcd[22]=0;

 RMcd[23]=0;

 RMcd[24]=0;

 RMcd[25]=0;

 RMcd[26]=0;

 RMcd[27]=0;

 RMcd[28]=0;

 RMcd[29]=0;

 RMcd[30]=0;

 RMcd[31]=0;
     
        return RMcd;
	}
    
    int getnF()
    {
        return nF;
    }
    


};
#endif

