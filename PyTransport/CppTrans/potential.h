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

#ifndef POTENTIAL_H  // Prevents the class being re-defined
#define POTENTIAL_H


#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;

// #Rewrite
// Potential file rewriten at Mon Mar 20 17:12:55 2023

class potential
{
private:
	int nF; // field number
	int nP; // params number which definFs potential
    
    
public:
	// flow constructor
	potential()
	{
// #FP
nF=2;
nP=3;

//        p.resize(nP);
        
// pdef

    }
	
    //void setP(vector<double> pin){
    //    p=pin;
    //}
	//calculates V()
	double V(vector<double> f, vector<double> p)
	{
		double sum ;
        
// Pot
  sum=(1.0/2.0)*p[0]*std::pow(p[2], 2)*(std::pow(f[0], 2) + std::pow(f[1], 2)*p[1]);
         return sum;
	}
	
	//calculates V'()
	vector<double> dV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF,0.0);
	
// dPot
  double x0 = p[0]*std::pow(p[2], 2);

 sum[0]=f[0]*x0;

 sum[1]=f[1]*p[1]*x0;
        
		return sum;
	}
    
	// calculates V''
	vector<double> dVV(vector<double> f, vector<double> p)
	{
		vector<double> sum(nF*nF,0.0);
		
// ddPot
  double x0 = p[0]*std::pow(p[2], 2);
  double x1 = std::pow(f[0], 2);
  double x2 = std::pow(f[1], 2);
  double x3 = x1 + x2 - 1;
  double x4 = 2/x3;
  double x5 = x0*x4;
  double x6 = x1*x5;
  double x7 = p[1]*x0;
  double x8 = x4*x7;
  double x9 = x2*x8;
  double x10 = f[0]*f[1];
  double x11 = x10*x5 + x10*x8;

 sum[0]=x0 + x6 - x9;

 sum[2]=x11;

 sum[1]=x11;

 sum[3]=-x6 + x7 + x9;
     
        return sum;
	}
    
	// calculates V'''
	vector<double> dVVV(vector<double> f, vector<double> p)
	{
        vector<double> sum(nF*nF*nF,0.0);
// dddPot
  double x0 = std::pow(f[0], 3);
  double x1 = std::pow(f[0], 4);
  double x2 = 3*x1;
  double x3 = std::pow(f[1], 4);
  double x4 = 3*x3;
  double x5 = std::pow(f[0], 2);
  double x6 = std::pow(f[1], 2);
  double x7 = x5*x6;
  double x8 = std::pow(f[0], 6) + std::pow(f[1], 6) + x2*x6 - x2 + x4*x5 - x4 + 3*x5 + 3*x6 - 6*x7 - 1;
  double x9 = 1.0/x8;
  double x10 = 2*x9;
  double x11 = f[0]*x10;
  double x12 = std::pow(f[0], 5);
  double x13 = 4*x9;
  double x14 = -4*f[0]*x6*x9 + x0*x13*x6 - 4*x0*x9 + x10*x12 + x11*x3 + x11;
  double x15 = p[0]*std::pow(p[2], 2);
  double x16 = -x14*x15;
  double x17 = 10*x9;
  double x18 = std::pow(x8, -2);
  double x19 = 6*f[0];
  double x20 = 12*f[0]*x6 - 12*x0*x6 + 12*x0 - 6*x12 - x19*x3 - x19;
  double x21 = x18*x20;
  double x22 = 2*f[0];
  double x23 = x21*x22;
  double x24 = 2*x21;
  double x25 = 4*x21;
  double x26 = x10 + 12*x7*x9;
  double x27 = -4*f[0]*x18*x20*x6 - 4*x0*x18*x20 + x0*x25*x6 + x1*x17 + x10*x3 + x12*x24 + x23*x3 + x23 + x26 - 12*x5*x9 - 4*x6*x9;
  double x28 = -x27;
  double x29 = f[0]*x15;
  double x30 = x5 + x6 - 1;
  double x31 = std::pow(x30, 2);
  double x32 = std::pow(x30, -3);
  double x33 = x31*x32;
  double x34 = x15*x33;
  double x35 = 2*x34*x5;
  double x36 = p[1]*x15;
  double x37 = 2*x33*x36*x6;
  double x38 = x15 + x35 - x37;
  double x39 = f[1]*x36;
  double x40 = 4*x33*(f[1]*x22*x34 + x22*x33*x39);
  double x41 = f[1]*x40;
  double x42 = std::pow(f[1], 3);
  double x43 = x25*x42;
  double x44 = f[1]*x24;
  double x45 = std::pow(f[1], 5);
  double x46 = f[1]*x5;
  double x47 = 8*f[0]*x9;
  double x48 = 8*f[1]*x0*x9 - f[1]*x47 + x42*x47;
  double x49 = x1*x44 + x24*x45 - x25*x46 + x43*x5 - x43 + x44 + x48;
  double x50 = f[1]*x10;
  double x51 = -4*f[1]*x5*x9 + x1*x50 + x10*x45 + x13*x42*x5 - 4*x42*x9 + x50;
  double x52 = -x15*x51;
  double x53 = -x49;
  double x54 = f[0]*x40;
  double x55 = -x54;
  double x56 = -x35 + x36 + x37;
  double x57 = 2*f[1]*x31*x32*x38 - 2*f[1]*x33*x56 - x28*x39 - x29*x53 - x52 - x55;
  double x58 = -x41;
  double x59 = 6*f[1];
  double x60 = 12*f[1]*x5 - x1*x59 - 12*x42*x5 + 12*x42 - 6*x45 - x59;
  double x61 = x18*x60;
  double x62 = x22*x61;
  double x63 = 2*x61;
  double x64 = 4*x61;
  double x65 = -4*f[0]*x18*x6*x60 - 4*x0*x18*x60 + x0*x6*x64 + x12*x63 + x3*x62 + x48 + x62;
  double x66 = -x65;
  double x67 = x42*x64;
  double x68 = f[1]*x63;
  double x69 = x1*x10 + x1*x68 - x13*x5 + x17*x3 + x26 + x45*x63 - x46*x64 + x5*x67 - 12*x6*x9 - x67 + x68;
  double x70 = -x69;
  double x71 = 2*f[0]*x31*x32*x56 - p[1]*x16 - x22*x33*x38 - x29*x70 - x39*x66 - x58;

 sum[0]=4*f[0]*x31*x32*x38 - x16 - x28*x29 - x39*x49 - x41;

 sum[4]=4*f[1]*x31*x32*x38 - x29*x66 - x36*x51 - x39*x69 - x55;

 sum[2]=x57;

 sum[6]=x71;

 sum[1]=x57;

 sum[5]=x71;

 sum[3]=4*f[0]*x31*x32*x56 - x14*x15 - x27*x29 - x39*x53 - x58;

 sum[7]=4*f[1]*x31*x32*x56 - p[1]*x52 - x29*x65 - x39*x70 - x54;
       
        return sum;
	}
    
    int getnF()
    {
        return nF;
    }
    
    int getnP()
    {
        return nP;
    }

};
#endif