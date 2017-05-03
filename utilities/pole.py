#!/usr/bin/env python

import sys;

""" this program read in the data from PoleData.dat and then convert them into a .cpp file 
    note that this the last line is for test purpose, if you do not want to test this program
    just comment it out. """

data = []

beta      = "double b"
npHead    = "int np"
equal     = " = "
methodHead= "int m"
betaEStr  = " beta ["
shiftStr  = "    w ["
semicolon = ";"
tab       = "    "
tab2      = "        "
tab3      = "            "
newline   = "\n"
complex_vector = "std::vector< std::complex<double> >" 
type_convert   = " std::complex<double> "

push_operations = tab3 + "method.push_back(m);"    + newline\
                + tab3 + "beta.push_back(b);"       + newline\
                + tab3 + "numPole.push_back(np);"   + newline\
                + tab3 + "zshift.push_back(zvec);"  + newline\
                + tab3 + "zweight.push_back(wvec);" + newline

include_header =  \
"#include \"pexsi/environment.hpp\" "+ newline + \
"#include \"pexsi/getPole.hpp\" "+ newline + \
""+ newline +\
"   poleClass::poleClass(){ "+ newline 
#"class poleClass { "+ newline + \
#"    private:"+ newline + \
#"        std::vector<int> method;"+ newline + \
#"        std::vector<int> numPole;"+ newline + \
#"        std::vector<double> beta;"+ newline + \
#"        std::vector< std::vector< std::complex<double> > >  zshift;"+ newline + \
#"        std::vector< std::vector< std::complex<double> > >  zweight;"+ newline + \
#"    "+ newline + \
#"    public:"+ newline + \
#"        poleClass::poleClass(){"+ newline + \

CPP_CONTENT_REST =  \
"   };"+ newline + \
"   poleClass::~poleClass(){};"+ newline + \
""+ newline + \
"   // give me a set of parameter, return a set of pole and zshift."+ newline + \
"   bool poleClass::getPole ( int inputMethod, int inputPole, double inputBeta,  " + \
"std::vector< std::complex<double> > &out_zshift, std::vector < std::complex<double> > &out_zweight ){"+ newline + \
"        for ( int i = 0; i < method.size(); i++ ){"+ newline + \
"            if(inputMethod == method[i] ){ "+ newline + \
"               if(numPole[i] >= inputPole && beta[i] >= inputBeta){"+ newline + \
"                 if(i < method.size() - 1)" + newline + \
"                   if((numPole[i+1] >= inputPole) && (beta[i+1] >= inputBeta)){" + newline + \
"                      // assume the numPole and beta are in assending or desending order " + newline + \
"                      if( (numPole[i+1] < numPole[i]) || (beta[i+1] < beta[i]))" + newline + \
"                          continue;" +newline +\
"                   }" + newline+ \
"                   out_zshift = zshift[i];"+ newline + \
"                   out_zweight = zweight[i];"+ newline + \
"                   return true;"+ newline + \
"               }"+ newline + \
"            }"+ newline + \
"        }"+ newline + \
"        return false;"+ newline + \
"   }"+ newline + \
"   "+ newline 
#"};"+ newline 

CPP_TEST_CODE = \
"int main()" + newline + \
"{" + newline + \
"    poleClass pole1;" + newline + \
"    std::vector< std::complex<double> > zshift;" + newline + \
"    std::vector< std::complex<double> > zweight;" + newline + \
"    pole1.getPole(2, 10,10, zshift, zweight);" + newline + \
"" + newline + \
"    std::cout << \" method is \" << 2 << \" beta is :\" << 10 << \" numPole \" << 10 << std::endl; " + newline  + \
"    std::cout << std::endl<< \" zshift returns:\" << std::endl;" + newline + \
"    for(int i = 0; i < zshift.size(); i++)" + newline + \
"       std::cout<< std::setw(10) << zshift[i].real() << std::setw(10) << zshift[i].imag() << std::endl;" + newline + \
"    pole1.getPole(1, 20,20, zshift, zweight);" + newline + \
"" + newline + \
"    std::cout << \" method is \" << 1 << \" beta is :\" << 20 << \" numPole \" << 20 << std::endl;" + newline + \
"    std::cout << std::endl<< \" zshift returns:\" << std::endl;" + newline + \
"    for(int i = 0; i < zshift.size(); i++)" + newline + \
"       std::cout<< std::setw(10) << zshift[i].real() << std::setw(10) << zshift[i].imag() << std::endl;" + newline + \
"    return 0;" + newline + \
"}  " 

def readnumbers():
    """ read the PoleData.dat and split them into list data """
    try:
        with open("PoleData.dat", "r") as infile:
            for line in infile:
                data.append(filter(None, line.strip('\n').split(' ')))
    except:
        pass
    return data


# this is the main function of the pole.py
def constructPoleCPP():
    """ construct the wx, wy, zx, zy  then write them to file"""
    for i in range(len(data)):
        if len(data[i]) == 3:  # control line
            method  = data[i][0]
            numPole = int(data[i][1])
            betaEnergy = data[i][2]
            methodStr  = tab3 + methodHead + equal + method + semicolon + newline
            betaStr    = tab3 + beta +  equal + betaEnergy + semicolon + newline
            numPoleStr = tab3 + npHead +  equal + str(numPole) + semicolon + newline

            wx = [];wy = [];zx = [];zy = []
            for j in range(numPole):
                wx.append(data[i+j+1][0])
                wy.append(data[i+j+1][1])
                zx.append(data[i+j+1][2])
                zy.append(data[i+j+1][3])

            #zshiftStr = tab3+ complex_vector + 'zvec'+' = {\n' + tab3
            zshiftStr = tab3+ type_convert + 'zvec1[] = {\n' + tab3
            for j in range(numPole):
                zshiftStr +=  type_convert + '(' + str(wx[j]) + ' , ' + str(wy[j]) + '), \n' + tab3
            zshiftStr += '}'+ semicolon + newline
            zshiftStr += tab3+'std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));' + newline

            #zweightStr = tab3+complex_vector + 'wvec'+' = {\n'+ tab3
            zweightStr = tab3+type_convert + 'wvec1[] = {\n'+ tab3

            for j in range(numPole):
                zweightStr +=  type_convert + '(' + str(zx[j]) + ' , ' + str(zy[j]) + '), \n' + tab3
            zweightStr += '}' + semicolon + newline
            zweightStr += tab3+'std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));' + newline
            print tab2+"{", newline,numPoleStr, methodStr, betaStr, zshiftStr, zweightStr, push_operations+ tab2+"}",newline
    
    print CPP_CONTENT_REST 

readnumbers()
print include_header
constructPoleCPP()       
#print CPP_TEST_CODE
