#!/usr/bin/env python

import sys;

""" this program read in the data from PoleData13.dat and then convert them into a .cpp file 
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
                + tab3 + "zweight.push_back(wvec1);" + newline\
                + tab3 + "fweight.push_back(wvec2);" + newline\
                + tab3 + "eweight.push_back(wvec3);" + newline

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
"std::vector< std::complex<double> > &out_zshift, std::vector < std::complex<double> > &out_zweight, " + \
"std::vector< std::complex<double> > &out_fweight, std::vector < std::complex<double> > &out_eweight  ){"+ newline + \
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
"                   if( inputMethod != 2) {  "+ newline + \
"                     out_fweight = fweight[i];"+ newline + \
"                     out_eweight = eweight[i];"+ newline + \
"                   }                        "+ newline + \
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
    """ read the PoleData13.dat and split them into list data """
    try:
        with open("PoleData13.dat", "r") as infile:
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

            wx = [];wy = [];zx = [];zy = []; zx2 =[]; zy2 =[]; zx3 = []; zy3= []
            for j in range(numPole):
                wx.append(data[i+j+1][0])
                wy.append(data[i+j+1][1])
                zx.append(data[i+j+1][2])
                zy.append(data[i+j+1][3])
                zx2.append(data[i+j+1][4])
                zy2.append(data[i+j+1][5])
                zx3.append(data[i+j+1][6])
                zy3.append(data[i+j+1][7])

            #zshiftStr = tab3+ complex_vector + 'zvec'+' = {\n' + tab3
            zshiftStr = tab3+ type_convert + 'vec_temp[] = {\n' + tab3
            for j in range(numPole):
                zshiftStr +=  type_convert + '(' + str(wx[j]) + ' , ' + str(wy[j]) + '), \n' + tab3
            zshiftStr += '}'+ semicolon + newline
            zshiftStr += tab3+'std::vector< std::complex<double> >zvec (vec_temp, vec_temp + sizeof(vec_temp)/ sizeof(std::complex<double>));' + newline

            #zweightStr = tab3+complex_vector + 'wvec'+' = {\n'+ tab3
            zweightStr = tab3+type_convert + 'vec_temp1[] = {\n'+ tab3

            for j in range(numPole):
                zweightStr +=  type_convert + '(' + str(zx[j]) + ' , ' + str(zy[j]) + '), \n' + tab3
            zweightStr += '}' + semicolon + newline
            zweightStr += tab3+'std::vector< std::complex<double> >wvec1 (vec_temp1, vec_temp1 + sizeof(vec_temp1)/ sizeof(std::complex<double>));' + newline

            fweight = tab3+ type_convert + 'vec_temp2[] = {\n' + tab3
            for j in range(numPole):
                fweight  +=  type_convert + '(' + str(zx2[j]) + ' , ' + str(zy2[j]) + '), \n' + tab3
            fweight  += '}'+ semicolon + newline
            fweight  += tab3+'std::vector< std::complex<double> >wvec2 (vec_temp2, vec_temp2 + sizeof(vec_temp2)/ sizeof(std::complex<double>));' + newline

            eweight = tab3+ type_convert + 'vec_temp3[] = {\n' + tab3
            for j in range(numPole):
                eweight  +=  type_convert + '(' + str(zx3[j]) + ' , ' + str(zy3[j]) + '), \n' + tab3
            eweight  += '}'+ semicolon + newline
            eweight  += tab3+'std::vector< std::complex<double> >wvec3 (vec_temp3, vec_temp3 + sizeof(vec_temp3)/ sizeof(std::complex<double>));' + newline




            print tab2+"{", newline,numPoleStr, methodStr, betaStr, zshiftStr, zweightStr, fweight, eweight, push_operations+ tab2+"}",newline
    
    print CPP_CONTENT_REST 

readnumbers()
print include_header
constructPoleCPP()       
#print CPP_TEST_CODE
