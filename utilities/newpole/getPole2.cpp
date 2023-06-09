#include "pexsi/environment.hpp" 
#include "pexsi/getPole.hpp" 

   poleClass::poleClass(){ 

        { 
            int np = 5;
            int m = 2;
            double b = 4.754276176478917e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.946861030831248e+06 , 3.466755326606405e+06), 
             std::complex<double> (-3.194821559260857e+03 , 1.610341217030958e+05), 
             std::complex<double> (-4.297502137712511e+00 , 5.907285762749364e+03), 
             std::complex<double> (-5.778795889148681e-03 , 2.166148672731455e+02), 
             std::complex<double> (-7.536111765206299e-06 , 7.656827074796709e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (3.318064098884963e+06 , -2.022541739446533e+06), 
             std::complex<double> (6.653510807494008e+03 , -1.676182108199611e+05), 
             std::complex<double> (8.952581631955500e+00 , -6.153043298116445e+03), 
             std::complex<double> (1.203899544460975e-02 , -2.256483267707010e+02), 
             std::complex<double> (1.613855688985912e-05 , -8.578097024537158e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 5;
            int m = 2;
            double b = 1.929359655736823e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.166555857792123e+05 , 1.643367718523635e+05), 
             std::complex<double> (-7.514143203926686e+02 , 1.615703144671552e+04), 
             std::complex<double> (-3.735604646002479e+00 , 1.140431725572104e+03), 
             std::complex<double> (-1.853137934183301e-02 , 8.033020123555853e+01), 
             std::complex<double> (-8.235334290060451e-05 , 5.300406238408231e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.505030088371715e+05 , -5.259160857469719e+04), 
             std::complex<double> (1.264161757745586e+03 , -1.356170840398750e+04), 
             std::complex<double> (6.293858648115211e+00 , -9.607057064381456e+02), 
             std::complex<double> (3.124454004433670e-02 , -6.770929837956629e+01), 
             std::complex<double> (1.517929579459690e-04 , -5.110489998965694e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 5;
            int m = 2;
            double b = 2.352699225192334e+04;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.862211631189100e+04 , 2.163767382592681e+04), 
             std::complex<double> (-3.034717531824400e+02 , 3.631659031140577e+03), 
             std::complex<double> (-3.608328744822498e+00 , 3.973734595337082e+02), 
             std::complex<double> (-4.255372563255384e-02 , 4.321103903220411e+01), 
             std::complex<double> (-3.938722873497464e-04 , 4.296822012346186e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.878579439536142e+04 , -2.830082326877334e+03), 
             std::complex<double> (4.264339341259439e+02 , -2.533758290280452e+03), 
             std::complex<double> (5.092612048899103e+00 , -2.803844017083945e+02), 
             std::complex<double> (6.029749726793857e-02 , -3.053576668190605e+01), 
             std::complex<double> (6.694655782523798e-04 , -3.701416529003495e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 5;
            int m = 2;
            double b = 1.765915008963304e+03;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.947377806674931e+03 , 1.674466524728663e+03), 
             std::complex<double> (-1.047929084361665e+02 , 5.865023929298897e+02), 
             std::complex<double> (-3.762497415441033e+00 , 1.128940039023474e+02), 
             std::complex<double> (-1.288451256862734e-01 , 2.122897595164211e+01), 
             std::complex<double> (-2.382903724382977e-03 , 3.517696741585170e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.288627566572213e+03 , 1.953118537029985e+02), 
             std::complex<double> (1.095143820548338e+02 , -2.966666335909262e+02), 
             std::complex<double> (4.006000061451483e+00 , -5.996285636016435e+01), 
             std::complex<double> (1.415274613924904e-01 , -1.130840867511528e+01), 
             std::complex<double> (3.931476553003528e-03 , -2.604882313795629e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 5;
            int m = 2;
            double b = 1.000579208201004e+02;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.547238643537001e+02 , 8.105637703696272e+01), 
             std::complex<double> (-3.406686599486403e+01 , 7.481751505704190e+01), 
             std::complex<double> (-4.504086450290978e+00 , 3.029656991427235e+01), 
             std::complex<double> (-4.145553155245011e-01 , 1.117658457886601e+01), 
             std::complex<double> (-6.486243326569170e-03 , 3.165000076234382e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (4.609342039725342e+01 , 3.201479941137274e+01), 
             std::complex<double> (2.041059397555925e+01 , -1.755755921707144e+01), 
             std::complex<double> (3.067084432942332e+00 , -9.516725097319602e+00), 
             std::complex<double> (3.731958350013732e-01 , -3.572807519232854e+00), 
             std::complex<double> (1.315295794894268e-02 , -2.048979911103968e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 5;
            int m = 2;
            double b = 1.209274548499165e+01;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.554090010277767e+01 , 8.638640679328239e+00), 
             std::complex<double> (-1.246440939244841e+01 , 1.609925983411248e+01), 
             std::complex<double> (-2.525722116114511e+00 , 1.336080612620282e+01), 
             std::complex<double> (1.081718564305512e-01 , 9.270389934959379e+00), 
             std::complex<double> (9.225429531107947e-05 , 3.141649048770424e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (3.033611448442282e+00 , 4.764487587422160e+00), 
             std::complex<double> (4.183144608690097e+00 , 1.482503127170896e-01), 
             std::complex<double> (2.076676551343867e+00 , -1.223030987690861e+00), 
             std::complex<double> (-5.341903753127276e-02 , -1.689775951588069e+00), 
             std::complex<double> (-3.135698211663160e-04 , -2.000162109432435e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 10;
            int m = 2;
            double b = 1.290985308157311e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.580972615908010e+06 , 1.203561258897960e+06), 
             std::complex<double> (-1.278275496987531e+05 , 5.503393895202582e+05), 
             std::complex<double> (-6.698848876743868e+03 , 1.291650688368753e+05), 
             std::complex<double> (-3.418913221857892e+02 , 2.921746870195993e+04), 
             std::complex<double> (-1.742544928216568e+01 , 6.596580351519240e+03), 
             std::complex<double> (-8.880680791697854e-01 , 1.489201251999560e+03), 
             std::complex<double> (-4.525215306773906e-02 , 3.361928044929759e+02), 
             std::complex<double> (-2.298812720870707e-03 , 7.590558678879101e+01), 
             std::complex<double> (-1.100909893889450e-04 , 1.716220242977646e+01), 
             std::complex<double> (-2.289831112613212e-06 , 3.360647540251634e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (8.701120180884706e+05 , 2.402813672182368e+05), 
             std::complex<double> (1.175901570695602e+05 , -2.394756773944366e+05), 
             std::complex<double> (6.337463238907602e+03 , -6.093413403216130e+04), 
             std::complex<double> (3.239120631959876e+02 , -1.383859994824320e+04), 
             std::complex<double> (1.651030185618987e+01 , -3.125047689860075e+03), 
             std::complex<double> (8.414388270961817e-01 , -7.054973718325768e+02), 
             std::complex<double> (4.288318148496536e-02 , -1.592667565848489e+02), 
             std::complex<double> (2.185485516150124e-03 , -3.595060328223087e+01), 
             std::complex<double> (1.110150185677044e-04 , -8.117993876151560e+00), 
             std::complex<double> (3.870864037752220e-06 , -2.371501741456056e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 10;
            int m = 2;
            double b = 1.011393241376369e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.428603078252365e+05 , 8.854625242936362e+04), 
             std::complex<double> (-2.096982075063316e+04 , 6.088418821626271e+04), 
             std::complex<double> (-1.928152539915390e+03 , 1.943088566651270e+04), 
             std::complex<double> (-1.686646753695693e+02 , 5.772668827112795e+03), 
             std::complex<double> (-1.468849355918116e+01 , 1.704216565544883e+03), 
             std::complex<double> (-1.278554065937448e+00 , 5.028496502931816e+02), 
             std::complex<double> (-1.111599977344547e-01 , 1.483790496092092e+02), 
             std::complex<double> (-9.537726369604812e-03 , 4.383014383939805e+01), 
             std::complex<double> (-7.037281719416144e-04 , 1.305469865817097e+01), 
             std::complex<double> (-1.421408210600723e-05 , 3.215478433048512e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (5.660149505805381e+04 , 2.811926983163520e+04), 
             std::complex<double> (1.533514053977184e+04 , -1.962129774604957e+04), 
             std::complex<double> (1.490363817974273e+03 , -7.435596941967240e+03), 
             std::complex<double> (1.310048542130718e+02 , -2.239951076563770e+03), 
             std::complex<double> (1.141375733465775e+01 , -6.620771714240800e+02), 
             std::complex<double> (9.936511611577747e-01 , -1.953707146259180e+02), 
             std::complex<double> (8.649872984867023e-02 , -5.763767097170560e+01), 
             std::complex<double> (7.529090987797448e-03 , -1.698662423589999e+01), 
             std::complex<double> (6.439158240058465e-04 , -5.009045346064647e+00), 
             std::complex<double> (2.611935108885642e-05 , -2.139547045657126e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 10;
            int m = 2;
            double b = 1.770766131099239e+04;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.735487926298620e+04 , 1.430831195405075e+04), 
             std::complex<double> (-6.052868944139718e+03 , 1.319995629014940e+04), 
             std::complex<double> (-8.332473897685190e+02 , 5.323088658739502e+03), 
             std::complex<double> (-1.064222880882963e+02 , 1.922595865092381e+03), 
             std::complex<double> (-1.345888836775062e+01 , 6.846647159118429e+02), 
             std::complex<double> (-1.699046933802864e+00 , 2.433981937143934e+02), 
             std::complex<double> (-2.135110433905035e-01 , 8.654964106249528e+01), 
             std::complex<double> (-2.589746005590427e-02 , 3.088642287785381e+01), 
             std::complex<double> (-2.357685107190892e-03 , 1.120947014647582e+01), 
             std::complex<double> (-3.704363969261008e-05 , 3.165754611818758e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (8.137790793302991e+03 , 5.650698761644630e+03), 
             std::complex<double> (3.603512235862979e+03 , -3.103027570633254e+03), 
             std::complex<double> (5.417765597007763e+02 , -1.688126686747975e+03), 
             std::complex<double> (6.997870141348011e+01 , -6.301622624525596e+02), 
             std::complex<double> (8.863287854229052e+00 , -2.253290314681059e+02), 
             std::complex<double> (1.119811581063546e+00 , -8.013481424618158e+01), 
             std::complex<double> (1.414348171048943e-01 , -2.846599076368053e+01), 
             std::complex<double> (1.785365037973103e-02 , -1.007802731807097e+01), 
             std::complex<double> (2.124014496053039e-03 , -3.602552308819546e+00), 
             std::complex<double> (7.498359129652677e-05 , -2.050381159183316e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 10;
            int m = 2;
            double b = 1.309935703693513e+03;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.274457093148581e+03 , 8.568804108514295e+02), 
             std::complex<double> (-9.187129229487347e+02 , 1.241847095125482e+03), 
             std::complex<double> (-2.417599090810289e+02 , 7.546960297935219e+02), 
             std::complex<double> (-5.614604098767667e+01 , 3.778615998140903e+02), 
             std::complex<double> (-1.263601231666925e+01 , 1.810151451948116e+02), 
             std::complex<double> (-2.804984046956645e+00 , 8.593503637047091e+01), 
             std::complex<double> (-6.020694504085211e-01 , 4.087468806042249e+01), 
             std::complex<double> (-1.101531568024575e-01 , 1.975481603753779e+01), 
             std::complex<double> (-9.249622639135294e-03 , 9.713004727279520e+00), 
             std::complex<double> (-5.098768728747187e-05 , 3.142848009151120e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (3.766119122544866e+02 , 4.288970428869044e+02), 
             std::complex<double> (3.500773456138319e+02 , -1.070948056163151e+02), 
             std::complex<double> (1.096326999220600e+02 , -1.535266579530764e+02), 
             std::complex<double> (2.648933059164146e+01 , -8.710117864869949e+01), 
             std::complex<double> (6.026108975876521e+00 , -4.281107786365286e+01), 
             std::complex<double> (1.352377667258937e+00 , -2.039245273804895e+01), 
             std::complex<double> (3.025591460325434e-01 , -9.609751665571704e+00), 
             std::complex<double> (6.679308594421993e-02 , -4.471462469775012e+00), 
             std::complex<double> (9.879394028859808e-03 , -2.332857725660222e+00), 
             std::complex<double> (1.261386207797346e-04 , -2.003158096401944e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 10;
            int m = 2;
            double b = 3.172305323816737e+02;
             std::complex<double> zvec1[] = {
             std::complex<double> (-5.809205613865953e+02 , 1.712116651384128e+02), 
             std::complex<double> (-3.192914878221479e+02 , 3.157981973922744e+02), 
             std::complex<double> (-1.228707736463758e+02 , 2.502456547118489e+02), 
             std::complex<double> (-4.073432649941893e+01 , 1.556500960342447e+02), 
             std::complex<double> (-1.273734591445475e+01 , 8.975596322505282e+01), 
             std::complex<double> (-3.813186932424677e+00 , 5.072137998496303e+01), 
             std::complex<double> (-1.027939809349396e+00 , 2.878927093668447e+01), 
             std::complex<double> (-1.846368604191310e-01 , 1.676449656605937e+01), 
             std::complex<double> (-8.323591044985035e-03 , 9.459938671536779e+00), 
             std::complex<double> (-1.505288530203722e-05 , 3.141645614487579e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (6.077593094379878e+01 , 9.420513079285334e+01), 
             std::complex<double> (8.348916402748404e+01 , 9.929637216513899e-01), 
             std::complex<double> (4.119366193170025e+01 , -3.171509381847607e+01), 
             std::complex<double> (1.479181796837569e+01 , -2.611492948096926e+01), 
             std::complex<double> (4.792005490178411e+00 , -1.616552764687287e+01), 
             std::complex<double> (1.503225379305250e+00 , -9.235793010792859e+00), 
             std::complex<double> (4.660631616769804e-01 , -5.101101538909798e+00), 
             std::complex<double> (1.274852078225833e-01 , -2.818387233161705e+00), 
             std::complex<double> (1.149867943559144e-02 , -2.051749410099756e+00), 
             std::complex<double> (4.246018475094892e-05 , -2.000150776505864e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 15;
            int m = 2;
            double b = 2.166546948895497e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-3.385879641866369e+06 , 1.726210205061707e+06), 
             std::complex<double> (-7.922212469180174e+05 , 1.658903021713648e+06), 
             std::complex<double> (-1.154176033095820e+05 , 6.921301461044513e+05), 
             std::complex<double> (-1.553772831992236e+04 , 2.569858018791626e+05), 
             std::complex<double> (-2.069069447518909e+03 , 9.392688282197727e+04), 
             std::complex<double> (-2.751255803032834e+02 , 3.425777964848087e+04), 
             std::complex<double> (-3.657655141499721e+01 , 1.249129065893094e+04), 
             std::complex<double> (-4.862534038900822e+00 , 4.554486110796908e+03), 
             std::complex<double> (-6.464213496231945e-01 , 1.660618643218002e+03), 
             std::complex<double> (-8.592669116100400e-02 , 6.054867959546815e+02), 
             std::complex<double> (-1.141395426989474e-02 , 2.207869595142982e+02), 
             std::complex<double> (-1.508164195324225e-03 , 8.055575273372025e+01), 
             std::complex<double> (-1.913738517796651e-04 , 2.951522815522122e+01), 
             std::complex<double> (-1.774919765884570e-05 , 1.101184090070446e+01), 
             std::complex<double> (-2.636634171792006e-07 , 3.161572777755021e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (9.635623239800008e+05 , 6.993662871740375e+05), 
             std::complex<double> (4.565627974088923e+05 , -3.690014177845073e+05), 
             std::complex<double> (7.306277918570081e+04 , -2.129776265811492e+05), 
             std::complex<double> (9.960607732572666e+03 , -8.207047339297933e+04), 
             std::complex<double> (1.328624427353701e+03 , -3.014229259682151e+04), 
             std::complex<double> (1.767075335463783e+02 , -1.100082553212721e+04), 
             std::complex<double> (2.349307381994124e+01 , -4.011534460335000e+03), 
             std::complex<double> (3.123217856103624e+00 , -1.462673362452836e+03), 
             std::complex<double> (4.152042060442371e-01 , -5.333068338400971e+02), 
             std::complex<double> (5.519767866959679e-02 , -1.944471431488452e+02), 
             std::complex<double> (7.338035018694221e-03 , -7.089113145949230e+01), 
             std::complex<double> (9.755186371577616e-04 , -2.583037302487413e+01), 
             std::complex<double> (1.296002154544079e-04 , -9.376270006378212e+00), 
             std::complex<double> (1.604941909550027e-05 , -3.446970226546104e+00), 
             std::complex<double> (5.421813582989691e-07 , -2.042299222653249e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 15;
            int m = 2;
            double b = 1.601767940704642e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.703524319835728e+05 , 1.123075595283472e+05), 
             std::complex<double> (-9.356434578914900e+04 , 1.445897671219019e+05), 
             std::complex<double> (-2.065312922187001e+04 , 7.823439465627717e+04), 
             std::complex<double> (-4.079526205677131e+03 , 3.572944274484625e+04), 
             std::complex<double> (-7.878027820667667e+02 , 1.578346927723723e+04), 
             std::complex<double> (-1.514669122122054e+02 , 6.927706044029221e+03), 
             std::complex<double> (-2.909703958808262e+01 , 3.036969796221203e+03), 
             std::complex<double> (-5.588538408212697e+00 , 1.331036389659291e+03), 
             std::complex<double> (-1.073189797099922e+00 , 5.833472157869685e+02), 
             std::complex<double> (-2.059454462438908e-01 , 2.556817422155146e+02), 
             std::complex<double> (-3.937860304792234e-02 , 1.121179493277981e+02), 
             std::complex<double> (-7.387205093537825e-03 , 4.928300299709128e+01), 
             std::complex<double> (-1.246295164609699e-03 , 2.191396271182843e+01), 
             std::complex<double> (-1.159700831130521e-04 , 9.965604514487502e+00), 
             std::complex<double> (-9.355796827656802e-07 , 3.145038541211152e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (5.350042922869226e+04 , 5.328211658562325e+04), 
             std::complex<double> (4.100212181764019e+04 , -1.841505731653905e+04), 
             std::complex<double> (1.047334334363492e+04 , -1.845416745864924e+04), 
             std::complex<double> (2.128119871225179e+03 , -9.197792120560949e+03), 
             std::complex<double> (4.132087648398298e+02 , -4.128963256574471e+03), 
             std::complex<double> (7.952887758707516e+01 , -1.817850196095971e+03), 
             std::complex<double> (1.528079074130672e+01 , -7.973766736195694e+02), 
             std::complex<double> (2.935119059861098e+00 , -3.495090168341490e+02), 
             std::complex<double> (5.637395684940451e-01 , -1.531744583664404e+02), 
             std::complex<double> (1.082744804936059e-01 , -6.712182572379126e+01), 
             std::complex<double> (2.079565354768343e-02 , -2.939924269841751e+01), 
             std::complex<double> (3.993955855684405e-03 , -1.284649312230424e+01), 
             std::complex<double> (7.634114708543606e-04 , -5.562647550836926e+00), 
             std::complex<double> (1.146407092880542e-04 , -2.572673341944408e+00), 
             std::complex<double> (2.184467006656409e-06 , -2.008216560524396e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 15;
            int m = 2;
            double b = 1.110657715232228e+03;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.102013583493623e+03 , 4.854486997637087e+02), 
             std::complex<double> (-1.426932638084718e+03 , 1.059889853957434e+03), 
             std::complex<double> (-7.373954335086265e+02 , 1.043612485056151e+03), 
             std::complex<double> (-3.290209767739146e+02 , 7.877117552200033e+02), 
             std::complex<double> (-1.372338807633277e+02 , 5.341329130011428e+02), 
             std::complex<double> (-5.560171963187359e+01 , 3.468534670173635e+02), 
             std::complex<double> (-2.223138070151641e+01 , 2.214122743727300e+02), 
             std::complex<double> (-8.810308831015886e+00 , 1.404504350390194e+02), 
             std::complex<double> (-3.447398125067049e+00 , 8.899174697876053e+01), 
             std::complex<double> (-1.309501923286421e+00 , 5.655490167853964e+01), 
             std::complex<double> (-4.573557444135257e-01 , 3.628965934123512e+01), 
             std::complex<double> (-1.210035395684593e-01 , 2.380282928497694e+01), 
             std::complex<double> (-1.317166632301192e-02 , 1.585258274196877e+01), 
             std::complex<double> (-1.988156522421698e-04 , 9.426469287518469e+00), 
             std::complex<double> (-1.166695585979868e-07 , 3.141593548270519e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.378991690184070e+02 , 2.826512095221750e+02), 
             std::complex<double> (2.485688765195191e+02 , 7.503355866571955e+01), 
             std::complex<double> (1.763109435064582e+02 , -6.244070835866666e+01), 
             std::complex<double> (8.900682012933653e+01 , -8.790616164206112e+01), 
             std::complex<double> (3.900650802628571e+01 , -7.081920632000730e+01), 
             std::complex<double> (1.613305718034594e+01 , -4.890285378312584e+01), 
             std::complex<double> (6.518824995941943e+00 , -3.193738790389188e+01), 
             std::complex<double> (2.609625251922411e+00 , -2.040404511970779e+01), 
             std::complex<double> (1.040914809451731e+00 , -1.290556079321114e+01), 
             std::complex<double> (4.148036489526483e-01 , -8.101505849609401e+00), 
             std::complex<double> (1.653484677659098e-01 , -5.026674034803303e+00), 
             std::complex<double> (6.246784231319495e-02 , -3.084578618700840e+00), 
             std::complex<double> (1.261172748280164e-02 , -2.153263048067060e+00), 
             std::complex<double> (3.515563501331491e-04 , -2.003083586937637e+00), 
             std::complex<double> (3.673502695048175e-07 , -2.000002830271054e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 15;
            int m = 2;
            double b = 2.635002559079762e+02;
             std::complex<double> zvec1[] = {
             std::complex<double> (-5.114065966622117e+02 , 8.936507958623726e+01), 
             std::complex<double> (-4.062268650550005e+02 , 2.217612088065892e+02), 
             std::complex<double> (-2.679562277212592e+02 , 2.639459577395004e+02), 
             std::complex<double> (-1.558649586430436e+02 , 2.412763662304596e+02), 
             std::complex<double> (-8.403314427814890e+01 , 1.940707686957161e+02), 
             std::complex<double> (-4.332772177203707e+01 , 1.464156872544164e+02), 
             std::complex<double> (-2.168279729232532e+01 , 1.070573072146486e+02), 
             std::complex<double> (-1.053793311572390e+01 , 7.721260254532646e+01), 
             std::complex<double> (-4.887168088158905e+00 , 5.553087645137526e+01), 
             std::complex<double> (-2.038222987249636e+00 , 4.020387551153076e+01), 
             std::complex<double> (-6.288661726578217e-01 , 2.961642336847731e+01), 
             std::complex<double> (-8.439255346704523e-02 , 2.211157840630456e+01), 
             std::complex<double> (-2.315587080161116e-03 , 1.571023862887294e+01), 
             std::complex<double> (-7.859026601898820e-06 , 9.424784108720845e+00), 
             std::complex<double> (-1.309314700314198e-09 , 3.141592654510955e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.945281468416921e+01 , 5.407265327360384e+01), 
             std::complex<double> (4.306123267366677e+01 , 2.781158521875630e+01), 
             std::complex<double> (4.168085348003856e+01 , 7.827216221003346e-01), 
             std::complex<double> (2.910231099008997e+01 , -1.292226320855806e+01), 
             std::complex<double> (1.722113834643294e+01 , -1.588039862394899e+01), 
             std::complex<double> (9.359133482356702e+00 , -1.404048337348538e+01), 
             std::complex<double> (4.871289478653694e+00 , -1.097661233897070e+01), 
             std::complex<double> (2.481899557929571e+00 , -8.104760341445957e+00), 
             std::complex<double> (1.253008775579792e+00 , -5.797195381904102e+00), 
             std::complex<double> (6.305089163750802e-01 , -4.045004019026272e+00), 
             std::complex<double> (2.940170968520652e-01 , -2.777421793167176e+00), 
             std::complex<double> (7.221298359039841e-02 , -2.119243333659646e+00), 
             std::complex<double> (3.370262932733660e-03 , -2.003530196218343e+00), 
             std::complex<double> (1.731008268421628e-05 , -2.000013858491859e+00), 
             std::complex<double> (4.615440311858388e-09 , -2.000000003262454e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 20;
            int m = 2;
            double b = 1.094344975543047e+07;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.843340376824855e+07 , 7.706592728737802e+06), 
             std::complex<double> (-6.310869391131245e+06 , 9.840582271795252e+06), 
             std::complex<double> (-1.376419175767611e+06 , 5.283211086231250e+06), 
             std::complex<double> (-2.688884224216870e+05 , 2.398035423845490e+06), 
             std::complex<double> (-5.137710500981526e+04 , 1.053541628178468e+06), 
             std::complex<double> (-9.775012002737843e+03 , 4.599845946876111e+05), 
             std::complex<double> (-1.858286238681387e+03 , 2.005951131898834e+05), 
             std::complex<double> (-3.532164486406742e+02 , 8.745801132727422e+04), 
             std::complex<double> (-6.713615726438483e+01 , 3.812942507138038e+04), 
             std::complex<double> (-1.276055454167964e+01 , 1.662330568924942e+04), 
             std::complex<double> (-2.425391136915597e+00 , 7.247261002033346e+03), 
             std::complex<double> (-4.609905141385088e-01 , 3.159588607322258e+03), 
             std::complex<double> (-8.761771661289432e-02 , 1.377490100893101e+03), 
             std::complex<double> (-1.665091004308733e-02 , 6.005559175283905e+02), 
             std::complex<double> (-3.162279819539571e-03 , 2.618516648008789e+02), 
             std::complex<double> (-5.985018100927118e-04 , 1.142223461855365e+02), 
             std::complex<double> (-1.112076804732382e-04 , 4.994066120095570e+01), 
             std::complex<double> (-1.863718554027241e-05 , 2.208159867595298e+01), 
             std::complex<double> (-1.741486416532828e-06 , 9.986633955836457e+00), 
             std::complex<double> (-1.437803566715412e-08 , 3.145258673585706e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (3.690251956475755e+06 , 3.641952140790702e+06), 
             std::complex<double> (2.790426848799685e+06 , -1.280800826922408e+06), 
             std::complex<double> (7.030168549697653e+05 , -1.257643244584517e+06), 
             std::complex<double> (1.411886485631274e+05 , -6.216677669302863e+05), 
             std::complex<double> (2.711976968446797e+04 , -2.773984519972565e+05), 
             std::complex<double> (5.164982872229454e+03 , -1.214699096563232e+05), 
             std::complex<double> (9.820801364308887e+02 , -5.300140666708574e+04), 
             std::complex<double> (1.866770648156537e+02 , -2.311067112948131e+04), 
             std::complex<double> (3.548211742885658e+01 , -1.007585568270761e+04), 
             std::complex<double> (6.744088373584011e+00 , -4.392792697917959e+03), 
             std::complex<double> (1.281846600307488e+00 , -1.915126280273689e+03), 
             std::complex<double> (2.436400520924584e-01 , -8.349365756215076e+02), 
             std::complex<double> (4.630856000365761e-02 , -3.640056707263396e+02), 
             std::complex<double> (8.801847957688942e-03 , -1.586922720759284e+02), 
             std::complex<double> (1.672963434280613e-03 , -6.917776119705825e+01), 
             std::complex<double> (3.179793311526134e-04 , -3.014283286278016e+01), 
             std::complex<double> (6.043608081467732e-05 , -1.310436867180438e+01), 
             std::complex<double> (1.143526225643599e-05 , -5.646412859723841e+00), 
             std::complex<double> (1.714068815178379e-06 , -2.591754571863651e+00), 
             std::complex<double> (3.344087840258191e-08 , -2.008710099030617e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 20;
            int m = 2;
            double b = 2.017713616400269e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-3.513344352983811e+06 , 1.309314866250858e+06), 
             std::complex<double> (-1.441997872322787e+06 , 1.921064125313279e+06), 
             std::complex<double> (-3.866462132480044e+05 , 1.182195672718637e+06), 
             std::complex<double> (-9.140064554676184e+04 , 5.978009732383400e+05), 
             std::complex<double> (-2.095264501593017e+04 , 2.887879885537566e+05), 
             std::complex<double> (-4.769195697278365e+03 , 1.380585697064890e+05), 
             std::complex<double> (-1.083797998470479e+03 , 6.584383610131998e+04), 
             std::complex<double> (-2.462020598654694e+02 , 3.138573301310041e+04), 
             std::complex<double> (-5.592404878843604e+01 , 1.495877903690297e+04), 
             std::complex<double> (-1.270272338377913e+01 , 7.129317781697649e+03), 
             std::complex<double> (-2.885301920713309e+00 , 3.397796128974548e+03), 
             std::complex<double> (-6.553551740698245e-01 , 1.619374189762830e+03), 
             std::complex<double> (-1.488416738066650e-01 , 7.717953298672829e+02), 
             std::complex<double> (-3.379145712496281e-02 , 3.678576953499330e+02), 
             std::complex<double> (-7.658782225027643e-03 , 1.753710201919258e+02), 
             std::complex<double> (-1.722960623943232e-03 , 8.369042250538658e+01), 
             std::complex<double> (-3.746827891786966e-04 , 4.011401309500958e+01), 
             std::complex<double> (-6.898404003858887e-05 , 1.955789105661458e+01), 
             std::complex<double> (-5.700188323553038e-06 , 9.692098428416061e+00), 
             std::complex<double> (-3.004041861052978e-08 , 3.142707342514222e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (5.702461395212540e+05 , 6.588271925886495e+05), 
             std::complex<double> (5.407969321635104e+05 , -1.572635932487692e+05), 
             std::complex<double> (1.730983776389920e+05 , -2.363231333492580e+05), 
             std::complex<double> (4.261051037490748e+04 , -1.360884057135358e+05), 
             std::complex<double> (9.858535159363310e+03 , -6.758191296856185e+04), 
             std::complex<double> (2.248688258417920e+03 , -3.250864819868572e+04), 
             std::complex<double> (5.112570644978483e+02 , -1.552595825450053e+04), 
             std::complex<double> (1.161527984543761e+02 , -7.403098068133852e+03), 
             std::complex<double> (2.638441102246710e+01 , -3.528650667788154e+03), 
             std::complex<double> (5.993060094777460e+00 , -1.681773289946056e+03), 
             std::complex<double> (1.361275769354281e+00 , -8.015261931489052e+02), 
             std::complex<double> (3.092023205184201e-01 , -3.820013851445188e+02), 
             std::complex<double> (7.023267030014961e-02 , -1.820566517998944e+02), 
             std::complex<double> (1.595275082353307e-02 , -8.676114496760820e+01), 
             std::complex<double> (3.623532170723744e-03 , -4.133745880859607e+01), 
             std::complex<double> (8.230585301372153e-04 , -1.967543159501197e+01), 
             std::complex<double> (1.869399833383931e-04 , -9.325199391087438e+00), 
             std::complex<double> (4.186515468561719e-05 , -4.370117040865639e+00), 
             std::complex<double> (6.147897251077524e-06 , -2.311820447484295e+00), 
             std::complex<double> (7.477918316363095e-08 , -2.002820348497667e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 20;
            int m = 2;
            double b = 1.069717540448873e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.959655367138746e+05 , 5.745751043909320e+04), 
             std::complex<double> (-1.082740367096186e+05 , 1.063897037408112e+05), 
             std::complex<double> (-4.198433931434520e+04 , 8.468822249692281e+04), 
             std::complex<double> (-1.405697423969717e+04 , 5.285732932751087e+04), 
             std::complex<double> (-4.473794332008163e+03 , 3.052971270797111e+04), 
             std::complex<double> (-1.400791782253232e+03 , 1.720882170413647e+04), 
             std::complex<double> (-4.363590151123119e+02 , 9.626645027053943e+03), 
             std::complex<double> (-1.357121813956949e+02 , 5.372427007557712e+03), 
             std::complex<double> (-4.218658961763070e+01 , 2.996034573285604e+03), 
             std::complex<double> (-1.311150263458092e+01 , 1.670418439345865e+03), 
             std::complex<double> (-4.074518391506804e+00 , 9.312747116316690e+02), 
             std::complex<double> (-1.265861664979148e+00 , 5.192024495356003e+02), 
             std::complex<double> (-3.929599751562658e-01 , 2.894970264151671e+02), 
             std::complex<double> (-1.216724824424809e-01 , 1.614789660905933e+02), 
             std::complex<double> (-3.735891695711960e-02 , 9.018175810892765e+01), 
             std::complex<double> (-1.115298673243625e-02 , 5.056121144817830e+01), 
             std::complex<double> (-3.007072028933906e-03 , 2.869171048694175e+01), 
             std::complex<double> (-5.376108560072794e-04 , 1.674248149639425e+01), 
             std::complex<double> (-2.390396050194380e-05 , 9.459017768247080e+00), 
             std::complex<double> (-4.255114254655321e-08 , 3.141643914629745e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (2.031064741891039e+04 , 3.165835297690931e+04), 
             std::complex<double> (2.807761526179717e+04 , 4.929726536221493e+02), 
             std::complex<double> (1.396364875411875e+04 , -1.062205355878130e+04), 
             std::complex<double> (5.049921187189751e+03 , -8.822916548341604e+03), 
             std::complex<double> (1.646266432090366e+03 , -5.496538220861484e+03), 
             std::complex<double> (5.193299461116089e+02 , -3.168863148169208e+03), 
             std::complex<double> (1.621523465312370e+02 , -1.784966835870925e+03), 
             std::complex<double> (5.046766447557334e+01 , -9.982859239650902e+02), 
             std::complex<double> (1.569173429813937e+01 , -5.570811086742511e+02), 
             std::complex<double> (4.877466861761349e+00 , -3.106576535387721e+02), 
             std::complex<double> (1.515918844363814e+00 , -1.732000574990181e+02), 
             std::complex<double> (4.711341767433641e-01 , -9.655376524215987e+01), 
             std::complex<double> (1.464230204445247e-01 , -5.381830214752993e+01), 
             std::complex<double> (4.550657708728133e-02 , -2.998626717088552e+01), 
             std::complex<double> (1.414331365119132e-02 , -1.668711427648223e+01), 
             std::complex<double> (4.396627657958475e-03 , -9.249929343129404e+00), 
             std::complex<double> (1.364269292114220e-03 , -5.070294118840708e+00), 
             std::complex<double> (3.728103157060478e-04 , -2.803521036589960e+00), 
             std::complex<double> (3.314797149644424e-05 , -2.050438868294707e+00), 
             std::complex<double> (1.202303182951202e-07 , -2.000146100753941e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 20;
            int m = 2;
            double b = 1.222015728799431e+04;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.309805899879820e+04 , 5.394337935974982e+03), 
             std::complex<double> (-1.555512023205595e+04 , 1.170161182960154e+04), 
             std::complex<double> (-7.944568726043814e+03 , 1.141913680938806e+04), 
             std::complex<double> (-3.502200012056795e+03 , 8.546398544466469e+03), 
             std::complex<double> (-1.444426366201753e+03 , 5.753002318023964e+03), 
             std::complex<double> (-5.793929308317200e+02 , 3.711786276856882e+03), 
             std::complex<double> (-2.298135741818555e+02 , 2.354833269772197e+03), 
             std::complex<double> (-9.074583308320911e+01 , 1.484058229105163e+03), 
             std::complex<double> (-3.576614499002655e+01 , 9.328439054036305e+02), 
             std::complex<double> (-1.408356915073149e+01 , 5.857796884490144e+02), 
             std::complex<double> (-5.540794642283571e+00 , 3.677254055179496e+02), 
             std::complex<double> (-2.176285362100652e+00 , 2.308589924572627e+02), 
             std::complex<double> (-8.513896566733504e-01 , 1.450124171586070e+02), 
             std::complex<double> (-3.296773333068714e-01 , 9.122618242677341e+01), 
             std::complex<double> (-1.242020417880811e-01 , 5.761277066152592e+01), 
             std::complex<double> (-4.320931302449264e-02 , 3.674009552907650e+01), 
             std::complex<double> (-1.153408637231899e-02 , 2.394424799537720e+01), 
             std::complex<double> (-1.303181251949702e-03 , 1.587059405402117e+01), 
             std::complex<double> (-2.085902319099095e-05 , 9.426789422361649e+00), 
             std::complex<double> (-1.289228047401315e-08 , 3.141593772410077e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.548492418633355e+03 , 3.134434847243796e+03), 
             std::complex<double> (2.762357706195695e+03 , 7.970070455891505e+02), 
             std::complex<double> (1.930566204173797e+03 , -7.158771784887764e+02), 
             std::complex<double> (9.605158103310569e+02 , -9.751626584059218e+02), 
             std::complex<double> (4.154711000213336e+02 , -7.752252657037695e+02), 
             std::complex<double> (1.698139353899212e+02 , -5.306775644443155e+02), 
             std::complex<double> (6.785740977726294e+01 , -3.443282231993282e+02), 
             std::complex<double> (2.687443406079781e+01 , -2.189017173941930e+02), 
             std::complex<double> (1.060595450100369e+01 , -1.380614226705834e+02), 
             std::complex<double> (4.179811117901654e+00 , -8.680178347784999e+01), 
             std::complex<double> (1.646364454077571e+00 , -5.450210718208769e+01), 
             std::complex<double> (6.483397505023826e-01 , -3.419679316545660e+01), 
             std::complex<double> (2.552986571074330e-01 , -2.143932211548994e+01), 
             std::complex<double> (1.005350552224406e-01 , -1.341924503247874e+01), 
             std::complex<double> (3.960856720787807e-02 , -8.365485017306149e+00), 
             std::complex<double> (1.561944493023762e-02 , -5.163119900059615e+00), 
             std::complex<double> (5.877642838882013e-03 , -3.150366045562627e+00), 
             std::complex<double> (1.225777886228267e-03 , -2.169601257628324e+00), 
             std::complex<double> (3.650367791888589e-05 , -2.003632537865038e+00), 
             std::complex<double> (4.039716013375616e-08 , -2.000003522523785e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 20;
            int m = 2;
            double b = 1.510880363666204e+03;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.927141208542207e+03 , 5.106538085828216e+02), 
             std::complex<double> (-2.325981028895545e+03 , 1.267125509004393e+03), 
             std::complex<double> (-1.535688008522509e+03 , 1.507944389415590e+03), 
             std::complex<double> (-8.950254109695736e+02 , 1.378003992149494e+03), 
             std::complex<double> (-4.844699425090780e+02 , 1.107700097778699e+03), 
             std::complex<double> (-2.518230561807228e+02 , 8.346384901737318e+02), 
             std::complex<double> (-1.281252156642089e+02 , 6.087211109238522e+02), 
             std::complex<double> (-6.445609392413515e+01 , 4.367766010645107e+02), 
             std::complex<double> (-3.221675285006782e+01 , 3.108944876292850e+02), 
             std::complex<double> (-1.602547749921077e+01 , 2.204618540587877e+02), 
             std::complex<double> (-7.927018945924983e+00 , 1.561190466467564e+02), 
             std::complex<double> (-3.884381939842143e+00 , 1.105899675614983e+02), 
             std::complex<double> (-1.867880188187992e+00 , 7.850932135841815e+01), 
             std::complex<double> (-8.615654932368548e-01 , 5.602289474000054e+01), 
             std::complex<double> (-3.582812215974934e-01 , 4.040488905038147e+01), 
             std::complex<double> (-1.106217794879829e-01 , 2.970136723593025e+01), 
             std::complex<double> (-1.512926250933984e-02 , 2.213131238487415e+01), 
             std::complex<double> (-4.354708159150232e-04 , 1.571106968396146e+01), 
             std::complex<double> (-1.546488206167689e-06 , 9.424787610884120e+00), 
             std::complex<double> (-2.653887044152492e-10 , 3.141592655157290e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.111834864143358e+02 , 3.089808067753514e+02), 
             std::complex<double> (2.461187480730546e+02 , 1.588750617806762e+02), 
             std::complex<double> (2.382290027542898e+02 , 4.370527897657581e+00), 
             std::complex<double> (1.663352309506400e+02 , -7.399339617669277e+01), 
             std::complex<double> (9.842684359488585e+01 , -9.094989637795143e+01), 
             std::complex<double> (5.348939252118156e+01 , -8.050606806339782e+01), 
             std::complex<double> (2.783550262120055e+01 , -6.309981344466345e+01), 
             std::complex<double> (1.417228302507842e+01 , -4.684001710946127e+01), 
             std::complex<double> (7.137129022089970e+00 , -3.387951906643313e+01), 
             std::complex<double> (3.574647928315040e+00 , -2.419122532916500e+01), 
             std::complex<double> (1.785548941189774e+00 , -1.715527650070549e+01), 
             std::complex<double> (8.907943095958681e-01 , -1.211180272374532e+01), 
             std::complex<double> (4.443369402086432e-01 , -8.513795443379966e+00), 
             std::complex<double> (2.219619650623886e-01 , -5.944411863209304e+00), 
             std::complex<double> (1.110522843950140e-01 , -4.099195632842489e+00), 
             std::complex<double> (5.151317985913287e-02 , -2.803878435380058e+00), 
             std::complex<double> (1.277456133110521e-02 , -2.132757003019982e+00), 
             std::complex<double> (6.266086797201381e-04 , -2.004658545693135e+00), 
             std::complex<double> (3.387774389032928e-06 , -2.000021408739579e+00), 
             std::complex<double> (9.336869462210888e-10 , -2.000000005527582e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 3.867378956616377e+10;
             std::complex<double> zvec1[] = {
             std::complex<double> (-6.077320369364378e+10 , 3.059678380237500e+10), 
             std::complex<double> (-1.460591138411945e+10 , 2.998846283807779e+10), 
             std::complex<double> (-2.187022290547302e+09 , 1.272079536885996e+10), 
             std::complex<double> (-3.019947565476583e+08 , 4.786854771686548e+09), 
             std::complex<double> (-4.122599256073815e+07 , 1.771663032000535e+09), 
             std::complex<double> (-5.619031450880168e+06 , 6.542260379070795e+08), 
             std::complex<double> (-7.657004317265709e+05 , 2.415129232401696e+08), 
             std::complex<double> (-1.043382765859096e+05 , 8.915272935196063e+07), 
             std::complex<double> (-1.421761338798866e+04 , 3.290988926955791e+07), 
             std::complex<double> (-1.937356333725410e+03 , 1.214836576612728e+07), 
             std::complex<double> (-2.639929209504956e+02 , 4.484450672581832e+06), 
             std::complex<double> (-3.597286697153279e+01 , 1.655391180191916e+06), 
             std::complex<double> (-4.901825221181357e+00 , 6.110714899422112e+05), 
             std::complex<double> (-6.679448315017542e-01 , 2.255710736027956e+05), 
             std::complex<double> (-9.101717711898168e-02 , 8.326735923839128e+04), 
             std::complex<double> (-1.240240937633931e-02 , 3.073733272919445e+04), 
             std::complex<double> (-1.690006930054038e-03 , 1.134638663765282e+04), 
             std::complex<double> (-2.302872728253521e-04 , 4.188408922003046e+03), 
             std::complex<double> (-3.137932262441905e-05 , 1.546113204863072e+03), 
             std::complex<double> (-4.275322047103569e-06 , 5.707405828619054e+02), 
             std::complex<double> (-5.820414759643836e-07 , 2.107049456895719e+02), 
             std::complex<double> (-7.878369771534275e-08 , 7.783799434035748e+01), 
             std::complex<double> (-1.021345387634282e-08 , 2.888503208043884e+01), 
             std::complex<double> (-9.547142584879242e-10 , 1.092131801702088e+01), 
             std::complex<double> (-1.377894061757376e-11 , 3.159754212802467e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.692123827916010e+10 , 1.254542897972015e+10), 
             std::complex<double> (8.283262825530298e+09 , -6.486299976716513e+09), 
             std::complex<double> (1.366304795103021e+09 , -3.856098788856364e+09), 
             std::complex<double> (1.911953096040010e+08 , -1.509267032556283e+09), 
             std::complex<double> (2.614798773212380e+07 , -5.615430484522009e+08), 
             std::complex<double> (3.564808418811441e+06 , -2.075107131353775e+08), 
             std::complex<double> (4.857896789660904e+05 , -7.661175394573624e+07), 
             std::complex<double> (6.619650313517894e+04 , -2.828104430008595e+07), 
             std::complex<double> (9.020245589203836e+03 , -1.043969968781972e+07), 
             std::complex<double> (1.229139584856770e+03 , -3.853714624625777e+06), 
             std::complex<double> (1.674881111582609e+02 , -1.422561194960687e+06), 
             std::complex<double> (2.282268602571365e+01 , -5.251245778205775e+05), 
             std::complex<double> (3.109922208626706e+00 , -1.938446103894729e+05), 
             std::complex<double> (4.237720376404838e-01 , -7.155584509053650e+04), 
             std::complex<double> (5.774509134515964e-02 , -2.641414147579403e+04), 
             std::complex<double> (7.868606908817414e-03 , -9.750522352166572e+03), 
             std::complex<double> (1.072211537678262e-03 , -3.599309987712317e+03), 
             std::complex<double> (1.461043361598581e-04 , -1.328649790290224e+03), 
             std::complex<double> (1.990881361665280e-05 , -4.904571604462841e+02), 
             std::complex<double> (2.712842261234037e-06 , -1.810449373218616e+02), 
             std::complex<double> (3.696608362272592e-07 , -6.682411872254279e+01), 
             std::complex<double> (5.037080204457085e-08 , -2.464917947356721e+01), 
             std::complex<double> (6.858622512117156e-09 , -9.055590449803987e+00), 
             std::complex<double> (8.653130651466150e-10 , -3.375077256324674e+00), 
             std::complex<double> (2.855449324229107e-11 , -2.038736577872309e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 1.767980261388775e+09;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.932655679354784e+09 , 1.284247920399536e+09), 
             std::complex<double> (-9.259043230571325e+08 , 1.542328752657431e+09), 
             std::complex<double> (-1.847767716302595e+08 , 7.820873231726789e+08), 
             std::complex<double> (-3.326166069214310e+07 , 3.393297005108014e+08), 
             std::complex<double> (-5.874197308565493e+06 , 1.431645022532985e+08), 
             std::complex<double> (-1.033905861368113e+06 , 6.010389351872045e+07), 
             std::complex<double> (-1.818670887732783e+05 , 2.521111837568082e+07), 
             std::complex<double> (-3.198759615877932e+04 , 1.057341233012433e+07), 
             std::complex<double> (-5.626017366237886e+03 , 4.434314977046109e+06), 
             std::complex<double> (-9.895075693209335e+02 , 1.859669861954770e+06), 
             std::complex<double> (-1.740351310589180e+02 , 7.799105689244438e+05), 
             std::complex<double> (-3.060939068999734e+01 , 3.270797784793235e+05), 
             std::complex<double> (-5.383595692458747e+00 , 1.371710847732878e+05), 
             std::complex<double> (-9.468696199055263e-01 , 5.752696339250897e+04), 
             std::complex<double> (-1.665359069650703e-01 , 2.412572263571203e+04), 
             std::complex<double> (-2.929040576991615e-02 , 1.011787323094865e+04), 
             std::complex<double> (-5.151596472610081e-03 , 4.243246538712087e+03), 
             std::complex<double> (-9.060506689300735e-04 , 1.779541127851121e+03), 
             std::complex<double> (-1.593420252114536e-04 , 7.463146315166932e+02), 
             std::complex<double> (-2.801054687150691e-05 , 3.130109576506529e+02), 
             std::complex<double> (-4.911889095193575e-06 , 1.313201614183852e+02), 
             std::complex<double> (-8.492918922036244e-07 , 5.518980769577839e+01), 
             std::complex<double> (-1.349971143507167e-07 , 2.340886474914185e+01), 
             std::complex<double> (-1.288850676466663e-08 , 1.015828114798678e+01), 
             std::complex<double> (-1.249074556375340e-10 , 3.147253059075107e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (6.381758991990042e+08 , 5.889233575065727e+08), 
             std::complex<double> (4.364889293543850e+08 , -2.325231407979714e+08), 
             std::complex<double> (9.936060296953343e+07 , -1.985396356749212e+08), 
             std::complex<double> (1.830882242995034e+07 , -9.249438682146960e+07), 
             std::complex<double> (3.246788914719936e+06 , -3.949836359824469e+07), 
             std::complex<double> (5.718752561212354e+05 , -1.661745031779712e+07), 
             std::complex<double> (1.006073689724109e+05 , -6.972928181302273e+06), 
             std::complex<double> (1.769567036288116e+04 , -2.924600926924957e+06), 
             std::complex<double> (3.112348327917191e+03 , -1.226543546410261e+06), 
             std::complex<double> (5.474022145395021e+02 , -5.143907655510101e+05), 
             std::complex<double> (9.627741187995788e+01 , -2.157258905768880e+05), 
             std::complex<double> (1.693332237101915e+01 , -9.047137417366235e+04), 
             std::complex<double> (2.978241635773591e+00 , -3.794198688127947e+04), 
             std::complex<double> (5.238147035808129e-01 , -1.591215298501723e+04), 
             std::complex<double> (9.212880486985725e-02 , -6.673256476615108e+03), 
             std::complex<double> (1.620366190585845e-02 , -2.798637601312119e+03), 
             std::complex<double> (2.849908244862314e-03 , -1.173695350699257e+03), 
             std::complex<double> (5.012432452484116e-04 , -4.922246643761710e+02), 
             std::complex<double> (8.815890249401233e-05 , -2.064273307158491e+02), 
             std::complex<double> (1.550543416014607e-05 , -8.656600503326850e+01), 
             std::complex<double> (2.727106346110066e-06 , -3.629056062550842e+01), 
             std::complex<double> (4.796333355889154e-07 , -1.518793607690816e+01), 
             std::complex<double> (8.411667129933547e-08 , -6.308038107375014e+00), 
             std::complex<double> (1.232839497736220e-08 , -2.743980107230492e+00), 
             std::complex<double> (2.823699846137673e-10 , -2.013099043213079e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 1.652235747470210e+08;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.846723006826614e+08 , 1.103262119615021e+08), 
             std::complex<double> (-1.101745521898821e+08 , 1.547126323552420e+08), 
             std::complex<double> (-2.759500506014759e+07 , 9.096121385542513e+07), 
             std::complex<double> (-6.124478582229281e+06 , 4.436020178283665e+07), 
             std::complex<double> (-1.322174627952584e+06 , 2.076466267252767e+07), 
             std::complex<double> (-2.837234202331722e+05 , 9.634261835367033e+06), 
             std::complex<double> (-6.080508321988305e+04 , 4.461580259382267e+06), 
             std::complex<double> (-1.302759281475449e+04 , 2.065297138332049e+06), 
             std::complex<double> (-2.791018174257961e+03 , 9.559574815322031e+05), 
             std::complex<double> (-5.979372718073802e+02 , 4.424727314750738e+05), 
             std::complex<double> (-1.280994903758838e+02 , 2.048013010982928e+05), 
             std::complex<double> (-2.744346359050791e+01 , 9.479349591996644e+04), 
             std::complex<double> (-5.879364366324149e+00 , 4.387572346506749e+04), 
             std::complex<double> (-1.259568422780648e+00 , 2.030813450210310e+04), 
             std::complex<double> (-2.698440906799862e-01 , 9.399739127143552e+03), 
             std::complex<double> (-5.780999785020665e-02 , 4.350725892187787e+03), 
             std::complex<double> (-1.238476568295157e-02 , 2.013763002010345e+03), 
             std::complex<double> (-2.653066622690015e-03 , 9.320910185761129e+02), 
             std::complex<double> (-5.681904106493250e-04 , 4.314434038509463e+02), 
             std::complex<double> (-1.215356407006027e-04 , 1.997385331683443e+02), 
             std::complex<double> (-2.584618280724964e-05 , 9.254166528875541e+01), 
             std::complex<double> (-5.346081594905893e-06 , 4.302918264823252e+01), 
             std::complex<double> (-9.593790903523458e-07 , 2.030895866732791e+01), 
             std::complex<double> (-8.363793627475115e-08 , 9.774042277543659e+00), 
             std::complex<double> (-5.174431041864318e-10 , 3.143296733062278e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (4.964920095470116e+07 , 5.443346978302167e+07), 
             std::complex<double> (4.374305791163472e+07 , -1.513786104680932e+07), 
             std::complex<double> (1.293038487943203e+07 , -1.934981401840410e+07), 
             std::complex<double> (2.974296942773023e+06 , -1.056624352272299e+07), 
             std::complex<double> (6.470544838033246e+05 , -5.060372986848258e+06), 
             std::complex<double> (1.390793096384244e+05 , -2.359277228624903e+06), 
             std::complex<double> (2.981676289467578e+04 , -1.093701092857774e+06), 
             std::complex<double> (6.388775336217243e+03 , -5.063941820539886e+05), 
             std::complex<double> (1.368746870901357e+03 , -2.344041809567730e+05), 
             std::complex<double> (2.932362166060562e+02 , -1.084969936728278e+05), 
             std::complex<double> (6.282170332442647e+01 , -5.021862881884091e+04), 
             std::complex<double> (1.345864382797343e+01 , -2.324400076133422e+04), 
             std::complex<double> (2.883319596899132e+00 , -1.075862254797321e+04), 
             std::complex<double> (6.177094454714178e-01 , -4.979691095708837e+03), 
             std::complex<double> (1.323352954808784e-01 , -2.304878789278241e+03), 
             std::complex<double> (2.835091912773195e-02 , -1.066826080659305e+03), 
             std::complex<double> (6.073773969700597e-03 , -4.937856256796169e+02), 
             std::complex<double> (1.301218416290769e-03 , -2.285493189900311e+02), 
             std::complex<double> (2.787673121752139e-04 , -1.057805639168893e+02), 
             std::complex<double> (5.972191203728079e-05 , -4.895075262854895e+01), 
             std::complex<double> (1.279459692542021e-05 , -2.263481099893045e+01), 
             std::complex<double> (2.740932652981146e-06 , -1.043002305859359e+01), 
             std::complex<double> (5.816028935902000e-07 , -4.754756830176578e+00), 
             std::complex<double> (8.717280879770618e-08 , -2.392986873161238e+00), 
             std::complex<double> (1.259211474614584e-09 , -2.004221935414873e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 2.532690578681951e+07;
             std::complex<double> zvec1[] = {
             std::complex<double> (-4.485855251667663e+07 , 1.559434905817375e+07), 
             std::complex<double> (-2.028458634865906e+07 , 2.466653547890673e+07), 
             std::complex<double> (-6.104811399024409e+06 , 1.642191258320905e+07), 
             std::complex<double> (-1.606919213769642e+06 , 8.843821967464814e+06), 
             std::complex<double> (-4.078531531727845e+05 , 4.510031034102946e+06), 
             std::complex<double> (-1.025568407292848e+05 , 2.268481154538866e+06), 
             std::complex<double> (-2.572793816316750e+04 , 1.137070996325078e+06), 
             std::complex<double> (-6.450437731618684e+03 , 5.694596620402246e+05), 
             std::complex<double> (-1.616996804271129e+03 , 2.851306011235656e+05), 
             std::complex<double> (-4.053339989421064e+02 , 1.427581789642742e+05), 
             std::complex<double> (-1.016044815887522e+02 , 7.147468177687380e+04), 
             std::complex<double> (-2.546898656996068e+01 , 3.578507859020046e+04), 
             std::complex<double> (-6.384253677279625e+00 , 1.791642507496807e+04), 
             std::complex<double> (-1.600325198837434e+00 , 8.970170319483859e+03), 
             std::complex<double> (-4.011484683363307e-01 , 4.491073548820071e+03), 
             std::complex<double> (-1.005535105564264e-01 , 2.248538451018544e+03), 
             std::complex<double> (-2.520404389470724e-02 , 1.125779086195757e+03), 
             std::complex<double> (-6.316361410832303e-03 , 5.636595074497301e+02), 
             std::complex<double> (-1.581827818295049e-03 , 2.822432614005124e+02), 
             std::complex<double> (-3.950322996892122e-04 , 1.413845708842510e+02), 
             std::complex<double> (-9.753921979638200e-05 , 7.093527042948136e+01), 
             std::complex<double> (-2.296503065539771e-05 , 3.580815253088542e+01), 
             std::complex<double> (-4.344331102989828e-06 , 1.845658665179792e+01), 
             std::complex<double> (-3.161894954300827e-07 , 9.584244111689816e+00), 
             std::complex<double> (-1.234329793473086e-09 , 3.142091104567090e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (6.401851326497968e+06 , 8.095001615585194e+06), 
             std::complex<double> (6.859297924436153e+06 , -1.350153746361203e+06), 
             std::complex<double> (2.515696021573704e+06 , -2.916002603762225e+06), 
             std::complex<double> (6.959952505782511e+05 , -1.852004657654446e+06), 
             std::complex<double> (1.788762416841343e+05 , -9.809166290685949e+05), 
             std::complex<double> (4.512086164736761e+04 , -4.980000304495365e+05), 
             std::complex<double> (1.132816840318087e+04 , -2.502014959046308e+05), 
             std::complex<double> (2.840727818641510e+03 , -1.253769441634810e+05), 
             std::complex<double> (7.121493584734540e+02 , -6.278585196461254e+04), 
             std::complex<double> (1.785173153859296e+02 , -3.143654533669901e+04), 
             std::complex<double> (4.474881452506213e+01 , -1.573946684005197e+04), 
             std::complex<double> (1.121710276388816e+01 , -7.880263701584214e+03), 
             std::complex<double> (2.811767162341049e+00 , -3.945393768471392e+03), 
             std::complex<double> (7.048194552619181e-01 , -1.975329886139517e+03), 
             std::complex<double> (1.766755228405779e-01 , -9.889826534762695e+02), 
             std::complex<double> (4.428685900994399e-02 , -4.951502668855513e+02), 
             std::complex<double> (1.110128792356559e-02 , -2.479034941959124e+02), 
             std::complex<double> (2.782735147080057e-03 , -1.241130577412773e+02), 
             std::complex<double> (6.975421393125956e-04 , -6.213112608684230e+01), 
             std::complex<double> (1.748516072984516e-04 , -3.109062229748772e+01), 
             std::complex<double> (4.383054944900958e-05 , -1.553354716832274e+01), 
             std::complex<double> (1.098585293879195e-05 , -7.715176851384284e+00), 
             std::complex<double> (2.677994895564640e-06 , -3.790763970489791e+00), 
             std::complex<double> (3.647200633727070e-07 , -2.198773295698555e+00), 
             std::complex<double> (3.194386060452414e-09 , -2.001307672244279e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 5.544898076864207e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.002609980890299e+07 , 3.160209167152561e+06), 
             std::complex<double> (-5.118501562828318e+06 , 5.497094676223950e+06), 
             std::complex<double> (-1.791213622776857e+06 , 4.066282654151815e+06), 
             std::complex<double> (-5.435137692224783e+05 , 2.386484029672952e+06), 
             std::complex<double> (-1.577144567287188e+05 , 1.309001607657194e+06), 
             std::complex<double> (-4.516940853951985e+04 , 7.041488679300844e+05), 
             std::complex<double> (-1.288793394445529e+04 , 3.766791560520380e+05), 
             std::complex<double> (-3.673292998230233e+03 , 2.011821279192351e+05), 
             std::complex<double> (-1.046633925383539e+03 , 1.074016062626152e+05), 
             std::complex<double> (-2.981920970195383e+02 , 5.732924080867106e+04), 
             std::complex<double> (-8.495454567140933e+01 , 3.060029856596788e+04), 
             std::complex<double> (-2.420326393147524e+01 , 1.633317477328620e+04), 
             std::complex<double> (-6.895408857120517e+00 , 8.717948761504394e+03), 
             std::complex<double> (-1.964466472173925e+00 , 4.653265450619026e+03), 
             std::complex<double> (-5.596606841221657e-01 , 2.483715664955833e+03), 
             std::complex<double> (-1.594372003932122e-01 , 1.325708780904054e+03), 
             std::complex<double> (-4.541515286271827e-02 , 7.076234812664014e+02), 
             std::complex<double> (-1.293073336466620e-02 , 3.777321098518221e+02), 
             std::complex<double> (-3.676052262374842e-03 , 2.016797505364832e+02), 
             std::complex<double> (-1.039423624394977e-03 , 1.077654749677000e+02), 
             std::complex<double> (-2.882349530758460e-04 , 5.774033972624650e+01), 
             std::complex<double> (-7.419907923414202e-05 , 3.122497845084447e+01), 
             std::complex<double> (-1.391940971779196e-05 , 1.732451124662663e+01), 
             std::complex<double> (-7.878771461676014e-07 , 9.494141891643114e+00), 
             std::complex<double> (-1.968299363710694e-09 , 3.141735081105699e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.191295470326250e+06 , 1.702008686015660e+06), 
             std::complex<double> (1.488761670258951e+06 , -1.063252048076552e+05), 
             std::complex<double> (6.538718188084075e+05 , -5.981695165057027e+05), 
             std::complex<double> (2.116812714528147e+05 , -4.406248411379426e+05), 
             std::complex<double> (6.257183716366295e+04 , -2.558980243186552e+05), 
             std::complex<double> (1.801537927643874e+04 , -1.398436426457292e+05), 
             std::complex<double> (5.147962518016867e+03 , -7.514238353250042e+04), 
             std::complex<double> (1.467890702366131e+03 , -4.018396083541004e+04), 
             std::complex<double> (4.182981776172394e+02 , -2.146004827387551e+04), 
             std::complex<double> (1.191797331436746e+02 , -1.145620308633133e+04), 
             std::complex<double> (3.395449253052530e+01 , -6.115090245261643e+03), 
             std::complex<double> (9.673550997198095e+00 , -3.264009351838906e+03), 
             std::complex<double> (2.755959642223176e+00 , -1.742191874164765e+03), 
             std::complex<double> (7.851620173567884e-01 , -9.299065042716626e+02), 
             std::complex<double> (2.236894886488015e-01 , -4.963426948200848e+02), 
             std::complex<double> (6.372822832602859e-02 , -2.649242145839163e+02), 
             std::complex<double> (1.815591390779050e-02 , -1.414014281589783e+02), 
             std::complex<double> (5.172546489847683e-03 , -7.546722122474711e+01), 
             std::complex<double> (1.473638615976540e-03 , -4.026856729391694e+01), 
             std::complex<double> (4.198379092446877e-04 , -2.147014017536402e+01), 
             std::complex<double> (1.196216715652353e-04 , -1.141623593004148e+01), 
             std::complex<double> (3.406515196956009e-05 , -6.017089152924379e+00), 
             std::complex<double> (9.055559978888164e-06 , -3.158204500196164e+00), 
             std::complex<double> (1.008667249298538e-06 , -2.095194761862561e+00), 
             std::complex<double> (5.364430684081745e-09 , -2.000392276848490e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 1.582461574498143e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.907075804495280e+06 , 8.381603087091710e+05), 
             std::complex<double> (-1.634001450709266e+06 , 1.573236096829913e+06), 
             std::complex<double> (-6.481945111191381e+05 , 1.273120040012515e+06), 
             std::complex<double> (-2.218708611921907e+05 , 8.058443946504710e+05), 
             std::complex<double> (-7.209133440994893e+04 , 4.709555467668507e+05), 
             std::complex<double> (-2.302688958308514e+04 , 2.682820629772965e+05), 
             std::complex<double> (-7.314845380342371e+03 , 1.515881939634075e+05), 
             std::complex<double> (-2.319622523760162e+03 , 8.543119507728654e+04), 
             std::complex<double> (-7.351722443070878e+02 , 4.810738041933499e+04), 
             std::complex<double> (-2.329617693710550e+02 , 2.708283752806788e+04), 
             std::complex<double> (-7.381692305476086e+01 , 1.524547298648565e+04), 
             std::complex<double> (-2.338940469058382e+01 , 8.581760254172368e+03), 
             std::complex<double> (-7.411032334520423e+00 , 4.830682043602524e+03), 
             std::complex<double> (-2.348191402667750e+00 , 2.719192705117721e+03), 
             std::complex<double> (-7.440043386251874e-01 , 1.530639923271575e+03), 
             std::complex<double> (-2.357098428814943e-01 , 8.616123131971789e+02), 
             std::complex<double> (-7.465429905341259e-02 , 4.850309094440976e+02), 
             std::complex<double> (-2.362308303327609e-02 , 2.730774972363135e+02), 
             std::complex<double> (-7.453581931343394e-03 , 1.538113815158709e+02), 
             std::complex<double> (-2.330131387467641e-03 , 8.675151319314284e+01), 
             std::complex<double> (-7.065596889720038e-04 , 4.913607539352810e+01), 
             std::complex<double> (-1.920346812872231e-04 , 2.818825776481373e+01), 
             std::complex<double> (-3.375259119311689e-05 , 1.663323083821144e+01), 
             std::complex<double> (-1.409585250856134e-06 , 9.453485964587877e+00), 
             std::complex<double> (-2.318360528765509e-09 , 3.141632582666794e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (2.918468413235057e+05 , 4.640486451952496e+05), 
             std::complex<double> (4.124144637368025e+05 , 1.563308582704653e+04), 
             std::complex<double> (2.108774597166502e+05 , -1.534094808921258e+05), 
             std::complex<double> (7.820053129700478e+04 , -1.312484603393439e+05), 
             std::complex<double> (2.606395510247544e+04 , -8.313994584891268e+04), 
             std::complex<double> (8.392610243572997e+03 , -4.853020778726826e+04), 
             std::complex<double> (2.672869743497869e+03 , -2.763093504184038e+04), 
             std::complex<double> (8.482854084470169e+02 , -1.560956719215150e+04), 
             std::complex<double> (2.689213901738477e+02 , -8.796633267121466e+03), 
             std::complex<double> (8.522289937459077e+01 , -4.953401931742339e+03), 
             std::complex<double> (2.700467611974488e+01 , -2.788581904704634e+03), 
             std::complex<double> (8.556700659832833e+00 , -1.569745363752446e+03), 
             std::complex<double> (2.711245565153868e+00 , -8.836169965353700e+02), 
             std::complex<double> (8.590725216115385e-01 , -4.973875470068554e+02), 
             std::complex<double> (2.722013705589704e-01 , -2.799772843202516e+02), 
             std::complex<double> (8.624831424832065e-02 , -1.575957201401409e+02), 
             std::complex<double> (2.732818868427299e-02 , -8.870482331686209e+01), 
             std::complex<double> (8.659068302080705e-03 , -4.992188767746226e+01), 
             std::complex<double> (2.743678002597711e-03 , -2.808331858458761e+01), 
             std::complex<double> (8.693829062662288e-04 , -1.577675282365806e+01), 
             std::complex<double> (2.755491967777021e-04 , -8.825588828273037e+00), 
             std::complex<double> (8.712583276945072e-05 , -4.879445916639436e+00), 
             std::complex<double> (2.382112992949155e-05 , -2.733110063500346e+00), 
             std::complex<double> (1.991169224268951e-06 , -2.042980526108393e+00), 
             std::complex<double> (6.604165683614133e-09 , -2.000114689334854e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 5.529312911523790e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.028457075629534e+06 , 2.732333550270321e+05), 
             std::complex<double> (-6.246229036149400e+05 , 5.455143588873432e+05), 
             std::complex<double> (-2.752861261270817e+05 , 4.767816465223085e+05), 
             std::complex<double> (-1.044963637691675e+05 , 3.227011061322295e+05), 
             std::complex<double> (-3.741725402929212e+04 , 1.994947001897425e+05), 
             std::complex<double> (-1.311770258749830e+04 , 1.194618935099800e+05), 
             std::complex<double> (-4.564670354209729e+03 , 7.074662103802926e+04), 
             std::complex<double> (-1.584286960719915e+03 , 4.173566801395968e+04), 
             std::complex<double> (-5.493723380910128e+02 , 2.458828259375870e+04), 
             std::complex<double> (-1.904424740662533e+02 , 1.447930359221225e+04), 
             std::complex<double> (-6.601054802528526e+01 , 8.525059485958051e+03), 
             std::complex<double> (-2.287943459110757e+01 , 5.019068928868332e+03), 
             std::complex<double> (-7.929904384581426e+00 , 2.954888306892509e+03), 
             std::complex<double> (-2.748388835547903e+00 , 1.739633121128664e+03), 
             std::complex<double> (-9.524835500011295e-01 , 1.024183703311025e+03), 
             std::complex<double> (-3.300268534112832e-01 , 6.029912142102426e+02), 
             std::complex<double> (-1.142849079392104e-01 , 3.550440967071883e+02), 
             std::complex<double> (-3.950925499188210e-02 , 2.091048331996693e+02), 
             std::complex<double> (-1.359205755669124e-02 , 1.232436185417129e+02), 
             std::complex<double> (-4.608809447029173e-03 , 7.279170914924730e+01), 
             std::complex<double> (-1.494287338736338e-03 , 4.325330552367545e+01), 
             std::complex<double> (-4.150752378264587e-04 , 2.612498895502299e+01), 
             std::complex<double> (-6.480508493301116e-05 , 1.621892036804682e+01), 
             std::complex<double> (-1.930211477569133e-06 , 9.436111337543549e+00), 
             std::complex<double> (-2.161926391177454e-09 , 3.141603672944558e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (8.834185860184884e+04 , 1.545254388405145e+05), 
             std::complex<double> (1.378967644343935e+05 , 1.873096685692346e+04), 
             std::complex<double> (8.023313843775926e+04 , -4.631720082688536e+04), 
             std::complex<double> (3.350263245523042e+04 , -4.630630848492371e+04), 
             std::complex<double> (1.240009051827291e+04 , -3.189342364562411e+04), 
             std::complex<double> (4.397430173681727e+03 , -1.978209290945825e+04), 
             std::complex<double> (1.536315230053567e+03 , -1.185590895556145e+04), 
             std::complex<double> (5.339547064627275e+02 , -7.022983247186423e+03), 
             std::complex<double> (1.852444829473326e+02 , -4.143426839519464e+03), 
             std::complex<double> (6.422654810531527e+01 , -2.441140859295768e+03), 
             std::complex<double> (2.226331046526604e+01 , -1.437528444693137e+03), 
             std::complex<double> (7.716711670362090e+00 , -8.463836611376719e+02), 
             std::complex<double> (2.674628657586064e+00 , -4.983017178630948e+02), 
             std::complex<double> (9.270236104038946e-01 , -2.933642462223322e+02), 
             std::complex<double> (3.213044496598842e-01 , -1.727087563285458e+02), 
             std::complex<double> (1.113633368968346e-01 , -1.016733660592169e+02), 
             std::complex<double> (3.859825603485854e-02 , -5.984961410423781e+01), 
             std::complex<double> (1.337808218583708e-02 , -3.522124589742737e+01), 
             std::complex<double> (4.636885710842437e-03 , -2.071229067722778e+01), 
             std::complex<double> (1.607337056359469e-03 , -1.215416610136053e+01), 
             std::complex<double> (5.574786128188092e-04 , -7.088712673039149e+00), 
             std::complex<double> (1.917176710114599e-04 , -4.077236896870208e+00), 
             std::complex<double> (5.048304262468848e-05 , -2.448381222044735e+00), 
             std::complex<double> (2.977085698122081e-06 , -2.018319335412383e+00), 
             std::complex<double> (6.399651460136722e-09 , -2.000032837341610e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 2.262579246114125e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-4.250278840434309e+05 , 1.046991821961638e+05), 
             std::complex<double> (-2.749692977014305e+05 , 2.198937319220287e+05), 
             std::complex<double> (-1.324816420703487e+05 , 2.053627241878025e+05), 
             std::complex<double> (-5.500604547661055e+04 , 1.475537274352964e+05), 
             std::complex<double> (-2.143354576121021e+04 , 9.593358031575846e+04), 
             std::complex<double> (-8.145192267487366e+03 , 6.004717384277641e+04), 
             std::complex<double> (-3.065878231628334e+03 , 3.705073070384436e+04), 
             std::complex<double> (-1.149851975217563e+03 , 2.273880295524879e+04), 
             std::complex<double> (-4.306662823203903e+02 , 1.392721470208671e+04), 
             std::complex<double> (-1.612200603286548e+02 , 8.523809485780688e+03), 
             std::complex<double> (-6.034114922252016e+01 , 5.215319157728231e+03), 
             std::complex<double> (-2.258259667377936e+01 , 3.190676062332682e+03), 
             std::complex<double> (-8.451109563895010e+00 , 1.951950280404210e+03), 
             std::complex<double> (-3.162463355414173e+00 , 1.194131162884469e+03), 
             std::complex<double> (-1.183238066346120e+00 , 7.305378177479554e+02), 
             std::complex<double> (-4.425356580765495e-01 , 4.469496912064553e+02), 
             std::complex<double> (-1.653367020516970e-01 , 2.734917055359470e+02), 
             std::complex<double> (-6.159816695603971e-02 , 1.674234687829023e+02), 
             std::complex<double> (-2.277458586423954e-02 , 1.026094240217175e+02), 
             std::complex<double> (-8.243505325082634e-03 , 6.307942377555342e+01), 
             std::complex<double> (-2.801784616164576e-03 , 3.909277773998902e+01), 
             std::complex<double> (-7.702577450462051e-04 , 2.470837545197407e+01), 
             std::complex<double> (-1.016831791651275e-04 , 1.597850209907064e+01), 
             std::complex<double> (-2.126044768331737e-06 , 9.429065404364550e+00), 
             std::complex<double> (-1.676847370376560e-09 , 3.141595654835235e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (3.158012519779156e+04 , 6.021035094222407e+04), 
             std::complex<double> (5.348223949869934e+04 , 1.205387040242063e+04), 
             std::complex<double> (3.475255070549959e+04 , -1.572573157493599e+04), 
             std::complex<double> (1.611024014908957e+04 , -1.860502037792518e+04), 
             std::complex<double> (6.541935250285595e+03 , -1.390959549700811e+04), 
             std::complex<double> (2.524803124283017e+03 , -9.135309727606527e+03), 
             std::complex<double> (9.558617553771969e+02 , -5.736182094660154e+03), 
             std::complex<double> (3.592719005135438e+02 , -3.543290313118555e+03), 
             std::complex<double> (1.346711796639663e+02 , -2.175464519777428e+03), 
             std::complex<double> (5.042958889873542e+01 , -1.332640875597609e+03), 
             std::complex<double> (1.887692486229664e+01 , -8.156545332838579e+02), 
             std::complex<double> (7.065051487895920e+00 , -4.990697908073975e+02), 
             std::complex<double> (2.644090428874665e+00 , -3.053257688940086e+02), 
             std::complex<double> (9.895292385658513e-01 , -1.867853055062293e+02), 
             std::complex<double> (3.703204905796083e-01 , -1.142628232787177e+02), 
             std::complex<double> (1.385880245092619e-01 , -6.989374582520954e+01), 
             std::complex<double> (5.186490403733011e-02 , -4.274651011803243e+01), 
             std::complex<double> (1.940993003412860e-02 , -2.613215361251221e+01), 
             std::complex<double> (7.264284944960565e-03 , -1.595682199862547e+01), 
             std::complex<double> (2.719436843847922e-03 , -9.713170234848066e+00), 
             std::complex<double> (1.018948105117762e-03 , -5.863917929860847e+00), 
             std::complex<double> (3.726987936631255e-04 , -3.490959588854945e+00), 
             std::complex<double> (8.815876853901260e-05 , -2.262430555805313e+00), 
             std::complex<double> (3.538735057081614e-06 , -2.007407588959979e+00), 
             std::complex<double> (5.134642880935755e-09 , -2.000009240486565e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 1.049710280230338e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.987758394264543e+05 , 4.564419809295447e+04), 
             std::complex<double> (-1.354791919315292e+05 , 9.997675741657047e+04), 
             std::complex<double> (-7.043025946631936e+04 , 9.888175008457511e+04), 
             std::complex<double> (-3.162454821989473e+04 , 7.495028302932758e+04), 
             std::complex<double> (-1.327371560704787e+04 , 5.100417721974160e+04), 
             std::complex<double> (-5.413907310483050e+03 , 3.322029011027178e+04), 
             std::complex<double> (-2.182346466919670e+03 , 2.125814778972344e+04), 
             std::complex<double> (-8.755356163477635e+02 , 1.350724469410716e+04), 
             std::complex<double> (-3.505865208112791e+02 , 8.558040642074013e+03), 
             std::complex<double> (-1.402760814830147e+02 , 5.416121288105456e+03), 
             std::complex<double> (-5.610947706435270e+01 , 3.426142291814335e+03), 
             std::complex<double> (-2.244026990193558e+01 , 2.166928646346737e+03), 
             std::complex<double> (-8.973866202946223e+00 , 1.370424814824114e+03), 
             std::complex<double> (-3.588184836570763e+00 , 8.666837062014156e+02), 
             std::complex<double> (-1.434323093537390e+00 , 5.481250153406387e+02), 
             std::complex<double> (-5.729519807994330e-01 , 3.466914784102993e+02), 
             std::complex<double> (-2.284741091552447e-01 , 2.193422825966488e+02), 
             std::complex<double> (-9.071046018459382e-02 , 1.388648734448946e+02), 
             std::complex<double> (-3.561363968411716e-02 , 8.806231837551802e+01), 
             std::complex<double> (-1.357289977969387e-02 , 5.607940011105132e+01), 
             std::complex<double> (-4.747585681149294e-03 , 3.608119185272472e+01), 
             std::complex<double> (-1.249881330185700e-03 , 2.373799408458941e+01), 
             std::complex<double> (-1.335253116894805e-04 , 1.584477848305304e+01), 
             std::complex<double> (-1.959678836597898e-06 , 9.426340414736684e+00), 
             std::complex<double> (-1.121976751108328e-09 , 3.141593461907410e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.289761855145293e+04 , 2.660309123228703e+04), 
             std::complex<double> (2.336904946328959e+04 , 7.211194817318362e+03), 
             std::complex<double> (1.669880474023177e+04 , -5.775275663001779e+03), 
             std::complex<double> (8.491817123598381e+03 , -8.271296545240717e+03), 
             std::complex<double> (3.745978686507347e+03 , -6.709509100954338e+03), 
             std::complex<double> (1.558577228923486e+03 , -4.654792886904642e+03), 
             std::complex<double> (6.332868269242400e+02 , -3.051902086783011e+03), 
             std::complex<double> (2.548793282070340e+02 , -1.957799765977790e+03), 
             std::complex<double> (1.021904363963120e+02 , -1.245168595520624e+03), 
             std::complex<double> (4.090929850393462e+01 , -7.892251533517771e+02), 
             std::complex<double> (1.636696782933405e+01 , -4.995509498307842e+02), 
             std::complex<double> (6.546485417272803e+00 , -3.160242072720284e+02), 
             std::complex<double> (2.618217456092200e+00 , -1.998771448259560e+02), 
             std::complex<double> (1.047095426852841e+00 , -1.264039274095929e+02), 
             std::complex<double> (4.187550492392016e-01 , -7.993264903249630e+01), 
             std::complex<double> (1.674677929130239e-01 , -5.053999047797362e+01), 
             std::complex<double> (6.697344397118465e-02 , -3.194676559831449e+01), 
             std::complex<double> (2.678435552835270e-02 , -2.018015331077571e+01), 
             std::complex<double> (1.071287196743573e-02 , -1.272570342076817e+01), 
             std::complex<double> (4.287233251765175e-03 , -7.990246897928382e+00), 
             std::complex<double> (1.717175574411375e-03 , -4.964257722358515e+00), 
             std::complex<double> (6.496138615619556e-04 , -3.053742926962065e+00), 
             std::complex<double> (1.289543417402177e-04 , -2.145995594158197e+00), 
             std::complex<double> (3.482077749431756e-06 , -2.002860123765703e+00), 
             std::complex<double> (3.540894197991953e-09 , -2.000002562706204e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 5.389957740118304e+04;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.027391706743321e+05 , 2.209351949384508e+04), 
             std::complex<double> (-7.313579536359299e+04 , 5.012525996047722e+04), 
             std::complex<double> (-4.058515250369693e+04 , 5.210906364146554e+04), 
             std::complex<double> (-1.952197754199242e+04 , 4.144050304126611e+04), 
             std::complex<double> (-8.752276757708827e+03 , 2.939698228470431e+04), 
             std::complex<double> (-3.800571723523838e+03 , 1.985146042081549e+04), 
             std::complex<double> (-1.627494289324524e+03 , 1.312598239096218e+04), 
             std::complex<double> (-6.927700715173618e+02 , 8.601538346432042e+03), 
             std::complex<double> (-2.941368607868795e+02 , 5.615218831326496e+03), 
             std::complex<double> (-1.247487674279758e+02 , 3.659790063981623e+03), 
             std::complex<double> (-5.288308085242623e+01 , 2.383687443241005e+03), 
             std::complex<double> (-2.241284078970750e+01 , 1.552096987346185e+03), 
             std::complex<double> (-9.497388851328285e+00 , 1.010510416394529e+03), 
             std::complex<double> (-4.023546713750131e+00 , 6.578907970458954e+02), 
             std::complex<double> (-1.703731528204116e+00 , 4.283409589698261e+02), 
             std::complex<double> (-7.206134749316105e-01 , 2.789316027054076e+02), 
             std::complex<double> (-3.039785999502871e-01 , 1.817123282844690e+02), 
             std::complex<double> (-1.274105386236773e-01 , 1.184938493189041e+02), 
             std::complex<double> (-5.257432696984863e-02 , 7.744813580840386e+01), 
             std::complex<double> (-2.084014127549728e-02 , 5.089673808011091e+01), 
             std::complex<double> (-7.366966547123454e-03 , 3.386996174024080e+01), 
             std::complex<double> (-1.794771689024560e-03 , 2.308235626019172e+01), 
             std::complex<double> (-1.496019498419593e-04 , 1.577395512479946e+01), 
             std::complex<double> (-1.559461779785682e-06 , 9.425329234877381e+00), 
             std::complex<double> (-6.653414111055195e-10 , 3.141592869201194e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (5.870416131556249e+03 , 1.301808911719944e+04), 
             std::complex<double> (1.125505273438091e+04 , 4.353951520513552e+03), 
             std::complex<double> (8.731299298615784e+03 , -2.205061951201345e+03), 
             std::complex<double> (4.821218547607376e+03 , -3.981546898712687e+03), 
             std::complex<double> (2.291301388881804e+03 , -3.506897132566718e+03), 
             std::complex<double> (1.019882951187834e+03 , -2.565937506878220e+03), 
             std::complex<double> (4.413420115415439e+02 , -1.752379631729090e+03), 
             std::complex<double> (1.887019051369650e+02 , -1.163871180943330e+03), 
             std::complex<double> (8.027067541448241e+01 , -7.640927517247154e+02), 
             std::complex<double> (3.407182450864039e+01 , -4.991938238083632e+02), 
             std::complex<double> (1.444888609352675e+01 , -3.254593517450245e+02), 
             std::complex<double> (6.124970336280453e+00 , -2.120036162859004e+02), 
             std::complex<double> (2.595982967402665e+00 , -1.380459364857038e+02), 
             std::complex<double> (1.100194160756966e+00 , -8.987156986420281e+01), 
             std::complex<double> (4.662556779406380e-01 , -5.850053510951644e+01), 
             std::complex<double> (1.975941883547968e-01 , -3.807225347093418e+01), 
             std::complex<double> (8.373855338311331e-02 , -2.476685751773691e+01), 
             std::complex<double> (3.548908858908655e-02 , -1.609534320390479e+01), 
             std::complex<double> (1.504403602114764e-02 , -1.043510050380000e+01), 
             std::complex<double> (6.383863441956811e-03 , -6.726701937056735e+00), 
             std::complex<double> (2.706983850738468e-03 , -4.281317986782434e+00), 
             std::complex<double> (1.021562492111798e-03 , -2.726418339001828e+00), 
             std::complex<double> (1.601653948149719e-04 , -2.076920072659586e+00), 
             std::complex<double> (2.931097181302017e-06 , -2.001061326872212e+00), 
             std::complex<double> (2.157737582301619e-09 , -2.000000701952047e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 3.007317109997245e+04;
             std::complex<double> zvec1[] = {
             std::complex<double> (-5.763716124685593e+04 , 1.165465401139088e+04), 
             std::complex<double> (-4.255967779642986e+04 , 2.723952035457892e+04), 
             std::complex<double> (-2.498789059908102e+04 , 2.957574180477994e+04), 
             std::complex<double> (-1.277306086810044e+04 , 2.455919200828622e+04), 
             std::complex<double> (-6.073676156877335e+03 , 1.809679059392678e+04), 
             std::complex<double> (-2.789292151987068e+03 , 1.263172610312310e+04), 
             std::complex<double> (-1.260530761537012e+03 , 8.604388028934538e+03), 
             std::complex<double> (-5.655193133145067e+02 , 5.797260802748918e+03), 
             std::complex<double> (-2.528817125539699e+02 , 3.886872101133649e+03), 
             std::complex<double> (-1.129134377699256e+02 , 2.600335518344097e+03), 
             std::complex<double> (-5.038216018327863e+01 , 1.737949060749353e+03), 
             std::complex<double> (-2.247252983584780e+01 , 1.161075313340091e+03), 
             std::complex<double> (-1.002085979433980e+01 , 7.755501031922179e+02), 
             std::complex<double> (-4.466686434856854e+00 , 5.180182421633351e+02), 
             std::complex<double> (-1.989406893927575e+00 , 3.460315999056143e+02), 
             std::complex<double> (-8.845260760434258e-01 , 2.312044701234011e+02), 
             std::complex<double> (-3.917461292853863e-01 , 1.545739527782731e+02), 
             std::complex<double> (-1.719568022575579e-01 , 1.034820169382285e+02), 
             std::complex<double> (-7.390638162734330e-02 , 6.948909605202205e+01), 
             std::complex<double> (-3.012597162166876e-02 , 4.698179659354058e+01), 
             std::complex<double> (-1.056161738205788e-02 , 3.223629250475251e+01), 
             std::complex<double> (-2.299040172398579e-03 , 2.264983681159068e+01), 
             std::complex<double> (-1.455146710255502e-04 , 1.573834944559516e+01), 
             std::complex<double> (-1.098170565241017e-06 , 9.424967099354957e+00), 
             std::complex<double> (-3.570426412894763e-10 , 3.141592710620479e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (2.921769481203409e+03 , 6.929288266620085e+03), 
             std::complex<double> (5.875580712897955e+03 , 2.709797573361160e+03), 
             std::complex<double> (4.895531659638326e+03 , -8.291183275307773e+02), 
             std::complex<double> (2.909442602462147e+03 , -2.040449365425031e+03), 
             std::complex<double> (1.479181319064367e+03 , -1.955418251667207e+03), 
             std::complex<double> (6.998810798998495e+02 , -1.507481565713051e+03), 
             std::complex<double> (3.205291865534780e+02 , -1.070482525760811e+03), 
             std::complex<double> (1.446583778196000e+02 , -7.343979805786191e+02), 
             std::complex<double> (6.485881886834673e+01 , -4.963276725193895e+02), 
             std::complex<double> (2.899485769148739e+01 , -3.332193030241388e+02), 
             std::complex<double> (1.294507424770513e+01 , -2.230559374174589e+02), 
             std::complex<double> (5.776098444013410e+00 , -1.491162266744689e+02), 
             std::complex<double> (2.576627275451514e+00 , -9.962597495123035e+01), 
             std::complex<double> (1.149259954268274e+00 , -6.654020323158292e+01), 
             std::complex<double> (5.125815541291086e-01 , -4.443170541123320e+01), 
             std::complex<double> (2.286124573766687e-01 , -2.965935315958090e+01), 
             std::complex<double> (1.019627235328294e-01 , -1.978585973333206e+01), 
             std::complex<double> (4.548023999730320e-02 , -1.318079363563392e+01), 
             std::complex<double> (2.029550000422393e-02 , -8.752759537532944e+00), 
             std::complex<double> (9.071864751949677e-03 , -5.769832343111419e+00), 
             std::complex<double> (4.030059683252116e-03 , -3.749466414165739e+00), 
             std::complex<double> (1.451229555743454e-03 , -2.484390839322632e+00), 
             std::complex<double> (1.712183038289819e-04 , -2.038362296781771e+00), 
             std::complex<double> (2.167261351065470e-06 , -2.000380680168940e+00), 
             std::complex<double> (1.186949757225600e-09 , -2.000000190218532e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 1.797557536969104e+04;
             std::complex<double> zvec1[] = {
             std::complex<double> (-3.460969817815411e+04 , 6.604163400620757e+03), 
             std::complex<double> (-2.636328557160672e+04 , 1.583115477760520e+04), 
             std::complex<double> (-1.625596855227097e+04 , 1.785618883079884e+04), 
             std::complex<double> (-8.770309253532130e+03 , 1.541661301167759e+04), 
             std::complex<double> (-4.396596187464927e+03 , 1.176331902232053e+04), 
             std::complex<double> (-2.123446863653772e+03 , 8.465172035237203e+03), 
             std::complex<double> (-1.007174471780252e+03 , 5.925641920122363e+03), 
             std::complex<double> (-4.736156401012078e+02 , 4.094465565486505e+03), 
             std::complex<double> (-2.218099685058050e+02 , 2.812033428867710e+03), 
             std::complex<double> (-1.036810929452091e+02 , 1.925803269574458e+03), 
             std::complex<double> (-4.841816111632444e+01 , 1.317137284141384e+03), 
             std::complex<double> (-2.259881773345828e+01 , 9.003042660478659e+02), 
             std::complex<double> (-1.054313304823508e+01 , 6.152326156680514e+02), 
             std::complex<double> (-4.915640162738932e+00 , 4.204050159381569e+02), 
             std::complex<double> (-2.289122786420149e+00 , 2.873085503224958e+02), 
             std::complex<double> (-1.063322073021556e+00 , 1.964204880644473e+02), 
             std::complex<double> (-4.912446948332309e-01 , 1.343953711152879e+02), 
             std::complex<double> (-2.242385078438374e-01 , 9.212186543170240e+01), 
             std::complex<double> (-9.957086303861180e-02 , 6.339017821090112e+01), 
             std::complex<double> (-4.127886824448182e-02 , 4.398191818846385e+01), 
             std::complex<double> (-1.406327199444138e-02 , 3.103112718030820e+01), 
             std::complex<double> (-2.641842790454905e-03 , 2.237366989291479e+01), 
             std::complex<double> (-1.248859612469691e-04 , 1.572135530159720e+01), 
             std::complex<double> (-6.977813433560928e-07 , 9.424841295061135e+00), 
             std::complex<double> (-1.762399828940226e-10 , 3.141592668563918e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.566848292974774e+03 , 3.956111564469049e+03), 
             std::complex<double> (3.281529879093910e+03 , 1.747052473007945e+03), 
             std::complex<double> (2.910141801015012e+03 , -2.736333912391690e+02), 
             std::complex<double> (1.847249415583363e+03 , -1.098124351410455e+03), 
             std::complex<double> (9.985198464134793e+02 , -1.149191098232400e+03), 
             std::complex<double> (4.995162872032331e+02 , -9.330104659811948e+02), 
             std::complex<double> (2.408478215764314e+02 , -6.880287253681282e+02), 
             std::complex<double> (1.141289969083505e+02 , -4.867158044793719e+02), 
             std::complex<double> (5.364310472868186e+01 , -3.379000044651108e+02), 
             std::complex<double> (2.511760104178017e+01 , -2.325668424800111e+02), 
             std::complex<double> (1.174004682124476e+01 , -1.594279395120610e+02), 
             std::complex<double> (5.482779020623415e+00 , -1.090846942380697e+02), 
             std::complex<double> (2.559547471449550e+00 , -7.457079969126342e+01), 
             std::complex<double> (1.194666971096660e+00 , -5.095194700777535e+01), 
             std::complex<double> (5.575631164482663e-01 , -3.480094155426643e+01), 
             std::complex<double> (2.602120672346646e-01 , -2.375808120660787e+01), 
             std::complex<double> (1.214423721937116e-01 , -1.620477599368805e+01), 
             std::complex<double> (5.668781407579560e-02 , -1.103208654541415e+01), 
             std::complex<double> (2.648230944082334e-02 , -7.479622938044910e+00), 
             std::complex<double> (1.239954779449508e-02 , -5.025240853463135e+00), 
             std::complex<double> (5.693487289688867e-03 , -3.327879562175432e+00), 
             std::complex<double> (1.860537802467127e-03 , -2.310260680395144e+00), 
             std::complex<double> (1.598197294753274e-04 , -2.018152106357870e+00), 
             std::complex<double> (1.437389085953849e-06 , -2.000132617367441e+00), 
             std::complex<double> (5.993178254446946e-10 , -2.000000051066325e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 1.138274331919166e+04;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.200140744024029e+04 , 3.974462749966869e+03), 
             std::complex<double> (-1.721129861413965e+04 , 9.736596746171608e+03), 
             std::complex<double> (-1.107683046705493e+04 , 1.135617784198418e+04), 
             std::complex<double> (-6.270456593459743e+03 , 1.015591397626164e+04), 
             std::complex<double> (-3.296802846550238e+03 , 8.002049672433395e+03), 
             std::complex<double> (-1.666587789549242e+03 , 5.923407880623897e+03), 
             std::complex<double> (-8.258198601294596e+02 , 4.252129119762406e+03), 
             std::complex<double> (-4.051584607219896e+02 , 3.006872174944886e+03), 
             std::complex<double> (-1.978041658306479e+02 , 2.110771324703107e+03), 
             std::complex<double> (-9.633645806488423e+01 , 1.476455451524824e+03), 
             std::complex<double> (-4.685987604832491e+01 , 1.030983850084752e+03), 
             std::complex<double> (-2.277632239555227e+01 , 7.193329136413156e+02), 
             std::complex<double> (-1.106309994783381e+01 , 5.017148012703821e+02), 
             std::complex<double> (-5.368608420831922e+00 , 3.499073737872292e+02), 
             std::complex<double> (-2.600722961499763e+00 , 2.440742966261574e+02), 
             std::complex<double> (-1.255471158259803e+00 , 1.703362234370004e+02), 
             std::complex<double> (-6.016531106911772e-01 , 1.190060952464079e+02), 
             std::complex<double> (-2.838368939061326e-01 , 8.333572503097746e+01), 
             std::complex<double> (-1.292429994736169e-01 , 5.863595433538533e+01), 
             std::complex<double> (-5.388850307076155e-02 , 4.166197489437655e+01), 
             std::complex<double> (-1.744760289143409e-02 , 3.014999937825728e+01), 
             std::complex<double> (-2.737183175197586e-03 , 2.220418688588807e+01), 
             std::complex<double> (-9.607833698541100e-05 , 1.571363315166138e+01), 
             std::complex<double> (-4.065179810660882e-07 , 9.424798719304857e+00), 
             std::complex<double> (-8.111167317574562e-11 , 3.141592657495105e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (8.948400255004358e+02 , 2.395959021318808e+03), 
             std::complex<double> (1.940568522442985e+03 , 1.166267085546447e+03), 
             std::complex<double> (1.817740498629313e+03 , -4.527450160958868e+01), 
             std::complex<double> (1.224147158574178e+03 , -6.134315606782586e+02), 
             std::complex<double> (6.997510751916801e+02 , -7.050700883448437e+02), 
             std::complex<double> (3.683907051618648e+02 , -6.028359017261160e+02), 
             std::complex<double> (1.861824845141888e+02 , -4.612320277751987e+02), 
             std::complex<double> (9.222591723168787e+01 , -3.359934975261960e+02), 
             std::complex<double> (4.523836496502550e+01 , -2.392227827861900e+02), 
             std::complex<double> (2.208433347127993e+01 , -1.684734636890618e+02), 
             std::complex<double> (1.075600442417885e+01 , -1.180242233851356e+02), 
             std::complex<double> (5.232704271656762e+00 , -8.246940754747892e+01), 
             std::complex<double> (2.544266960624304e+00 , -5.755083865962865e+01), 
             std::complex<double> (1.236754346410644e+00 , -4.013233041302669e+01), 
             std::complex<double> (6.011041566788931e-01 , -2.797018277385432e+01), 
             std::complex<double> (2.921443382225316e-01 , -1.948030300227659e+01), 
             std::complex<double> (1.419937333043957e-01 , -1.355088645076997e+01), 
             std::complex<double> (6.903797504032956e-02 , -9.403150330218590e+00), 
             std::complex<double> (3.361141043927978e-02 , -6.491072225419416e+00), 
             std::complex<double> (1.640320159366630e-02 , -4.432289692738207e+00), 
             std::complex<double> (7.640741457627663e-03 , -2.990804154053763e+00), 
             std::complex<double> (2.151485487405385e-03 , -2.189855236246026e+00), 
             std::complex<double> (1.323173537777747e-04 , -2.008182286479190e+00), 
             std::complex<double> (8.699166744696476e-07 , -2.000045041084302e+00), 
             std::complex<double> (2.816290368664371e-10 , -2.000000013594450e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 7.568486488496546e+03;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.467760939146840e+04 , 2.517250040114760e+03), 
             std::complex<double> (-1.174862011184371e+04 , 6.283604305282818e+03), 
             std::complex<double> (-7.850635752005036e+03 , 7.548873761978218e+03), 
             std::complex<double> (-4.639386664936884e+03 , 6.969532465280259e+03), 
             std::complex<double> (-2.546755726761585e+03 , 5.656337791774859e+03), 
             std::complex<double> (-1.341983527394147e+03 , 4.298365091706187e+03), 
             std::complex<double> (-6.919668892184336e+02 , 3.158600221046387e+03), 
             std::complex<double> (-3.528119077706226e+02 , 2.281824362836700e+03), 
             std::complex<double> (-1.788528477590552e+02 , 1.634270162413233e+03), 
             std::complex<double> (-9.039651409964291e+01 , 1.165405670515718e+03), 
             std::complex<double> (-4.561438300652406e+01 , 8.292495986235490e+02), 
             std::complex<double> (-2.299318516720757e+01 , 5.894296807432193e+02), 
             std::complex<double> (-1.157922597523006e+01 , 4.187706374830202e+02), 
             std::complex<double> (-5.823394665244928e+00 , 2.974939287813858e+02), 
             std::complex<double> (-2.921672132969186e+00 , 2.113865223849131e+02), 
             std::complex<double> (-1.458989022742513e+00 , 1.503007904478883e+02), 
             std::complex<double> (-7.216764574349218e-01 , 1.070181904946213e+02), 
             std::complex<double> (-3.498997534832062e-01 , 7.641812482962595e+01), 
             std::complex<double> (-1.622208465026699e-01 , 5.488096489165680e+01), 
             std::complex<double> (-6.723123054533840e-02 , 3.986020028729303e+01), 
             std::complex<double> (-2.018898535585851e-02 , 2.951631519746510e+01), 
             std::complex<double> (-2.569447008291321e-03 , 2.210474313643918e+01), 
             std::complex<double> (-6.722069177672878e-05 , 1.571027852231603e+01), 
             std::complex<double> (-2.199444639199661e-07 , 9.424784636434616e+00), 
             std::complex<double> (-3.514385821565780e-11 , 3.141592654602162e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (5.391886832116896e+02 , 1.525715272787130e+03), 
             std::complex<double> (1.205027547524299e+03 , 8.042918020751428e+02), 
             std::complex<double> (1.184513458919501e+03 , 4.644485737686573e+01), 
             std::complex<double> (8.413512386330929e+02 , -3.519254019726017e+02), 
             std::complex<double> (5.061656624147366e+02 , -4.481384869169310e+02), 
             std::complex<double> (2.792982255840163e+02 , -4.036851425953848e+02), 
             std::complex<double> (1.474007367453372e+02 , -3.202565633707214e+02), 
             std::complex<double> (7.604284769304232e+01 , -2.400028134129853e+02), 
             std::complex<double> (3.878004209292357e+01 , -1.750194130940599e+02), 
             std::complex<double> (1.966186363697081e+01 , -1.259279757456473e+02), 
             std::complex<double> (9.939396063222505e+00 , -9.000042106543675e+01), 
             std::complex<double> (5.017018692892879e+00 , -6.410498463890151e+01), 
             std::complex<double> (2.530468977524313e+00 , -4.557930222776973e+01), 
             std::complex<double> (1.275820342683220e+00 , -3.237391500566422e+01), 
             std::complex<double> (6.431294800847643e-01 , -2.297613490246886e+01), 
             std::complex<double> (3.241789590653019e-01 , -1.629085192805470e+01), 
             std::complex<double> (1.634276754601557e-01 , -1.153224301855493e+01), 
             std::complex<double> (8.243930123182515e-02 , -8.138174614809062e+00), 
             std::complex<double> (4.167260527403722e-02 , -5.706198008980571e+00), 
             std::complex<double> (2.108790107167119e-02 , -3.950939354890322e+00), 
             std::complex<double> (9.712589211716362e-03 , -2.721868835745906e+00), 
             std::complex<double> (2.244804969055464e-03 , -2.110597529401569e+00), 
             std::complex<double> (9.866673914841044e-05 , -2.003530954312899e+00), 
             std::complex<double> (4.870594346865956e-07 , -2.000014959693225e+00), 
             std::complex<double> (1.244127935221651e-10 , -2.000000003591815e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 5.246506559650061e+03;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.020379184122255e+04 , 1.665718296083240e+03), 
             std::complex<double> (-8.331893883993624e+03 , 4.226375054838867e+03), 
             std::complex<double> (-5.755037726090139e+03 , 5.212300281293105e+03), 
             std::complex<double> (-3.534709563532957e+03 , 4.953282830018878e+03), 
             std::complex<double> (-2.017871240732303e+03 , 4.131282031101155e+03), 
             std::complex<double> (-1.104379196371829e+03 , 3.217167909117524e+03), 
             std::complex<double> (-5.905474438453264e+02 , 2.416256307890350e+03), 
             std::complex<double> (-3.118677656744005e+02 , 1.780554340499586e+03), 
             std::complex<double> (-1.636070849938727e+02 , 1.299115725599460e+03), 
             std::complex<double> (-8.552375116492996e+01 , 9.429526861414693e+02), 
             std::complex<double> (-4.461663197920470e+01 , 6.826044209372221e+02), 
             std::complex<double> (-2.324449480231984e+01 , 4.934734540951367e+02), 
             std::complex<double> (-1.209431005358674e+01 , 3.565298069888012e+02), 
             std::complex<double> (-6.281263753520784e+00 , 2.575560396188275e+02), 
             std::complex<double> (-3.251742396166525e+00 , 1.861116000891881e+02), 
             std::complex<double> (-1.673091560113771e+00 , 1.345980658884251e+02), 
             std::complex<double> (-8.504302711413242e-01 , 9.751417158504238e+01), 
             std::complex<double> (-4.215340472539931e-01 , 7.089283008661154e+01), 
             std::complex<double> (-1.975530955691623e-01 , 5.188698902841897e+01), 
             std::complex<double> (-8.031551977901168e-02 , 3.846192462598919e+01), 
             std::complex<double> (-2.178909320570640e-02 , 2.907131174387766e+01), 
             std::complex<double> (-2.197320846097559e-03 , 2.204913563487819e+01), 
             std::complex<double> (-4.331291559333815e-05 , 1.570887873472113e+01), 
             std::complex<double> (-1.116591399229271e-07 , 9.424780071515466e+00), 
             std::complex<double> (-1.448140937043927e-11 , 3.141592653850790e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (3.402134027865599e+02 , 1.014273262603915e+03), 
             std::complex<double> (7.804893207887276e+02 , 5.713882236897982e+02), 
             std::complex<double> (8.006371577385115e+02 , 7.944837260986665e+01), 
             std::complex<double> (5.967026115395553e+02 , -2.051675700245182e+02), 
             std::complex<double> (3.762083691395865e+02 , -2.932212403031658e+02), 
             std::complex<double> (2.167984604006476e+02 , -2.785474571115707e+02), 
             std::complex<double> (1.190876253106622e+02 , -2.290496942028836e+02), 
             std::complex<double> (6.377944688635845e+01 , -1.764520603595772e+02), 
             std::complex<double> (3.370726766194443e+01 , -1.316609551501948e+02), 
             std::complex<double> (1.769091633035518e+01 , -9.666547164092943e+01), 
             std::complex<double> (9.251386977020044e+00 , -7.038394923585382e+01), 
             std::complex<double> (4.828921870815317e+00 , -5.102537341016783e+01), 
             std::complex<double> (2.518123109491837e+00 , -3.690419038662468e+01), 
             std::complex<double> (1.312490172778282e+00 , -2.665334993586288e+01), 
             std::complex<double> (6.839458933948733e-01 , -1.922876299119833e+01), 
             std::complex<double> (3.563941011719078e-01 , -1.385466320461801e+01), 
             std::complex<double> (1.857490649552708e-01 , -9.961968656937502e+00), 
             std::complex<double> (9.689947577428010e-02 , -7.135167179860161e+00), 
             std::complex<double> (5.069362060898600e-02 , -5.070869353041386e+00), 
             std::complex<double> (2.642099148952734e-02 , -3.554319620235563e+00), 
             std::complex<double> (1.163824127286680e-02 , -2.510264722885624e+00), 
             std::complex<double> (2.118094708600942e-03 , -2.061247874648716e+00), 
             std::complex<double> (6.720150544399583e-05 , -2.001466144053758e+00), 
             std::complex<double> (2.550583408283563e-07 , -2.000004871632734e+00), 
             std::complex<double> (5.218622440454387e-11 , -2.000000000942600e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 25;
            int m = 2;
            double b = 3.769320093606350e+03;
             std::complex<double> zvec1[] = {
             std::complex<double> (-7.349137635536470e+03 , 1.144610275291591e+03), 
             std::complex<double> (-6.106209195403381e+03 , 2.945908606179457e+03), 
             std::complex<double> (-4.343272122031992e+03 , 3.718881760554523e+03), 
             std::complex<double> (-2.761819916629966e+03 , 3.627997699062357e+03), 
             std::complex<double> (-1.633929200490645e+03 , 3.103261322358168e+03), 
             std::complex<double> (-9.258825245736423e+02 , 2.472428695538470e+03), 
             std::complex<double> (-5.119239934276529e+02 , 1.895259411566562e+03), 
             std::complex<double> (-2.792039689110888e+02 , 1.422779869783170e+03), 
             std::complex<double> (-1.511368189771692e+02 , 1.056118878267329e+03), 
             std::complex<double> (-8.147004939085437e+01 , 7.792226459025721e+02), 
             std::complex<double> (-4.380725551985272e+01 , 5.730792947877991e+02), 
             std::complex<double> (-2.351419103912193e+01 , 4.207738280011841e+02), 
             std::complex<double> (-1.259973001923183e+01 , 3.087117708799910e+02), 
             std::complex<double> (-6.735100629379208e+00 , 2.264570621414299e+02), 
             std::complex<double> (-3.585382823343307e+00 , 1.661789241917392e+02), 
             std::complex<double> (-1.894006887673428e+00 , 1.220721982555267e+02), 
             std::complex<double> (-9.855586380914142e-01 , 8.986395348198425e+01), 
             std::complex<double> (-4.971788166440539e-01 , 6.642625967767339e+01), 
             std::complex<double> (-2.339369550657568e-01 , 4.948374102186163e+01), 
             std::complex<double> (-9.198207201044944e-02 , 3.738249144648803e+01), 
             std::complex<double> (-2.195370938567570e-02 , 2.876801835231873e+01), 
             std::complex<double> (-1.726212578816507e-03 , 2.201951166967633e+01), 
             std::complex<double> (-2.604155477773438e-05 , 1.570831479203366e+01), 
             std::complex<double> (-5.376690136381725e-08 , 9.424778617555379e+00), 
             std::complex<double> (-5.719436824667930e-12 , 3.141592653656703e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (2.233825342648173e+02 , 6.997558576362865e+02), 
             std::complex<double> (5.243206284476087e+02 , 4.169410534379146e+02), 
             std::complex<double> (5.586255964813120e+02 , 8.706860210024267e+01), 
             std::complex<double> (4.348509883081007e+02 , -1.200817181055730e+02), 
             std::complex<double> (2.862248321565207e+02 , -1.964373688943619e+02), 
             std::complex<double> (1.717110689432833e+02 , -1.970910866697924e+02), 
             std::complex<double> (9.789130005802932e+01 , -1.679605081114850e+02), 
             std::complex<double> (5.427796693472952e+01 , -1.329339945876383e+02), 
             std::complex<double> (2.964650279123934e+01 , -1.014094371657507e+02), 
             std::complex<double> (1.606224507349390e+01 , -7.589980212477514e+01), 
             std::complex<double> (8.664614769687560e+00 , -5.623732681565048e+01), 
             std::complex<double> (4.663122228494752e+00 , -4.144303669529222e+01), 
             std::complex<double> (2.506450404835361e+00 , -3.044836675279958e+01), 
             std::complex<double> (1.346345905218170e+00 , -2.232889958620853e+01), 
             std::complex<double> (7.229884064051585e-01 , -1.635070224156870e+01), 
             std::complex<double> (3.882566146409886e-01 , -1.195315379384877e+01), 
             std::complex<double> (2.086028783026887e-01 , -8.715660270887964e+00), 
             std::complex<double> (1.122623477384150e-01 , -6.324797219658176e+00), 
             std::complex<double> (6.064720648773948e-02 , -4.547561780135394e+00), 
             std::complex<double> (3.227589977984663e-02 , -3.224098175845442e+00), 
             std::complex<double> (1.308601967709231e-02 , -2.347867283048108e+00), 
             std::complex<double> (1.818763283925938e-03 , -2.032252584839498e+00), 
             std::complex<double> (4.241752531159779e-05 , -2.000588093530181e+00), 
             std::complex<double> (1.263489618030886e-07 , -2.000001557577005e+00), 
             std::complex<double> (2.096735891755558e-11 , -2.000000000245687e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 30;
            int m = 2;
            double b = 1.362802852101742e+11;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.260563430113050e+11 , 9.899300161597076e+10), 
             std::complex<double> (-7.137099207454512e+10 , 1.188865095897478e+11), 
             std::complex<double> (-1.424304992986035e+10 , 6.028522251575875e+10), 
             std::complex<double> (-2.563890957766498e+09 , 2.615637141197538e+10), 
             std::complex<double> (-4.527976369835516e+08 , 1.103547342985601e+10), 
             std::complex<double> (-7.969601739563417e+07 , 4.632956560581381e+09), 
             std::complex<double> (-1.401876438849539e+07 , 1.943335272312015e+09), 
             std::complex<double> (-2.465682916825068e+06 , 8.150247372464416e+08), 
             std::complex<double> (-4.336673142192835e+05 , 3.418079505644055e+08), 
             std::complex<double> (-7.627368744880818e+04 , 1.433479460818319e+08), 
             std::complex<double> (-1.341505775291618e+04 , 6.011743292117146e+07), 
             std::complex<double> (-2.359447432738017e+03 , 2.521211716431283e+07), 
             std::complex<double> (-4.149808500901972e+02 , 1.057348600316693e+07), 
             std::complex<double> (-7.298704903278420e+01 , 4.434320411262261e+06), 
             std::complex<double> (-1.283700045347998e+01 , 1.859670262788693e+06), 
             std::complex<double> (-2.257778267474294e+00 , 7.799105984904652e+05), 
             std::complex<double> (-3.970991745646501e-01 , 3.270797806601724e+05), 
             std::complex<double> (-6.984200945055032e-02 , 1.371710849341482e+05), 
             std::complex<double> (-1.228385149331165e-02 , 5.752696340436584e+04), 
             std::complex<double> (-2.160488920004809e-03 , 2.412572263657624e+04), 
             std::complex<double> (-3.799858292012274e-04 , 1.011787323100246e+04), 
             std::complex<double> (-6.683059995445239e-05 , 4.243246538709969e+03), 
             std::complex<double> (-1.175324304299886e-05 , 1.779541127848084e+03), 
             std::complex<double> (-2.066709078414816e-06 , 7.463146315152027e+02), 
             std::complex<double> (-3.632561502681531e-07 , 3.130109576500223e+02), 
             std::complex<double> (-6.369356947468524e-08 , 1.313201614181625e+02), 
             std::complex<double> (-1.101427091490508e-08 , 5.518980769571331e+01), 
             std::complex<double> (-1.751054480068750e-09 , 2.340886474912669e+01), 
             std::complex<double> (-1.674126285839405e-10 , 1.015828114798508e+01), 
             std::complex<double> (-1.635747317515906e-12 , 3.147253059075087e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (4.919217451488852e+10 , 4.539566695439857e+10), 
             std::complex<double> (3.364564474083897e+10 , -1.792345799212862e+10), 
             std::complex<double> (7.658960683588418e+09 , -1.530392548277468e+10), 
             std::complex<double> (1.411289253112843e+09 , -7.129695784316122e+09), 
             std::complex<double> (2.502705086579407e+08 , -3.044631421551719e+09), 
             std::complex<double> (4.408155719307550e+07 , -1.280914113259307e+09), 
             std::complex<double> (7.755064486344404e+06 , -5.374905263660849e+08), 
             std::complex<double> (1.364025977297759e+06 , -2.254354628001529e+08), 
             std::complex<double> (2.399074961948096e+05 , -9.454500594733365e+07), 
             std::complex<double> (4.219511473185286e+04 , -3.965051068154848e+07), 
             std::complex<double> (7.421300706560492e+03 , -1.662868445845063e+07), 
             std::complex<double> (1.305262312540593e+03 , -6.973756973368203e+06), 
             std::complex<double> (2.295702203474029e+02 , -2.924662061392994e+06), 
             std::complex<double> (4.037693036399335e+01 , -1.226548055774383e+06), 
             std::complex<double> (7.101515694977040e+00 , -5.143910981667502e+05), 
             std::complex<double> (1.249018314044843e+00 , -2.157259151110064e+05), 
             std::complex<double> (2.196779705719005e-01 , -9.047137598333800e+04), 
             std::complex<double> (3.863707057795756e-02 , -3.794198701476835e+04), 
             std::complex<double> (6.795509471537484e-03 , -1.591215299486408e+04), 
             std::complex<double> (1.195198476973741e-03 , -6.673256477342017e+03), 
             std::complex<double> (2.102120245046098e-04 , -2.798637601364831e+03), 
             std::complex<double> (3.697201442680757e-05 , -1.173695350701868e+03), 
             std::complex<double> (6.502482611233163e-06 , -4.922246643756272e+02), 
             std::complex<double> (1.143537794983518e-06 , -2.064273307154715e+02), 
             std::complex<double> (2.010977199820400e-07 , -8.656600503307962e+01), 
             std::complex<double> (3.536354387973731e-08 , -3.629056062542763e+01), 
             std::complex<double> (6.219484889780642e-09 , -1.518793607688082e+01), 
             std::complex<double> (1.090944126469887e-09 , -6.308038107366983e+00), 
             std::complex<double> (1.599569759616599e-10 , -2.743980107228885e+00), 
             std::complex<double> (3.690635197941082e-12 , -2.013099043213039e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 30;
            int m = 2;
            double b = 7.777560169487656e+09;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.340036341986429e+10 , 5.193379656137459e+09), 
             std::complex<double> (-5.186240584073159e+09 , 7.282779161291145e+09), 
             std::complex<double> (-1.298978142551649e+09 , 4.281812174394330e+09), 
             std::complex<double> (-2.882972405851151e+08 , 2.088165317958621e+09), 
             std::complex<double> (-6.223865292177729e+07 , 9.774538142130129e+08), 
             std::complex<double> (-1.335569682272292e+07 , 4.535130729843550e+08), 
             std::complex<double> (-2.862274309657053e+06 , 2.100197200753430e+08), 
             std::complex<double> (-6.132471539992122e+05 , 9.721961763471982e+07), 
             std::complex<double> (-1.313814437510622e+05 , 4.499973350262288e+07), 
             std::complex<double> (-2.814666805905915e+04 , 2.082849434449822e+07), 
             std::complex<double> (-6.030020222776345e+03 , 9.640600284193022e+06), 
             std::complex<double> (-1.291844652912483e+03 , 4.462208970071673e+06), 
             std::complex<double> (-2.767590046687632e+02 , 2.065359485395700e+06), 
             std::complex<double> (-5.929160573872303e+01 , 9.559636639607314e+05), 
             std::complex<double> (-1.270236692955135e+01 , 4.424733445272422e+05), 
             std::complex<double> (-2.721297873507507e+00 , 2.048013618887120e+05), 
             std::complex<double> (-5.829985617043429e-01 , 9.479350194798034e+04), 
             std::complex<double> (-1.248989590203771e-01 , 4.387572406281675e+04), 
             std::complex<double> (-2.675777797803507e-02 , 2.030813456137690e+04), 
             std::complex<double> (-5.732458834817584e-03 , 9.399739133018220e+03), 
             std::complex<double> (-1.228090375508765e-03 , 4.350725892767989e+03), 
             std::complex<double> (-2.630946076463827e-04 , 2.013763002066167e+03), 
             std::complex<double> (-5.635895693398096e-05 , 9.320910185808137e+02), 
             std::complex<double> (-1.206968716420404e-05 , 4.314434038509972e+02), 
             std::complex<double> (-2.581638244107004e-06 , 1.997385331681191e+02), 
             std::complex<double> (-5.490055153730659e-07 , 9.254166528862122e+01), 
             std::complex<double> (-1.135505618228694e-07 , 4.302918264817856e+01), 
             std::complex<double> (-2.037507137607511e-08 , 2.030895866731232e+01), 
             std::complex<double> (-1.776287578310367e-09 , 9.774042277542037e+00), 
             std::complex<double> (-1.099461804086121e-11 , 3.143296733062273e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (2.337134082612616e+09 , 2.562343703782293e+09), 
             std::complex<double> (2.059114538739392e+09 , -7.125836933970790e+08), 
             std::complex<double> (6.086712902069863e+08 , -9.108527219967749e+08), 
             std::complex<double> (1.400089150081459e+08 , -4.973841952581840e+08), 
             std::complex<double> (3.045875982471164e+07 , -2.382066569219896e+08), 
             std::complex<double> (6.546872627964534e+06 , -1.110581261191356e+08), 
             std::complex<double> (1.403562825822862e+06 , -5.148373087887733e+07), 
             std::complex<double> (3.007384670227370e+05 , -2.383746524327619e+07), 
             std::complex<double> (6.443094553630956e+04 , -1.103409500840049e+07), 
             std::complex<double> (1.380349216431272e+04 , -5.107272965336684e+06), 
             std::complex<double> (2.957202558297040e+03 , -2.363938729470136e+06), 
             std::complex<double> (6.335379950957779e+02 , -1.094163559852584e+06), 
             std::complex<double> (1.357263438120013e+02 , -5.064400452784285e+05), 
             std::complex<double> (2.907740259125273e+01 , -2.344087288652768e+05), 
             std::complex<double> (6.229412046489145e+00 , -1.084974446461929e+05), 
             std::complex<double> (1.334561210617000e+00 , -5.021867353745985e+04), 
             std::complex<double> (2.859103826622169e-01 , -2.324400519564840e+04), 
             std::complex<double> (6.125214739452008e-02 , -1.075862298768337e+04), 
             std::complex<double> (1.312238096865155e-02 , -4.979691139313308e+03), 
             std::complex<double> (2.811279011398451e-03 , -2.304878793602031e+03), 
             std::complex<double> (6.022757187639150e-04 , -1.066826081087913e+03), 
             std::complex<double> (1.290287982182687e-04 , -4.937856257218228e+02), 
             std::complex<double> (2.764221827286943e-05 , -2.285493189940098e+02), 
             std::complex<double> (5.921774023805938e-06 , -1.057805639172012e+02), 
             std::complex<double> (1.268618871986980e-06 , -4.895075262853618e+01), 
             std::complex<double> (2.717802630711590e-07 , -2.263481099889794e+01), 
             std::complex<double> (5.822038762103108e-08 , -1.043002305857571e+01), 
             std::complex<double> (1.235258522685602e-08 , -4.754756830168871e+00), 
             std::complex<double> (1.851266663506792e-09 , -2.392986873159544e+00), 
             std::complex<double> (2.675316853230272e-11 , -2.004221935414854e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 30;
            int m = 2;
            double b = 8.050742173805323e+08;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.425932735889674e+09 , 4.957024150292232e+08), 
             std::complex<double> (-6.447924439297898e+08 , 7.840828213828845e+08), 
             std::complex<double> (-1.940555352739612e+08 , 5.220084337206816e+08), 
             std::complex<double> (-5.107963994105916e+07 , 2.811213145805761e+08), 
             std::complex<double> (-1.296455480427407e+07 , 1.433617566893554e+08), 
             std::complex<double> (-3.260006136783800e+06 , 7.210891474402954e+07), 
             std::complex<double> (-8.178219582224345e+05 , 3.614442878099144e+07), 
             std::complex<double> (-2.050420669439653e+05 , 1.810159107396463e+07), 
             std::complex<double> (-5.139997951468418e+04 , 9.063534939492136e+06), 
             std::complex<double> (-1.288447772381847e+04 , 4.537898553064034e+06), 
             std::complex<double> (-3.229733260016483e+03 , 2.271987893797936e+06), 
             std::complex<double> (-8.095906457000415e+02 , 1.137511390681061e+06), 
             std::complex<double> (-2.029383165299460e+02 , 5.695149413375852e+05), 
             std::complex<double> (-5.087009643409584e+01 , 2.851375390434747e+05), 
             std::complex<double> (-1.275149340409899e+01 , 1.427590496966748e+05), 
             std::complex<double> (-3.196388268890898e+00 , 7.147479105827559e+04), 
             std::complex<double> (-8.012314361188517e-01 , 3.578509230676397e+04), 
             std::complex<double> (-2.008428430190844e-01 , 1.791642679721231e+04), 
             std::complex<double> (-5.034478457210373e-02 , 8.970170536023263e+03), 
             std::complex<double> (-1.261977463155531e-02 , 4.491073576192549e+03), 
             std::complex<double> (-3.163329895190986e-03 , 2.248538454549567e+03), 
             std::complex<double> (-7.929019028998871e-04 , 1.125779086684483e+03), 
             std::complex<double> (-1.987098730255687e-04 , 5.636595075324266e+02), 
             std::complex<double> (-4.976385096472798e-05 , 2.822432614208342e+02), 
             std::complex<double> (-1.242744655450400e-05 , 1.413845708914785e+02), 
             std::complex<double> (-3.068333489803973e-06 , 7.093527043253458e+01), 
             std::complex<double> (-7.223158726723104e-07 , 3.580815253213384e+01), 
             std::complex<double> (-1.366075084778147e-07 , 1.845658665214465e+01), 
             std::complex<double> (-9.938834294812918e-09 , 9.584244111714973e+00), 
             std::complex<double> (-3.876361316436709e-11 , 3.142091104567148e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (2.034976356727266e+08 , 2.573183295748415e+08), 
             std::complex<double> (2.180386327005781e+08 , -4.291775631403627e+07), 
             std::complex<double> (7.996721047450256e+07 , -9.269188008511361e+07), 
             std::complex<double> (2.212381711263887e+07 , -5.887024703655051e+07), 
             std::complex<double> (5.685994629279939e+06 , -3.118070142885414e+07), 
             std::complex<double> (1.434270837624724e+06 , -1.583008157969254e+07), 
             std::complex<double> (3.600920060476071e+05 , -7.953232629329551e+06), 
             std::complex<double> (9.029909714178318e+04 , -3.985395849507198e+06), 
             std::complex<double> (2.263731276692567e+04 , -1.995793369572163e+06), 
             std::complex<double> (5.674585329031583e+03 , -9.992832272603641e+05), 
             std::complex<double> (1.422444458828822e+03 , -5.003153203802009e+05), 
             std::complex<double> (3.565615275682567e+02 , -2.504923915818458e+05), 
             std::complex<double> (8.937851583105515e+01 , -1.254134621198605e+05), 
             std::complex<double> (2.240431490228944e+01 , -6.279043533435158e+04), 
             std::complex<double> (5.616039658388154e+00 , -3.143712056635074e+04), 
             std::complex<double> (1.407760110187167e+00 , -1.573953903304089e+04), 
             std::complex<double> (3.528800735085265e-01 , -7.880272762253429e+03), 
             std::complex<double> (8.845565555765417e-02 , -3.945394905775825e+03), 
             std::complex<double> (2.217297994279907e-02 , -1.975330028960530e+03), 
             std::complex<double> (5.558051315024736e-03 , -9.889826714447177e+02), 
             std::complex<double> (1.393224794864747e-03 , -4.951502691628797e+02), 
             std::complex<double> (3.492373147348032e-04 , -2.479034944926362e+02), 
             std::complex<double> (8.754308641287641e-05 , -1.241130577837041e+02), 
             std::complex<double> (2.194447516961261e-05 , -6.213112609457483e+01), 
             std::complex<double> (5.500831781069279e-06 , -3.109062229958284e+01), 
             std::complex<double> (1.378890739493853e-06 , -1.553354716914265e+01), 
             std::complex<double> (3.455819827470462e-07 , -7.715176851779923e+00), 
             std::complex<double> (8.422268247454932e-08 , -3.790763970679506e+00), 
             std::complex<double> (1.146586071259843e-08 , -2.198773295729515e+00), 
             std::complex<double> (1.003402882138138e-10 , -2.001307672244444e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 30;
            int m = 2;
            double b = 1.279905757901758e+08;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.314282913221650e+08 , 7.294579364443612e+07), 
             std::complex<double> (-1.181482424285914e+08 , 1.268871497765626e+08), 
             std::complex<double> (-4.134583896343520e+07 , 9.386031104946351e+07), 
             std::complex<double> (-1.254570225004170e+07 , 5.508621814877839e+07), 
             std::complex<double> (-3.640457200547342e+06 , 3.021513960279294e+07), 
             std::complex<double> (-1.042626668484624e+06 , 1.625357541256095e+07), 
             std::complex<double> (-2.974868185433647e+05 , 8.694728270603180e+06), 
             std::complex<double> (-8.478909459452594e+04 , 4.643803370879969e+06), 
             std::complex<double> (-2.415901581426754e+04 , 2.479106597653165e+06), 
             std::complex<double> (-6.883044409137209e+03 , 1.323307011191927e+06), 
             std::complex<double> (-1.960970704278958e+03 , 7.063339496389678e+05), 
             std::complex<double> (-5.586741342503982e+02 , 3.770118541898234e+05), 
             std::complex<double> (-1.591641124661413e+02 , 2.012327380491327e+05), 
             std::complex<double> (-4.534521828939148e+01 , 1.074093030896666e+05), 
             std::complex<double> (-1.291866927400398e+01 , 5.733041126095514e+04), 
             std::complex<double> (-3.680476388876830e+00 , 3.060047655183592e+04), 
             std::complex<double> (-1.048552553397212e+00 , 1.633320183847794e+04), 
             std::complex<double> (-2.987280900539186e-01 , 8.717952877066295e+03), 
             std::complex<double> (-8.510612711665420e-02 , 4.653266076428588e+03), 
             std::complex<double> (-2.424607329743559e-02 , 2.483715760121476e+03), 
             std::complex<double> (-6.907273854928963e-03 , 1.325708795379122e+03), 
             std::complex<double> (-1.967514156922237e-03 , 7.076234834716391e+02), 
             std::complex<double> (-5.601966533124439e-04 , 3.777321101915988e+02), 
             std::complex<double> (-1.592579215613054e-04 , 2.016797505915737e+02), 
             std::complex<double> (-4.503170580998576e-05 , 1.077654749782195e+02), 
             std::complex<double> (-1.248785715582555e-05 , 5.774033972883079e+01), 
             std::complex<double> (-3.214871876352756e-06 , 3.122497845157856e+01), 
             std::complex<double> (-6.031303456583466e-07 , 1.732451124679044e+01), 
             std::complex<double> (-3.414088864119126e-08 , 9.494141891654081e+00), 
             std::complex<double> (-8.530747958690035e-11 , 3.141735081105732e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (2.749817779605933e+07 , 3.928675851365912e+07), 
             std::complex<double> (3.436446635272976e+07 , -2.454260473454576e+06), 
             std::complex<double> (1.509305120206410e+07 , -1.380729812793133e+07), 
             std::complex<double> (4.886150735896152e+06 , -1.017075992077017e+07), 
             std::complex<double> (1.444319689288647e+06 , -5.906787649253402e+06), 
             std::complex<double> (4.158415059837546e+05 , -3.227952632427697e+06), 
             std::complex<double> (1.188282774014820e+05 , -1.734480382269202e+06), 
             std::complex<double> (3.388271048795941e+04 , -9.275496530209908e+05), 
             std::complex<double> (9.655402838304772e+03 , -4.953533683153779e+05), 
             std::complex<double> (2.750976207793482e+03 , -2.644387718213535e+05), 
             std::complex<double> (7.837574289190582e+02 , -1.411520870983752e+05), 
             std::complex<double> (2.232905547392273e+02 , -7.534177512849862e+04), 
             std::complex<double> (6.361467081224986e+01 , -4.021429968376477e+04), 
             std::complex<double> (1.812356882920825e+01 , -2.146466254756372e+04), 
             std::complex<double> (5.163331720760453e+00 , -1.145690479052291e+04), 
             std::complex<double> (1.471012341135676e+00 , -6.115196951391106e+03), 
             std::complex<double> (4.190854753553419e-01 , -3.264025578102081e+03), 
             std::complex<double> (1.193957562144194e-01 , -1.742194341571693e+03), 
             std::complex<double> (3.401537322103352e-02 , -9.299068794628328e+02), 
             std::complex<double> (9.690848873605763e-03 , -4.963427518708056e+02), 
             std::complex<double> (2.760887128542860e-03 , -2.649242232595692e+02), 
             std::complex<double> (7.865661139762456e-04 , -1.414014294778943e+02), 
             std::complex<double> (2.240892746231011e-04 , -7.546722142548825e+01), 
             std::complex<double> (6.384215836620829e-05 , -4.026856732483116e+01), 
             std::complex<double> (1.818860139193252e-05 , -2.147014018048246e+01), 
             std::complex<double> (5.182454136748546e-06 , -1.141623593110925e+01), 
             std::complex<double> (1.475901739479886e-06 , -6.017089153223023e+00), 
             std::complex<double> (3.923684598711579e-07 , -3.158204500293917e+00), 
             std::complex<double> (4.370748052730614e-08 , -2.095194761876006e+00), 
             std::complex<double> (2.324864606821677e-10 , -2.000392276848578e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 30;
            int m = 2;
            double b = 2.800157732752085e+07;
             std::complex<double> zvec1[] = {
             std::complex<double> (-5.144055897068965e+07 , 1.483120416683455e+07), 
             std::complex<double> (-2.891357282658651e+07 , 2.783833296884530e+07), 
             std::complex<double> (-1.146976901179962e+07 , 2.252779455957355e+07), 
             std::complex<double> (-3.925993641212572e+06 , 1.425937570589937e+07), 
             std::complex<double> (-1.275652508562986e+06 , 8.333534521330995e+06), 
             std::complex<double> (-4.074596447736801e+05 , 4.747237500773826e+06), 
             std::complex<double> (-1.294358191095930e+05 , 2.682345403979623e+06), 
             std::complex<double> (-4.104560386435923e+04 , 1.511700660968927e+06), 
             std::complex<double> (-1.300883645931853e+04 , 8.512576494746523e+05), 
             std::complex<double> (-4.122247536582769e+03 , 4.792294228016091e+05), 
             std::complex<double> (-1.306187534343982e+03 , 2.697678464343942e+05), 
             std::complex<double> (-4.138751497003890e+02 , 1.518537680012092e+05), 
             std::complex<double> (-1.311386597164553e+02 , 8.547859461758003e+04), 
             std::complex<double> (-4.155194443347297e+01 , 4.811583630649425e+04), 
             std::complex<double> (-1.316593432979022e+01 , 2.708434580128943e+04), 
             std::complex<double> (-4.171688168539736e+00 , 1.524574200443283e+04), 
             std::complex<double> (-1.321817841941753e+00 , 8.581808236077301e+03), 
             std::complex<double> (-4.188226631971046e-01 , 4.830690601576561e+03), 
             std::complex<double> (-1.327042533774207e-01 , 2.719194231431160e+03), 
             std::complex<double> (-4.204622760276631e-02 , 1.530640195421175e+03), 
             std::complex<double> (-1.332077651308736e-02 , 8.616123616726908e+02), 
             std::complex<double> (-4.218971060427111e-03 , 4.850309180498789e+02), 
             std::complex<double> (-1.335019918776156e-03 , 2.730774987518769e+02), 
             std::complex<double> (-4.212269856076930e-04 , 1.538113817795715e+02), 
             std::complex<double> (-1.316843559579567e-04 , 8.675151323843581e+01), 
             std::complex<double> (-3.993082059949652e-05 , 4.913607540097872e+01), 
             std::complex<double> (-1.085287442972345e-05 , 2.818825776588787e+01), 
             std::complex<double> (-1.907543362141894e-06 , 1.663323083836115e+01), 
             std::complex<double> (-7.966303417965077e-08 , 9.453485964598679e+00), 
             std::complex<double> (-1.310238991487586e-10 , 3.141632582666812e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (5.164215060255514e+06 , 8.211317249128222e+06), 
             std::complex<double> (7.297653026028864e+06 , 2.766266611335894e+05), 
             std::complex<double> (3.731465958037045e+06 , -2.714572986446508e+06), 
             std::complex<double> (1.383754436514395e+06 , -2.322434851407710e+06), 
             std::complex<double> (4.612003641399817e+05 , -1.471157126686850e+06), 
             std::complex<double> (1.485068127955537e+05 , -8.587395665244827e+05), 
             std::complex<double> (4.729629461902505e+04 , -4.889279951055130e+05), 
             std::complex<double> (1.501036730919981e+04 , -2.762105008620486e+05), 
             std::complex<double> (4.758550376624863e+03 , -1.556559811192289e+05), 
             std::complex<double> (1.508014891090192e+03 , -8.765019862413383e+04), 
             std::complex<double> (4.778463773462889e+02 , -4.934382001082131e+04), 
             std::complex<double> (1.514103848109962e+02 , -2.777657278153364e+04), 
             std::complex<double> (4.797535328602482e+01 , -1.563557244579445e+04), 
             std::complex<double> (1.520124524293161e+01 , -8.801273171341929e+03), 
             std::complex<double> (4.816589773094258e+00 , -4.954229587063105e+03), 
             std::complex<double> (1.526159776592680e+00 , -2.788729528500177e+03), 
             std::complex<double> (4.835710885424361e-01 , -1.569771693928325e+03), 
             std::complex<double> (1.532218412029837e-01 , -8.836216927823823e+02), 
             std::complex<double> (4.854908948110069e-02 , -4.973883846318289e+02), 
             std::complex<double> (1.538302053410261e-02 , -2.799774337151554e+02), 
             std::complex<double> (4.874188439827747e-03 , -1.575957467787032e+02), 
             std::complex<double> (1.544411034443838e-03 , -8.870482806045658e+01), 
             std::complex<double> (4.893536585166794e-04 , -4.992188851806088e+01), 
             std::complex<double> (1.550541743038486e-04 , -2.808331873186118e+01), 
             std::complex<double> (4.913162730551911e-05 , -1.577675284921960e+01), 
             std::complex<double> (1.557232157210188e-05 , -8.825588832667364e+00), 
             std::complex<double> (4.923897018755638e-06 , -4.879445917330541e+00), 
             std::complex<double> (1.346262210116094e-06 , -2.733110063591393e+00), 
             std::complex<double> (1.125315837075026e-07 , -2.042980526122038e+00), 
             std::complex<double> (3.732368722903438e-10 , -2.000114689334910e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 30;
            int m = 2;
            double b = 1.037878524103058e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.965353452355160e+06 , 4.512972094019555e+05), 
             std::complex<double> (-1.339521437976738e+06 , 9.884987256279119e+05), 
             std::complex<double> (-6.963640820183916e+05 , 9.776720747111859e+05), 
             std::complex<double> (-3.126809423792953e+05 , 7.410548299523090e+05), 
             std::complex<double> (-1.312410243006304e+05 , 5.042928468151261e+05), 
             std::complex<double> (-5.352885708089831e+04 , 3.284584768792502e+05), 
             std::complex<double> (-2.157749299970477e+04 , 2.101853565679360e+05), 
             std::complex<double> (-8.656681258064169e+03 , 1.335499565627082e+05), 
             std::complex<double> (-3.466359745779085e+03 , 8.461574890451488e+04), 
             std::complex<double> (-1.386960394429364e+03 , 5.355067201267625e+04), 
             std::complex<double> (-5.547811413513231e+02 , 3.387514431116203e+04), 
             std::complex<double> (-2.218840742315764e+02 , 2.142488030488860e+04), 
             std::complex<double> (-8.873789662753697e+01 , 1.354952573356077e+04), 
             std::complex<double> (-3.548813107626696e+01 , 8.568745148482347e+03), 
             std::complex<double> (-1.419228989154608e+01 , 5.418829625762190e+03), 
             std::complex<double> (-5.675672269839032e+00 , 3.426827336355957e+03), 
             std::complex<double> (-2.269727200130654e+00 , 2.167101895751966e+03), 
             std::complex<double> (-9.076324772316735e-01 , 1.370468624168203e+03), 
             std::complex<double> (-3.629088983542017e-01 , 8.666947820190341e+02), 
             std::complex<double> (-1.450656880485980e-01 , 5.481278145665016e+02), 
             std::complex<double> (-5.794704980437172e-02 , 3.466921854736551e+02), 
             std::complex<double> (-2.310707295735892e-02 , 2.193424610314883e+02), 
             std::complex<double> (-9.174013991361829e-03 , 1.388649183955535e+02), 
             std::complex<double> (-3.601751611839669e-03 , 8.806232967527349e+01), 
             std::complex<double> (-1.372686139969153e-03 , 5.607940296954874e+01), 
             std::complex<double> (-4.801555835849495e-04 , 3.608119261510220e+01), 
             std::complex<double> (-1.264144330401454e-04 , 2.373799431571446e+01), 
             std::complex<double> (-1.350559549059643e-05 , 1.584477852588604e+01), 
             std::complex<double> (-1.982247330964216e-07 , 9.426340415704457e+00), 
             std::complex<double> (-1.134825816149471e-10 , 3.141593461908054e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.275224346846662e+05 , 2.630323535462741e+05), 
             std::complex<double> (2.310564600253989e+05 , 7.129914050357543e+04), 
             std::complex<double> (1.651058471994774e+05 , -5.710179874835252e+04), 
             std::complex<double> (8.396101888307597e+04 , -8.178066962260584e+04), 
             std::complex<double> (3.703756010136355e+04 , -6.633883186510229e+04), 
             std::complex<double> (1.541009775557300e+04 , -4.602326693548385e+04), 
             std::complex<double> (6.261487544679016e+03 , -3.017502912565647e+04), 
             std::complex<double> (2.520064645120270e+03 , -1.935732868379216e+04), 
             std::complex<double> (1.010386012192963e+03 , -1.231134315053274e+04), 
             std::complex<double> (4.044819153259888e+02 , -7.803303706091599e+03), 
             std::complex<double> (1.618248853034352e+02 , -4.939217683880694e+03), 
             std::complex<double> (6.472697082963181e+01 , -3.124645077766205e+03), 
             std::complex<double> (2.588706274368349e+01 , -1.976279628132655e+03), 
             std::complex<double> (1.035292994262468e+01 , -1.249850625740999e+03), 
             std::complex<double> (4.140349244223376e+00 , -7.904101005909387e+02), 
             std::complex<double> (1.655800415150153e+00 , -4.998507095528129e+02), 
             std::complex<double> (6.621828806799700e-01 , -3.161000241038133e+02), 
             std::complex<double> (2.648178182993499e-01 , -1.998963186743309e+02), 
             std::complex<double> (1.059047618454167e-01 , -1.264087758121324e+02), 
             std::complex<double> (4.235286678511035e-02 , -7.993387477234691e+01), 
             std::complex<double> (1.693751070435709e-02 , -5.054030024683711e+01), 
             std::complex<double> (6.773555400403035e-03 , -3.194684384172625e+01), 
             std::complex<double> (2.708878744118172e-03 , -2.018017305479761e+01), 
             std::complex<double> (1.083442860894291e-03 , -1.272570838680152e+01), 
             std::complex<double> (4.335808047196679e-04 , -7.990248136963910e+00), 
             std::complex<double> (1.736642471388445e-04 , -4.964258027946280e+00), 
             std::complex<double> (6.570072990076973e-05 , -3.053743016083243e+00), 
             std::complex<double> (1.304302882610784e-05 , -2.145995627766039e+00), 
             std::complex<double> (3.522151453232783e-07 , -2.002860125401003e+00), 
             std::complex<double> (3.581589429484659e-10 , -2.000002562708235e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 30;
            int m = 2;
            double b = 1.208214796412422e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.326264678299416e+05 , 4.438938330825218e+04), 
             std::complex<double> (-1.771988370881783e+05 , 1.064078944852527e+05), 
             std::complex<double> (-1.092632823369572e+05 , 1.200190004967771e+05), 
             std::complex<double> (-5.894900053435823e+04 , 1.036215723787185e+05), 
             std::complex<double> (-2.955142741499655e+04 , 7.906622848242997e+04), 
             std::complex<double> (-1.427263254153668e+04 , 5.689797406581194e+04), 
             std::complex<double> (-6.769697011325071e+03 , 3.982869873385919e+04), 
             std::complex<double> (-3.183425409766492e+03 , 2.752056002813206e+04), 
             std::complex<double> (-1.490935296312497e+03 , 1.890075075820784e+04), 
             std::complex<double> (-6.969428117749482e+02 , 1.294396431187710e+04), 
             std::complex<double> (-3.254984138179635e+02 , 8.852777516504178e+03), 
             std::complex<double> (-1.519564203745581e+02 , 6.050942491823870e+03), 
             std::complex<double> (-7.092547488342092e+01 , 4.134674697506451e+03), 
             std::complex<double> (-3.310091859981494e+01 , 2.824892302144719e+03), 
             std::complex<double> (-1.544709742818827e+01 , 1.929909114183509e+03), 
             std::complex<double> (-7.208079804683848e+00 , 1.318447454656883e+03), 
             std::complex<double> (-3.363071306277876e+00 , 9.007222998838214e+02), 
             std::complex<double> (-1.568706122456839e+00 , 6.153660422067364e+02), 
             std::complex<double> (-7.313294247014348e-01 , 4.204476388961213e+02), 
             std::complex<double> (-3.405507505142688e-01 , 2.873221889822486e+02), 
             std::complex<double> (-1.581854313205626e-01 , 1.964248666660590e+02), 
             std::complex<double> (-7.307924608434435e-02 , 1.343967860915221e+02), 
             std::complex<double> (-3.335824975853839e-02 , 9.212232859642917e+01), 
             std::complex<double> (-1.481202758961026e-02 , 6.339033385844660e+01), 
             std::complex<double> (-6.140014740625283e-03 , 4.398197357455712e+01), 
             std::complex<double> (-2.091355610712495e-03 , 3.103114862932594e+01), 
             std::complex<double> (-3.927080573219232e-04 , 2.237367633263601e+01), 
             std::complex<double> (-1.855158538082567e-05 , 1.572135579160974e+01), 
             std::complex<double> (-1.035269068760825e-07 , 9.424841298580073e+00), 
             std::complex<double> (-2.611236198126789e-11 , 3.141592668564905e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.053145124619291e+04 , 2.659069957290136e+04), 
             std::complex<double> (2.205655130852242e+04 , 1.174267778359135e+04), 
             std::complex<double> (1.956029415937687e+04 , -1.839207648879896e+03), 
             std::complex<double> (1.241614413667250e+04 , -7.380960904855198e+03), 
             std::complex<double> (6.711473305160343e+03 , -7.724203057249810e+03), 
             std::complex<double> (3.357459168373244e+03 , -6.271163067387643e+03), 
             std::complex<double> (1.618839045492734e+03 , -4.624539412736422e+03), 
             std::complex<double> (7.671084116412120e+02 , -3.271433417594468e+03), 
             std::complex<double> (3.605573723580402e+02 , -2.271184270331180e+03), 
             std::complex<double> (1.688255832535028e+02 , -1.563202500791270e+03), 
             std::complex<double> (7.890952583545085e+01 , -1.071615145499707e+03), 
             std::complex<double> (3.685189994399427e+01 , -7.332517168667483e+02), 
             std::complex<double> (1.720371022955850e+01 , -5.012909265052306e+02), 
             std::complex<double> (8.029816967337558e+00 , -3.425704466343458e+02), 
             std::complex<double> (3.747586198508032e+00 , -2.340594201571213e+02), 
             std::complex<double> (1.748955049976860e+00 , -1.599043562120944e+02), 
             std::complex<double> (8.161969895483653e-01 , -1.092366875668505e+02), 
             std::complex<double> (3.808954824460178e-01 , -7.461928818647574e+01), 
             std::complex<double> (1.777524083034638e-01 , -5.096742055962000e+01), 
             std::complex<double> (8.295184119236429e-02 , -3.480588288135431e+01), 
             std::complex<double> (3.871144926903150e-02 , -2.375966119911232e+01), 
             std::complex<double> (1.806636670447087e-02 , -1.620528252382361e+01), 
             std::complex<double> (8.433120720010314e-03 , -1.103224973368146e+01), 
             std::complex<double> (3.939700145429172e-03 , -7.479675907127952e+00), 
             std::complex<double> (1.844663046666714e-03 , -5.025258284456999e+00), 
             std::complex<double> (8.468920150845316e-04 , -3.327886031225074e+00), 
             std::complex<double> (2.766320167510206e-04 , -2.310264106555836e+00), 
             std::complex<double> (2.374624815065725e-05 , -2.018152676233364e+00), 
             std::complex<double> (2.133110970243484e-07 , -2.000132624439284e+00), 
             std::complex<double> (8.879757402359506e-11 , -2.000000051069655e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 30;
            int m = 2;
            double b = 1.786300538817640e+04;
             std::complex<double> zvec1[] = {
             std::complex<double> (-3.482835827820777e+04 , 5.424190012136590e+03), 
             std::complex<double> (-2.893810661534187e+04 , 1.396062841592816e+04), 
             std::complex<double> (-2.058347397381449e+04 , 1.762381479312181e+04), 
             std::complex<double> (-1.308886772633729e+04 , 1.719314303770886e+04), 
             std::complex<double> (-7.743687979452162e+03 , 1.470640607039742e+04), 
             std::complex<double> (-4.388164051580203e+03 , 1.171684906286218e+04), 
             std::complex<double> (-2.426344947905732e+03 , 8.981582692213115e+03), 
             std::complex<double> (-1.323438244666042e+03 , 6.742424415794953e+03), 
             std::complex<double> (-7.165038457904284e+02 , 5.004717352696811e+03), 
             std::complex<double> (-3.863464424585692e+02 , 3.692369813835811e+03), 
             std::complex<double> (-2.078659897858887e+02 , 2.715269912412166e+03), 
             std::complex<double> (-1.117039754454796e+02 , 1.993242169016990e+03), 
             std::complex<double> (-5.998670510754253e+01 , 1.461843158803807e+03), 
             std::complex<double> (-3.219928979683724e+01 , 1.071587995063758e+03), 
             std::complex<double> (-1.727712165046628e+01 , 7.853239667778385e+02), 
             std::complex<double> (-9.266083302374119e+00 , 5.754785508933257e+02), 
             std::complex<double> (-4.965957745297645e+00 , 4.217138332047496e+02), 
             std::complex<double> (-2.657890615469179e+00 , 3.090776606648884e+02), 
             std::complex<double> (-1.419128066110732e+00 , 2.265977561595959e+02), 
             std::complex<double> (-7.543804277059780e-01 , 1.662318650570458e+02), 
             std::complex<double> (-3.977262414022457e-01 , 1.220913778800418e+02), 
             std::complex<double> (-2.063665757446078e-01 , 8.987045871137030e+01), 
             std::complex<double> (-1.036556446532321e-01 , 6.642826418922323e+01), 
             std::complex<double> (-4.847304122988369e-02 , 4.948435976001259e+01), 
             std::complex<double> (-1.891070221918447e-02 , 3.738281221346241e+01), 
             std::complex<double> (-4.474243501553810e-03 , 2.876820897784674e+01), 
             std::complex<double> (-3.484281138185630e-04 , 2.201954255181188e+01), 
             std::complex<double> (-5.192442862488305e-06 , 1.570831546873151e+01), 
             std::complex<double> (-1.055839973835105e-08 , 9.424778619222671e+00), 
             std::complex<double> (-1.110227986340854e-12 , 3.141592653656908e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.058603073719127e+03 , 3.316169472438535e+03), 
             std::complex<double> (2.484770612913337e+03 , 1.975909551234731e+03), 
             std::complex<double> (2.647356830977340e+03 , 4.126335653275221e+02), 
             std::complex<double> (2.060793469291490e+03 , -5.690649691072831e+02), 
             std::complex<double> (1.356451257143194e+03 , -9.309256348256666e+02), 
             std::complex<double> (8.137644492551511e+02 , -9.340313657531439e+02), 
             std::complex<double> (4.639253995481013e+02 , -7.959868647191480e+02), 
             std::complex<double> (2.572343945037466e+02 , -6.300020229589354e+02), 
             std::complex<double> (1.404993165911721e+02 , -4.806162661507047e+02), 
             std::complex<double> (7.611904812153365e+01 , -3.597397667018288e+02), 
             std::complex<double> (4.105967092526156e+01 , -2.665776486608805e+02), 
             std::complex<double> (2.209633064296089e+01 , -1.964910978792202e+02), 
             std::complex<double> (1.187626013446552e+01 , -1.444182333881764e+02), 
             std::complex<double> (6.378880639037950e+00 , -1.059820124919450e+02), 
             std::complex<double> (3.424886138614437e+00 , -7.770961578069281e+01), 
             std::complex<double> (1.838480190683470e+00 , -5.695132906374673e+01), 
             std::complex<double> (9.868104016774549e-01 , -4.172417812342363e+01), 
             std::complex<double> (5.296525880106261e-01 , -3.055882763227439e+01), 
             std::complex<double> (2.842490722753033e-01 , -2.237209798332090e+01), 
             std::complex<double> (1.525166904702508e-01 , -1.636743313808943e+01), 
             std::complex<double> (8.182625183402431e-02 , -1.195951826587348e+01), 
             std::complex<double> (4.390937177772417e-02 , -8.717992458092597e+00), 
             std::complex<double> (2.357742506580863e-02 , -6.325576041433425e+00), 
             std::complex<double> (1.267582568904524e-02 , -4.547763452767348e+00), 
             std::complex<double> (6.686257427377829e-03 , -3.224132966353317e+00), 
             std::complex<double> (2.679830388945505e-03 , -2.347926640374447e+00), 
             std::complex<double> (3.682868882091529e-04 , -2.032279877905038e+00), 
             std::complex<double> (8.480852568846298e-06 , -2.000589142478385e+00), 
             std::complex<double> (2.486229004620719e-08 , -2.000001561429920e+00), 
             std::complex<double> (4.079761101264092e-12 , -2.000000000246414e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 35;
            int m = 2;
            double b = 3.661126584467638e+11;
             std::complex<double> zvec1[] = {
             std::complex<double> (-6.307945639596056e+11 , 2.444676724836913e+11), 
             std::complex<double> (-2.441316153383444e+11 , 3.428218594922180e+11), 
             std::complex<double> (-6.114672605093970e+10 , 2.015575069779266e+11), 
             std::complex<double> (-1.357100001422721e+10 , 9.829609018404015e+10), 
             std::complex<double> (-2.929756656504599e+09 , 4.601162917818328e+10), 
             std::complex<double> (-6.286919756152569e+08 , 2.134819572875746e+10), 
             std::complex<double> (-1.347356797636413e+08 , 9.886246633575878e+09), 
             std::complex<double> (-2.886734924069148e+07 , 4.576413668275360e+09), 
             std::complex<double> (-6.184511402221180e+06 , 2.118269959097176e+09), 
             std::complex<double> (-1.324946538595911e+06 , 9.804585589535079e+08), 
             std::complex<double> (-2.838508051697291e+05 , 4.538114424290624e+08), 
             std::complex<double> (-6.081093083117034e+04 , 2.100493153366877e+08), 
             std::complex<double> (-1.302786096428326e+04 , 9.722255249376732e+07), 
             std::complex<double> (-2.791030369309058e+03 , 4.500002452740357e+07), 
             std::complex<double> (-5.979377877390295e+02 , 2.082852320253444e+07), 
             std::complex<double> (-1.280994980774351e+02 , 9.640603145718211e+06), 
             std::complex<double> (-2.744345684132097e+01 , 4.462209253800593e+06), 
             std::complex<double> (-5.879360297411511e+00 , 2.065359513521473e+06), 
             std::complex<double> (-1.259566230413336e+00 , 9.559636667460443e+05), 
             std::complex<double> (-2.698431011062096e-01 , 4.424733448020630e+05), 
             std::complex<double> (-5.780967105366942e-02 , 2.048013619154774e+05), 
             std::complex<double> (-1.238477200671828e-02 , 9.479350195048715e+04), 
             std::complex<double> (-2.653224600100998e-03 , 4.387572406301765e+04), 
             std::complex<double> (-5.684019297029851e-04 , 2.030813456136477e+04), 
             std::complex<double> (-1.217665056769685e-04 , 9.399739132995324e+03), 
             std::complex<double> (-2.608425170878246e-05 , 4.350725892754523e+03), 
             std::complex<double> (-5.586125381983436e-06 , 2.013763002060904e+03), 
             std::complex<double> (-1.195437279580542e-06 , 9.320910185793575e+02), 
             std::complex<double> (-2.555646067401665e-07 , 4.314434038506678e+02), 
             std::complex<double> (-5.453656271165292e-08 , 1.997385331680464e+02), 
             std::complex<double> (-1.154526652348691e-08 , 9.254166528861525e+01), 
             std::complex<double> (-2.370098625984520e-09 , 4.302918264819222e+01), 
             std::complex<double> (-4.251346942266734e-10 , 2.030895866731932e+01), 
             std::complex<double> (-3.740047797430902e-11 , 9.774042277542533e+00), 
             std::complex<double> (-2.302410606659064e-13 , 3.143296733062276e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.100157830328794e+11 , 1.206170630375351e+11), 
             std::complex<double> (9.692858446525174e+10 , -3.354341267328455e+10), 
             std::complex<double> (2.865194988151737e+10 , -4.287651965861333e+10), 
             std::complex<double> (6.590631889983819e+09 , -2.341333863403447e+10), 
             std::complex<double> (1.433783511683214e+09 , -1.121308874822120e+10), 
             std::complex<double> (3.081805720490319e+08 , -5.227833010532295e+09), 
             std::complex<double> (6.606983503491054e+07 , -2.423490807918841e+09), 
             std::complex<double> (1.415664517139390e+07 , -1.122099679161702e+09), 
             std::complex<double> (3.032954322743739e+06 , -5.194073422832803e+08), 
             std::complex<double> (6.497710195047724e+05 , -2.404143770039436e+08), 
             std::complex<double> (1.392042312822887e+05 , -1.112775566882604e+08), 
             std::complex<double> (2.982249863297660e+04 , -5.150550056145880e+07), 
             std::complex<double> (6.389038547154360e+03 , -2.383962416005900e+07), 
             std::complex<double> (1.368758921187781e+03 , -1.103430909174908e+07), 
             std::complex<double> (2.932367553968926e+02 , -5.107294193942631e+06), 
             std::complex<double> (6.282172333381487e+01 , -2.363940834492869e+06), 
             std::complex<double> (1.345864313887621e+01 , -1.094163768581928e+06), 
             std::complex<double> (2.883318862290756e+00 , -5.064400659736095e+05), 
             std::complex<double> (6.177089623900379e-01 , -2.344087309163525e+05), 
             std::complex<double> (1.323350083110048e-01 , -1.084974448491311e+05), 
             std::complex<double> (2.835077156968269e-02 , -5.021867355741675e+04), 
             std::complex<double> (6.073701819300884e-03 , -2.324400519756631e+04), 
             std::complex<double> (1.301190012240960e-03 , -1.075862298786252e+04), 
             std::complex<double> (2.787578902710746e-04 , -4.979691139328100e+03), 
             std::complex<double> (5.971813761761080e-05 , -2.304878793599981e+03), 
             std::complex<double> (1.279344065777640e-05 , -1.066826081084737e+03), 
             std::complex<double> (2.740677098808649e-06 , -4.937856257199231e+02), 
             std::complex<double> (5.868590687873311e-07 , -2.285493189933585e+02), 
             std::complex<double> (1.255515936911348e-07 , -1.057805639170392e+02), 
             std::complex<double> (2.685043176250633e-08 , -4.895075262850140e+01), 
             std::complex<double> (5.739176017944376e-09 , -2.263481099888547e+01), 
             std::complex<double> (1.220194697817891e-09 , -1.043002305857490e+01), 
             std::complex<double> (2.567993505249665e-10 , -4.754756830172241e+00), 
             std::complex<double> (3.889406163732241e-11 , -2.392986873160181e+00), 
             std::complex<double> (5.623340971642698e-13 , -2.004221935414860e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 35;
            int m = 2;
            double b = 2.559114409292506e+10;
             std::complex<double> zvec1[] = {
             std::complex<double> (-4.532656657464156e+10 , 1.575704656434800e+10), 
             std::complex<double> (-2.049621759823529e+10 , 2.492388407197852e+10), 
             std::complex<double> (-6.168503546190593e+09 , 1.659324414653366e+10), 
             std::complex<double> (-1.623684372909160e+09 , 8.936090504032839e+09), 
             std::complex<double> (-4.121083284828110e+08 , 4.557084668253875e+09), 
             std::complex<double> (-1.036268271367227e+08 , 2.292148460050887e+09), 
             std::complex<double> (-2.599636046812138e+07 , 1.148934179140789e+09), 
             std::complex<double> (-6.517735830320256e+06 , 5.754008953303035e+08), 
             std::complex<double> (-1.633867102271246e+06 , 2.881053990083175e+08), 
             std::complex<double> (-4.095628915006536e+05 , 1.442475901532158e+08), 
             std::complex<double> (-1.026645329126901e+05 , 7.222038454714096e+07), 
             std::complex<double> (-2.573470883983175e+04 , 3.615842773009486e+07), 
             std::complex<double> (-6.450863159756488e+03 , 1.810334825278041e+07), 
             std::complex<double> (-1.617023527588249e+03 , 9.063755476664023e+06), 
             std::complex<double> (-4.053356846650801e+02 , 4.537926230981850e+06), 
             std::complex<double> (-1.016045938190524e+02 , 2.271991367391491e+06), 
             std::complex<double> (-2.546899750743354e+01 , 1.137511826608874e+06), 
             std::complex<double> (-6.384257109469424e+00 , 5.695149960417598e+05), 
             std::complex<double> (-1.600327696759544e+00 , 2.851375459065801e+05), 
             std::complex<double> (-4.011507501446508e-01 , 1.427590505568059e+05), 
             std::complex<double> (-1.005557260614090e-01 , 7.147479116552553e+04), 
             std::complex<double> (-2.520621157020497e-02 , 3.578509231981450e+04), 
             std::complex<double> (-6.318439949869721e-03 , 1.791642679861770e+04), 
             std::complex<double> (-1.583833232967032e-03 , 8.970170536082336e+03), 
             std::complex<double> (-3.970108999862773e-04 , 4.491073576147814e+03), 
             std::complex<double> (-9.951662711028358e-05 , 2.248538454521619e+03), 
             std::complex<double> (-2.494541209018015e-05 , 1.125779086671986e+03), 
             std::complex<double> (-6.252063242514667e-06 , 5.636595075277912e+02), 
             std::complex<double> (-1.565813725523499e-06 , 2.822432614196582e+02), 
             std::complex<double> (-3.910221927625447e-07 , 1.413845708914203e+02), 
             std::complex<double> (-9.650832960367997e-08 , 7.093527043268038e+01), 
             std::complex<double> (-2.268160142307858e-08 , 3.580815253223169e+01), 
             std::complex<double> (-4.277740800393754e-09 , 1.845658665217416e+01), 
             std::complex<double> (-3.103360619177208e-10 , 9.584244111717663e+00), 
             std::complex<double> (-1.207170778039565e-12 , 3.142091104567167e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (6.468642523996503e+09 , 8.179457630157963e+09), 
             std::complex<double> (6.930861710529076e+09 , -1.364240044619521e+09), 
             std::complex<double> (2.541942546031832e+09 , -2.946425569663490e+09), 
             std::complex<double> (7.032566431126790e+08 , -1.871326819524766e+09), 
             std::complex<double> (1.807424765469372e+08 , -9.911506367359675e+08), 
             std::complex<double> (4.559161239337274e+07 , -5.031957178101393e+08), 
             std::complex<double> (1.144635638627455e+07 , -2.528118747650265e+08), 
             std::complex<double> (2.870365432884742e+06 , -1.266850151831013e+08), 
             std::complex<double> (7.195792879144998e+05 , -6.344090345726652e+07), 
             std::complex<double> (1.803798051115711e+05 , -3.176452618474746e+07), 
             std::complex<double> (4.521568337544602e+04 , -1.590367841802815e+07), 
             std::complex<double> (1.133413185470298e+04 , -7.962479419785627e+06), 
             std::complex<double> (2.841102584052593e+03 , -3.986556656781173e+06), 
             std::complex<double> (7.121728979379262e+02 , -1.995939062421094e+06), 
             std::complex<double> (1.785187930889132e+02 , -9.993015121963143e+05), 
             std::complex<double> (4.474890822538592e+01 , -5.003176151656802e+05), 
             std::complex<double> (1.121710919927146e+01 , -2.504926795770592e+05), 
             std::complex<double> (2.811767799331259e+00 , -1.254134982621653e+05), 
             std::complex<double> (7.048195746723961e-01 , -6.279043986969013e+04), 
             std::complex<double> (1.766755603097578e-01 , -3.143712113531950e+04), 
             std::complex<double> (4.428687029620062e-02 , -1.573953910433104e+04), 
             std::complex<double> (1.110131089866292e-02 , -7.880272771129580e+03), 
             std::complex<double> (2.782764617619278e-03 , -3.945394906844859e+03), 
             std::complex<double> (6.975576671268311e-04 , -1.975330029066225e+03), 
             std::complex<double> (1.748545901521675e-04 , -9.889826714440525e+02), 
             std::complex<double> (4.382943183209879e-05 , -4.951502691566606e+02), 
             std::complex<double> (1.098677789825252e-05 , -2.479034944890813e+02), 
             std::complex<double> (2.754250788960751e-06 , -1.241130577820450e+02), 
             std::complex<double> (6.904662787171891e-07 , -6.213112609393782e+01), 
             std::complex<double> (1.730919441341018e-07 , -3.109062229942667e+01), 
             std::complex<double> (4.339541973531861e-08 , -1.553354716913738e+01), 
             std::complex<double> (1.086946993048000e-08 , -7.715176851804723e+00), 
             std::complex<double> (2.641280332353083e-09 , -3.790763970695156e+00), 
             std::complex<double> (3.582905464904482e-10 , -2.198773295732380e+00), 
             std::complex<double> (3.127939277216378e-12 , -2.001307672244486e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 35;
            int m = 2;
            double b = 2.954353221201707e+09;
             std::complex<double> zvec1[] = {
             std::complex<double> (-5.341962982081656e+09 , 1.683777411943043e+09), 
             std::complex<double> (-2.727166734118958e+09 , 2.928883297586369e+09), 
             std::complex<double> (-9.543688022029916e+08 , 2.166538517214723e+09), 
             std::complex<double> (-2.895872263850663e+08 , 1.271532259487842e+09), 
             std::complex<double> (-8.403115921867932e+07 , 6.974434989486202e+08), 
             std::complex<double> (-2.406651767276453e+07 , 3.751745202953250e+08), 
             std::complex<double> (-6.866764478129898e+06 , 2.006967959498261e+08), 
             std::complex<double> (-1.957151396234491e+06 , 1.071909815411994e+08), 
             std::complex<double> (-5.576525104650599e+05 , 5.722418636706638e+07), 
             std::complex<double> (-1.588784499699595e+05 , 3.054534528695559e+07), 
             std::complex<double> (-4.526427056272518e+04 , 1.630401274721097e+07), 
             std::complex<double> (-1.289564234942270e+04 , 8.702407804586332e+06), 
             std::complex<double> (-3.673918829546373e+03 , 4.644971584112166e+06), 
             std::complex<double> (-1.046684726652748e+03 , 2.479284260645661e+06), 
             std::complex<double> (-2.981962293470992e+02 , 1.323334028430303e+06), 
             std::complex<double> (-8.495488997079961e+01 , 7.063380581079188e+05), 
             std::complex<double> (-2.420330026097238e+01 , 3.770124789740098e+05), 
             std::complex<double> (-6.895420189086758e+00 , 2.012328330707724e+05), 
             std::complex<double> (-1.964475913194190e+00 , 1.074093175463131e+05), 
             std::complex<double> (-5.596704963621225e-01 , 5.733041346353997e+04), 
             std::complex<double> (-1.594475845061644e-01 , 3.060047688934252e+04), 
             std::complex<double> (-4.542589560486287e-02 , 1.633320189125173e+04), 
             std::complex<double> (-1.294161882336601e-02 , 8.717952885833845e+03), 
             std::complex<double> (-3.686974178705422e-03 , 4.653266078116745e+03), 
             std::complex<double> (-1.050364146464201e-03 , 2.483715760551319e+03), 
             std::complex<double> (-2.992148248955254e-04 , 1.325708795536287e+03), 
             std::complex<double> (-8.522587105986485e-05 , 7.076234835417355e+02), 
             std::complex<double> (-2.426612237888168e-05 , 3.777321102216015e+02), 
             std::complex<double> (-6.899186303050835e-06 , 2.016797506027174e+02), 
             std::complex<double> (-1.950868264678864e-06 , 1.077654749813868e+02), 
             std::complex<double> (-5.410102738937334e-07 , 5.774033972928748e+01), 
             std::complex<double> (-1.393652647174156e-07 , 3.122497845142529e+01), 
             std::complex<double> (-2.618504058369163e-08 , 1.732451124668063e+01), 
             std::complex<double> (-1.484093553939854e-09 , 9.494141891641695e+00), 
             std::complex<double> (-3.705382757286793e-12 , 3.141735081105678e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (6.347289997914248e+08 , 9.068399047523047e+08), 
             std::complex<double> (7.932206823482014e+08 , -5.665067366296517e+07), 
             std::complex<double> (3.483866226767799e+08 , -3.187081193295093e+08), 
             std::complex<double> (1.127849849502518e+08 , -2.347674205574897e+08), 
             std::complex<double> (3.333863059277647e+07 , -1.363439222798468e+08), 
             std::complex<double> (9.598696502990393e+06 , -7.450948789216682e+07), 
             std::complex<double> (2.742863696671803e+06 , -4.003628917695790e+07), 
             std::complex<double> (7.821005109659435e+05 , -2.141024281060456e+07), 
             std::complex<double> (2.228716468239017e+05 , -1.143403575026430e+07), 
             std::complex<double> (6.349963924286003e+04 , -6.103930172251933e+06), 
             std::complex<double> (1.809114662010968e+04 , -3.258154912007595e+06), 
             std::complex<double> (5.154122945860079e+03 , -1.739082855803733e+06), 
             std::complex<double> (1.468390949528229e+03 , -9.282499521752106e+05), 
             std::complex<double> (4.183387805658296e+02 , -4.954598777473632e+05), 
             std::complex<double> (1.191830281614660e+02 , -2.644549690019308e+05), 
             std::complex<double> (3.395476025476729e+01 , -1.411545501702359e+05), 
             std::complex<double> (9.673572759458848e+00 , -7.534214968356599e+04), 
             std::complex<double> (2.755961404606907e+00 , -4.021435664363626e+04), 
             std::complex<double> (7.851620498502022e-01 , -2.146467121050541e+04), 
             std::complex<double> (2.236892972004485e-01 , -1.145690610851495e+04), 
             std::complex<double> (6.372807375067199e-02 , -6.115197152241593e+03), 
             std::complex<double> (1.815583903815911e-02 , -3.264025608931459e+03), 
             std::complex<double> (5.172524007385113e-03 , -1.742194346428488e+03), 
             std::complex<double> (1.473633045031648e-03 , -9.299068802857494e+02), 
             std::complex<double> (4.198298378333071e-04 , -4.963427520330058e+02), 
             std::complex<double> (1.196040993535331e-04 , -2.649242233025061e+02), 
             std::complex<double> (3.407214254330199e-05 , -1.414014294957790e+02), 
             std::complex<double> (9.706426295970842e-06 , -7.546722143415607e+01), 
             std::complex<double> (2.765516387984685e-06 , -4.026856732870812e+01), 
             std::complex<double> (7.879854524567858e-07 , -2.147014018196367e+01), 
             std::complex<double> (2.244874108599590e-07 , -1.141623593152252e+01), 
             std::complex<double> (6.392933429834181e-08 , -6.017089153267768e+00), 
             std::complex<double> (1.702141884886936e-08 , -3.158204500251845e+00), 
             std::complex<double> (1.899614157011302e-09 , -2.095194761862603e+00), 
             std::complex<double> (1.010255487720709e-11 , -2.000392276848434e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 35;
            int m = 2;
            double b = 1.105465055323651e+08;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.056174747314652e+08 , 5.462702708340009e+07), 
             std::complex<double> (-1.248796737705579e+08 , 1.090636523230593e+08), 
             std::complex<double> (-5.503743367224923e+07 , 9.532205133049019e+07), 
             std::complex<double> (-2.089176005564509e+07 , 6.451702078131849e+07), 
             std::complex<double> (-7.480760734355396e+06 , 3.988459744282800e+07), 
             std::complex<double> (-2.622597454324994e+06 , 2.388379000683395e+07), 
             std::complex<double> (-9.126059256863295e+05 , 1.414423784180475e+07), 
             std::complex<double> (-3.167435211317990e+05 , 8.344132913198477e+06), 
             std::complex<double> (-1.098350016877171e+05 , 4.915888544543148e+06), 
             std::complex<double> (-3.807483410311422e+04 , 2.894819376252319e+06), 
             std::complex<double> (-1.319739487578971e+04 , 1.704398381237641e+06), 
             std::complex<double> (-4.574273071902612e+03 , 1.003451904734638e+06), 
             std::complex<double> (-1.585441815397109e+03 , 5.907634591168201e+05), 
             std::complex<double> (-5.495112214000010e+02 , 3.477985615529468e+05), 
             std::complex<double> (-1.904592741570209e+02 , 2.047580245603100e+05), 
             std::complex<double> (-6.601268526257599e+01 , 1.205462574670654e+05), 
             std::complex<double> (-2.287981408831644e+01 , 7.096862847584837e+04), 
             std::complex<double> (-7.930076614012533e+00 , 4.178102148884641e+04), 
             std::complex<double> (-2.748539587332887e+00 , 2.459754074270638e+04), 
             std::complex<double> (-9.526342175204983e-01 , 1.448119299389701e+04), 
             std::complex<double> (-3.301789648667092e-01 , 8.525445041452595e+03), 
             std::complex<double> (-1.144379758322808e-01 , 5.019147604508195e+03), 
             std::complex<double> (-3.966295597230947e-02 , 2.954904361447640e+03), 
             std::complex<double> (-1.374630778028932e-02 , 1.739636397369029e+03), 
             std::complex<double> (-4.763791474181150e-03 , 1.024184371975565e+03), 
             std::complex<double> (-1.650547643809279e-03 , 6.029913507360656e+02), 
             std::complex<double> (-5.715417791689269e-04 , 3.550441246139415e+02), 
             std::complex<double> (-1.975819291687172e-04 , 2.091048389162337e+02), 
             std::complex<double> (-6.797529705010001e-05 , 1.232436197132160e+02), 
             std::complex<double> (-2.305139212584215e-05 , 7.279170938651387e+01), 
             std::complex<double> (-7.475219227090609e-06 , 4.325330557046042e+01), 
             std::complex<double> (-2.077374999214285e-06 , 2.612498896489193e+01), 
             std::complex<double> (-3.245917253521437e-07 , 1.621892037029814e+01), 
             std::complex<double> (-9.678377189543531e-09 , 9.436111337664489e+00), 
             std::complex<double> (-1.084626152091708e-11 , 3.141603672944733e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.766202043514169e+07 , 3.089397823327621e+07), 
             std::complex<double> (2.756943886827726e+07 , 3.744846708951503e+06), 
             std::complex<double> (1.604085936285369e+07 , -9.260110219198819e+06), 
             std::complex<double> (6.698117833047478e+06 , -9.257932532015027e+06), 
             std::complex<double> (2.479126603495896e+06 , -6.376391772078541e+06), 
             std::complex<double> (8.791698827751188e+05 , -3.954996365712693e+06), 
             std::complex<double> (3.071525929525471e+05 , -2.370329429024518e+06), 
             std::complex<double> (1.067525526867399e+05 , -1.404091756706270e+06), 
             std::complex<double> (3.703557845169985e+04 , -8.283875140501200e+05), 
             std::complex<double> (1.284069206953309e+04 , -4.880527279064512e+05), 
             std::complex<double> (4.451061486500535e+03 , -2.874024415211375e+05), 
             std::complex<double> (1.542787546046310e+03 , -1.692160567828360e+05), 
             std::complex<double> (5.347334423072297e+02 , -9.962481547717983e+04), 
             std::complex<double> (1.853380763359534e+02 , -5.865226256159375e+04), 
             std::complex<double> (6.423779123409220e+01 , -3.453019157957943e+04), 
             std::complex<double> (2.226465997820563e+01 , -2.032881941477426e+04), 
             std::complex<double> (7.716872548092494e+00 , -1.196809323964989e+04), 
             std::complex<double> (2.674646787812702e+00 , -7.045918922174672e+03), 
             std::complex<double> (9.270248294013168e-01 , -4.148110026379695e+03), 
             std::complex<double> (3.213039751237012e-01 , -2.442096690521033e+03), 
             std::complex<double> (1.113629512680212e-01 , -1.437723499166115e+03), 
             std::complex<double> (3.859801743899460e-02 , -8.464234636198585e+02), 
             std::complex<double> (1.337790582194390e-02 , -4.983098398334885e+02), 
             std::complex<double> (4.636698951696993e-03 , -2.933659035918621e+02), 
             std::complex<double> (1.607036093282423e-03 , -1.727090945446318e+02), 
             std::complex<double> (5.569787072250253e-04 , -1.016734350860752e+02), 
             std::complex<double> (1.930387493218836e-04 , -5.984962819894576e+01), 
             std::complex<double> (6.690226720461467e-05 , -3.522124878051003e+01), 
             std::complex<double> (2.318794312521300e-05 , -2.071229126928188e+01), 
             std::complex<double> (8.038343157591242e-06 , -1.215416622295455e+01), 
             std::complex<double> (2.788118524494845e-06 , -7.088712697259929e+00), 
             std::complex<double> (9.590905579726293e-07 , -4.077236901379703e+00), 
             std::complex<double> (2.527622583919645e-07 , -2.448381223348704e+00), 
             std::complex<double> (1.492360183542152e-08 , -2.018319335581660e+00), 
             std::complex<double> (3.211729867164913e-11 , -2.000032837342118e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 35;
            int m = 2;
            double b = 1.026180134646449e+07;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.943201088908280e+07 , 4.462104438955036e+06), 
             std::complex<double> (-1.324423101045890e+07 , 9.773569241784515e+06), 
             std::complex<double> (-6.885150541906087e+06 , 9.666523038062979e+06), 
             std::complex<double> (-3.091565761846407e+06 , 7.327020753586905e+06), 
             std::complex<double> (-1.297617476472537e+06 , 4.986087403947366e+06), 
             std::complex<double> (-5.292550891860635e+05 , 3.247562768508746e+06), 
             std::complex<double> (-2.133428312601008e+05 , 2.078162647894131e+06), 
             std::complex<double> (-8.559107833917180e+04 , 1.320446559479368e+06), 
             std::complex<double> (-3.427288907645283e+04 , 8.366200707557073e+05), 
             std::complex<double> (-1.371327481963341e+04 , 5.294707805626044e+05), 
             std::complex<double> (-5.485281225019179e+03 , 3.349332149148736e+05), 
             std::complex<double> (-2.193832853418786e+03 , 2.118338930282334e+05), 
             std::complex<double> (-8.773784592422566e+02 , 1.339680039582493e+05), 
             std::complex<double> (-3.508826890904584e+02 , 8.472158903395882e+04), 
             std::complex<double> (-1.403245122891171e+02 , 5.357745107243057e+04), 
             std::complex<double> (-5.611819952818349e+01 , 3.388191809148429e+04), 
             std::complex<double> (-2.244259932000173e+01 , 2.142659359086906e+04), 
             std::complex<double> (-8.975154577868876e+00 , 1.354995907587761e+04), 
             std::complex<double> (-3.589301615805944e+00 , 8.568854764721897e+03), 
             std::complex<double> (-1.435412728945383e+00 , 5.418857360447386e+03), 
             std::complex<double> (-5.740383691926080e-01 , 3.426834357187925e+03), 
             std::complex<double> (-2.295608086877896e-01 , 2.167103674968478e+03), 
             std::complex<double> (-9.179840456822293e-02 , 1.370469075897189e+03), 
             std::complex<double> (-3.670471483280777e-02 , 8.666948967738975e+02), 
             std::complex<double> (-1.467174550828827e-02 , 5.481278434769908e+02), 
             std::complex<double> (-5.860460438223919e-03 , 3.466921925043665e+02), 
             std::complex<double> (-2.336793980118080e-03 , 2.193424625856353e+02), 
             std::complex<double> (-9.276987187545067e-04 , 1.388649186711725e+02), 
             std::complex<double> (-3.642018135552275e-04 , 8.806232969321276e+01), 
             std::complex<double> (-1.388090773146201e-04 , 5.607940295433601e+01), 
             std::complex<double> (-4.856440384813166e-05 , 3.608119261180187e+01), 
             std::complex<double> (-1.279071290789698e-05 , 2.373799432083292e+01), 
             std::complex<double> (-1.367146050728631e-06 , 1.584477852780033e+01), 
             std::complex<double> (-2.007386092974969e-08 , 9.426340415767395e+00), 
             std::complex<double> (-1.146809910772601e-11 , 3.141593461908137e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.260850738289951e+06 , 2.600675999295314e+06), 
             std::complex<double> (2.284521212948746e+06 , 7.049549601707552e+05), 
             std::complex<double> (1.632448665698916e+06 , -5.645817960072458e+05), 
             std::complex<double> (8.301465719983791e+05 , -8.085888386401753e+05), 
             std::complex<double> (3.662009332128101e+05 , -6.559109776654616e+05), 
             std::complex<double> (1.523640362045052e+05 , -4.550451848594168e+05), 
             std::complex<double> (6.190911451065743e+04 , -2.983491313779661e+05), 
             std::complex<double> (2.491659811687785e+04 , -1.913914376911295e+05), 
             std::complex<double> (9.989974618093216e+03 , -1.217257665934458e+05), 
             std::complex<double> (3.999228015927376e+03 , -7.715349340728085e+04), 
             std::complex<double> (1.600008817641251e+03 , -4.883545738235494e+04), 
             std::complex<double> (6.399740708396686e+02 , -3.089426058849887e+04), 
             std::complex<double> (2.559528237526280e+02 , -1.954004478811741e+04), 
             std::complex<double> (1.023624149738926e+02 , -1.235763598569351e+04), 
             std::complex<double> (4.093684686792244e+01 , -7.815019782804298e+03), 
             std::complex<double> (1.637139220941916e+01 , -4.942181589875851e+03), 
             std::complex<double> (6.547202168239145e+00 , -3.125394754837007e+03), 
             std::complex<double> (2.618335341756563e+00 , -1.976469237866021e+03), 
             std::complex<double> (1.047114701994884e+00 , -1.249898583026991e+03), 
             std::complex<double> (4.187577693442188e-01 , -7.904222315637359e+02), 
             std::complex<double> (1.674679348265500e-01 , -4.998537788182201e+02), 
             std::complex<double> (6.697319346403730e-02 , -3.161008010446330e+02), 
             std::complex<double> (2.678375967902768e-02 , -1.998965156610936e+02), 
             std::complex<double> (1.071132639059350e-02 , -1.264088259248508e+02), 
             std::complex<double> (4.283636765708093e-03 , -7.993388755917833e+01), 
             std::complex<double> (1.713062928279200e-03 , -5.054030348770195e+01), 
             std::complex<double> (6.850497187306558e-04 , -3.194684462764545e+01), 
             std::complex<double> (2.739465468429037e-04 , -2.018017322446559e+01), 
             std::complex<double> (1.095573762818831e-04 , -1.272570841597254e+01), 
             std::complex<double> (4.383957780249302e-05 , -7.990248136320462e+00), 
             std::complex<double> (1.756033196474276e-05 , -4.964258023328600e+00), 
             std::complex<double> (6.645886263047244e-06 , -3.053743015777302e+00), 
             std::complex<double> (1.320107146649731e-06 , -2.145995629049842e+00), 
             std::complex<double> (3.566779459360094e-08 , -2.002860125498104e+00), 
             std::complex<double> (3.621436857082297e-11 , -2.000002562708464e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 35;
            int m = 2;
            double b = 1.704319355157387e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-3.266436830905796e+06 , 6.604972004282690e+05), 
             std::complex<double> (-2.411959605195913e+06 , 1.543729023197506e+06), 
             std::complex<double> (-1.416124177854151e+06 , 1.676128288857639e+06), 
             std::complex<double> (-7.238804011981136e+05 , 1.391828324270356e+06), 
             std::complex<double> (-3.442101501833206e+05 , 1.025588458781179e+06), 
             std::complex<double> (-1.580761986944787e+05 , 7.158700522845953e+05), 
             std::complex<double> (-7.143760351525380e+04 , 4.876310249975937e+05), 
             std::complex<double> (-3.204962574611533e+04 , 3.285442299409541e+05), 
             std::complex<double> (-1.433169683430693e+04 , 2.202776876120392e+05), 
             std::complex<double> (-6.399356258240637e+03 , 1.473662316202510e+05), 
             std::complex<double> (-2.855557896216299e+03 , 9.849221363013831e+04), 
             std::complex<double> (-1.273852214445582e+03 , 6.579863408259478e+04), 
             std::complex<double> (-5.681865051762112e+02 , 4.394885717406782e+04), 
             std::complex<double> (-2.534183766365905e+02 , 2.935220376339105e+04), 
             std::complex<double> (-1.130250598934417e+02 , 1.960275407182156e+04), 
             std::complex<double> (-5.040890277960460e+01 , 1.309139801661447e+04), 
             std::complex<double> (-2.248218509340209e+01 , 8.742823437044553e+03), 
             std::complex<double> (-1.002697078424605e+01 , 5.838698968769202e+03), 
             std::complex<double> (-4.471995107153908e+00 , 3.899242174934777e+03), 
             std::complex<double> (-1.994498201724353e+00 , 2.604023865623547e+03), 
             std::complex<double> (-8.895428644675663e-01 , 1.739048307981923e+03), 
             std::complex<double> (-3.967302637594754e-01 , 1.161402907110942e+03), 
             std::complex<double> (-1.769264974936149e-01 , 7.756477565261306e+02), 
             std::complex<double> (-7.888329279436990e-02 , 5.180473699593130e+02), 
             std::complex<double> (-3.514755206612464e-02 , 3.460402986840670e+02), 
             std::complex<double> (-1.563575324528675e-02 , 2.312070741801582e+02), 
             std::complex<double> (-6.929701165283079e-03 , 1.545747356354551e+02), 
             std::complex<double> (-3.044194375228442e-03 , 1.034822536994857e+02), 
             std::complex<double> (-1.309432345454434e-03 , 6.948916820219362e+01), 
             std::complex<double> (-5.341526178979395e-04 , 4.698181902932901e+01), 
             std::complex<double> (-1.873934464100954e-04 , 3.223630007369167e+01), 
             std::complex<double> (-4.082279616871332e-05 , 2.264983927559361e+01), 
             std::complex<double> (-2.586628539483479e-06 , 1.573834971462750e+01), 
             std::complex<double> (-1.955685403185753e-08 , 9.424967102185864e+00), 
             std::complex<double> (-6.368901266460307e-12 , 3.141592710621642e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.655836959451874e+05 , 3.926994172378575e+05), 
             std::complex<double> (3.329832835416716e+05 , 1.535707370457346e+05), 
             std::complex<double> (2.774415507888164e+05 , -4.698813618638822e+04), 
             std::complex<double> (1.648851107990223e+05 , -1.156371814621338e+05), 
             std::complex<double> (8.382876594338447e+04 , -1.108182694574222e+05), 
             std::complex<double> (3.966394762505158e+04 , -8.543262918067907e+04), 
             std::complex<double> (1.816516281278227e+04 , -6.066685587765294e+04), 
             std::complex<double> (8.198140000003805e+03 , -4.162015613630528e+04), 
             std::complex<double> (3.675706145542460e+03 , -2.812816861208291e+04), 
             std::complex<double> (1.643208637822709e+03 , -1.888446893604327e+04), 
             std::complex<double> (7.336285069588095e+02 , -1.264131162859335e+04), 
             std::complex<double> (3.273452762380070e+02 , -8.451068401704257e+03), 
             std::complex<double> (1.460235569378905e+02 , -5.646483454555449e+03), 
             std::complex<double> (6.513128435425182e+01 , -3.771652683876342e+03), 
             std::complex<double> (2.904921319977776e+01 , -2.519039694133170e+03), 
             std::complex<double> (1.295596346040507e+01 , -1.682348085982877e+03), 
             std::complex<double> (5.778319146702528e+00 , -1.123535025101075e+03), 
             std::complex<double> (2.577107251309330e+00 , -7.503307148065359e+02), 
             std::complex<double> (1.149379556667323e+00 , -5.010908039082480e+02), 
             std::complex<double> (5.126187133708107e-01 , -3.346403700442774e+02), 
             std::complex<double> (2.286271532513289e-01 , -2.234795394017029e+02), 
             std::complex<double> (1.019698875027650e-01 , -1.492424563745788e+02), 
             std::complex<double> (4.548203389662622e-02 , -9.966359010821046e+01), 
             std::complex<double> (2.028814983846166e-02 , -6.655141504778342e+01), 
             std::complex<double> (9.050727051621337e-03 , -4.443504906707873e+01), 
             std::complex<double> (4.038075455994924e-03 , -2.966035142532295e+01), 
             std::complex<double> (1.801978135170209e-03 , -1.978615854203796e+01), 
             std::complex<double> (8.043635351767344e-04 , -1.318088350773374e+01), 
             std::complex<double> (3.592587556446812e-04 , -8.752786701384906e+00), 
             std::complex<double> (1.607299490014284e-04 , -5.769840505467695e+00), 
             std::complex<double> (7.146320116795132e-05 , -3.749468885020797e+00), 
             std::complex<double> (2.575742480466152e-05 , -2.484391952530531e+00), 
             std::complex<double> (3.042382407782342e-06 , -2.038362574123834e+00), 
             std::complex<double> (3.858186595723820e-08 , -2.000380685542973e+00), 
             std::complex<double> (2.117179035236950e-11 , -2.000000190222287e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 35;
            int m = 2;
            double b = 2.327278926058514e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-4.513323494576385e+05 , 7.740515564565791e+04), 
             std::complex<double> (-3.612667648042175e+05 , 1.932192496435580e+05), 
             std::complex<double> (-2.414049608313716e+05 , 2.321257309013744e+05), 
             std::complex<double> (-1.426600943555093e+05 , 2.143110017543970e+05), 
             std::complex<double> (-7.831239808870041e+04 , 1.739305084404605e+05), 
             std::complex<double> (-4.126608256876442e+04 , 1.321730924609136e+05), 
             std::complex<double> (-2.127833217778567e+04 , 9.712546710446510e+04), 
             std::complex<double> (-1.084948964241150e+04 , 7.016460937004353e+04), 
             std::complex<double> (-5.500374696087707e+03 , 5.025214270944060e+04), 
             std::complex<double> (-2.780433327627710e+03 , 3.583417806395017e+04), 
             std::complex<double> (-1.403440529048081e+03 , 2.549674913579196e+04), 
             std::complex<double> (-7.078688011863339e+02 , 1.812134185932561e+04), 
             std::complex<double> (-3.569015959533069e+02 , 1.287220836337846e+04), 
             std::complex<double> (-1.799128959698423e+02 , 9.140995109593387e+03), 
             std::complex<double> (-9.068512975840655e+01 , 6.490414203806098e+03), 
             std::complex<double> (-4.570803223460333e+01 , 4.608086298860901e+03), 
             std::complex<double> (-2.303791519735180e+01 , 3.271551325124475e+03), 
             std::complex<double> (-1.161147458870107e+01 , 2.322630950314166e+03), 
             std::complex<double> (-5.852053839419206e+00 , 1.648941525013890e+03), 
             std::complex<double> (-2.948981065267999e+00 , 1.170667273800167e+03), 
             std::complex<double> (-1.485670470621111e+00 , 8.311339325152760e+02), 
             std::complex<double> (-7.481085408723449e-01 , 5.901038999373613e+02), 
             std::complex<double> (-3.763950579318537e-01 , 4.190117337688313e+02), 
             std::complex<double> (-1.891136622768643e-01 , 2.975801694023966e+02), 
             std::complex<double> (-9.479577018048617e-02 , 2.114174611673811e+02), 
             std::complex<double> (-4.731220309300446e-02 , 1.503119887833896e+02), 
             std::complex<double> (-2.340465067149067e-02 , 1.070223269491178e+02), 
             std::complex<double> (-1.135691479094760e-02 , 7.641971535966709e+01), 
             std::complex<double> (-5.273377763145217e-03 , 5.488161913371739e+01), 
             std::complex<double> (-2.190474830861128e-03 , 3.986049130139978e+01), 
             std::complex<double> (-6.596869434928744e-04 , 2.951643804388525e+01), 
             std::complex<double> (-8.423457719383459e-05 , 2.210476863340843e+01), 
             std::complex<double> (-2.212380502593744e-06 , 1.571027947537927e+01), 
             std::complex<double> (-7.270203109281250e-09 , 9.424784640092799e+00), 
             std::complex<double> (-1.147655977277364e-12 , 3.141592654602734e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.657992182296191e+04 , 4.691523967763896e+04), 
             std::complex<double> (3.705423249236587e+04 , 2.473167433472940e+04), 
             std::complex<double> (3.642341077572560e+04 , 1.428112417648422e+03), 
             std::complex<double> (2.587127219683498e+04 , -1.082164882573312e+04), 
             std::complex<double> (1.556441930584309e+04 , -1.378017291449684e+04), 
             std::complex<double> (8.588315757585746e+03 , -1.241325300187122e+04), 
             std::complex<double> (4.532506726326130e+03 , -9.847864643971843e+03), 
             std::complex<double> (2.338271143690269e+03 , -7.380112805668017e+03), 
             std::complex<double> (1.192451907721230e+03 , -5.381933708418713e+03), 
             std::complex<double> (6.045787921756873e+02 , -3.872443502705089e+03), 
             std::complex<double> (3.056228914503818e+02 , -2.767759052015589e+03), 
             std::complex<double> (1.542671600241132e+02 , -1.971589314170016e+03), 
             std::complex<double> (7.780992143688849e+01 , -1.402081926875408e+03), 
             std::complex<double> (3.923117743489036e+01 , -9.962361849428841e+02), 
             std::complex<double> (1.977624407276434e+01 , -7.075642655925167e+02), 
             std::complex<double> (9.968152106637827e+00 , -5.024304930434237e+02), 
             std::complex<double> (5.024243941911435e+00 , -3.567290819880934e+02), 
             std::complex<double> (2.532400529909814e+00 , -2.532654560397969e+02), 
             std::complex<double> (1.276461974793216e+00 , -1.798037090348335e+02), 
             std::complex<double> (6.434118098953416e-01 , -1.276468109257120e+02), 
             std::complex<double> (3.243110432266079e-01 , -9.061663872455640e+01), 
             std::complex<double> (1.634571137646655e-01 , -6.432563381383615e+01), 
             std::complex<double> (8.236802912664512e-02 , -4.565822907538877e+01), 
             std::complex<double> (4.149020431188256e-02 , -3.240211424656403e+01), 
             std::complex<double> (2.089018079105966e-02 , -2.298620020568147e+01), 
             std::complex<double> (1.051646305523476e-02 , -1.629444707264225e+01), 
             std::complex<double> (5.296714115734110e-03 , -1.153353334917547e+01), 
             std::complex<double> (2.671521607649165e-03 , -8.138643780094219e+00), 
             std::complex<double> (1.351487635853184e-03 , -5.706374378403558e+00), 
             std::complex<double> (6.852382800892796e-04 , -3.951012256872180e+00), 
             std::complex<double> (3.167134776983506e-04 , -2.721909819451954e+00), 
             std::complex<double> (7.349704153853667e-05 , -2.110616073004013e+00), 
             std::complex<double> (3.244285388698720e-06 , -2.003532278419438e+00), 
             std::complex<double> (1.609167353623791e-08 , -2.000014967696955e+00), 
             std::complex<double> (4.083791125352976e-12 , -2.000000003593889e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 35;
            int m = 2;
            double b = 1.368346306979379e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.661284134218509e+05 , 4.344239887318401e+04), 
             std::complex<double> (-2.173073947926793e+05 , 1.102281011511850e+05), 
             std::complex<double> (-1.500998999121676e+05 , 1.359428229888415e+05), 
             std::complex<double> (-9.219094413570182e+04 , 1.291876425709090e+05), 
             std::complex<double> (-5.262986921484789e+04 , 1.077488849638818e+05), 
             std::complex<double> (-2.880484951293436e+04 , 8.390750971169297e+04), 
             std::complex<double> (-1.540348104515853e+04 , 6.301848662782524e+04), 
             std::complex<double> (-8.135164670759678e+03 , 4.643825073488006e+04), 
             std::complex<double> (-4.268327846531819e+03 , 3.388129044376720e+04), 
             std::complex<double> (-2.231796513383334e+03 , 2.459152418572268e+04), 
             std::complex<double> (-1.164862398718876e+03 , 1.780051409800733e+04), 
             std::complex<double> (-6.074223063401635e+02 , 1.286664295006284e+04), 
             std::complex<double> (-3.165850933192811e+02 , 9.293470364555118e+03), 
             std::complex<double> (-1.649536660020866e+02 , 6.710026444260639e+03), 
             std::complex<double> (-8.592946523961575e+01 , 4.843780845710698e+03), 
             std::complex<double> (-4.475519084405402e+01 , 3.496235736822668e+03), 
             std::complex<double> (-2.330593469342423e+01 , 2.523452346503214e+03), 
             std::complex<double> (-1.213412528275187e+01 , 1.821293596956446e+03), 
             std::complex<double> (-6.316563554236127e+00 , 1.314507486690182e+03), 
             std::complex<double> (-3.287979992668263e+00 , 9.487486246454343e+02), 
             std::complex<double> (-1.711724586449948e+00 , 6.847831117660782e+02), 
             std::complex<double> (-8.914373654139039e-01 , 4.942921108350927e+02), 
             std::complex<double> (-4.644591221152064e-01 , 3.568378829509836e+02), 
             std::complex<double> (-2.420156555284626e-01 , 2.576724769350066e+02), 
             std::complex<double> (-1.259402454237831e-01 , 1.861559728624649e+02), 
             std::complex<double> (-6.524140136727251e-02 , 1.346152047834656e+02), 
             std::complex<double> (-3.343082639174252e-02 , 9.752092994533440e+01), 
             std::complex<double> (-1.672157514354674e-02 , 7.089558394297437e+01), 
             std::complex<double> (-7.915074827835652e-03 , 5.188817335723763e+01), 
             std::complex<double> (-3.254043737308055e-03 , 3.846246990671016e+01), 
             std::complex<double> (-8.945782700727192e-04 , 2.907153492536575e+01), 
             std::complex<double> (-9.181508576232362e-05 , 2.204917271557439e+01), 
             std::complex<double> (-1.864347368121481e-06 , 1.570887975663578e+01), 
             std::complex<double> (-5.118614276716521e-09 , 9.424780074695237e+00), 
             std::complex<double> (-8.029598596237081e-13 , 3.141592653851281e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (8.873014296195062e+03 , 2.645366341243227e+04), 
             std::complex<double> (2.035605148911698e+04 , 1.490271654851431e+04), 
             std::complex<double> (2.088160936904645e+04 , 2.072259002925070e+03), 
             std::complex<double> (1.556276422441292e+04 , -5.350970799542912e+03), 
             std::complex<double> (9.811990846597333e+03 , -7.647595137664073e+03), 
             std::complex<double> (5.654373298771113e+03 , -7.264929713015909e+03), 
             std::complex<double> (3.105952974994619e+03 , -5.973994766473129e+03), 
             std::complex<double> (1.663450219856786e+03 , -4.602211588076935e+03), 
             std::complex<double> (8.791331683858780e+02 , -3.434043099568624e+03), 
             std::complex<double> (4.614085710234661e+02 , -2.521373833784440e+03), 
             std::complex<double> (2.412963932733338e+02 , -1.835998798462590e+03), 
             std::complex<double> (1.259532722173400e+02 , -1.331213816821988e+03), 
             std::complex<double> (6.568340717920867e+01 , -9.630695634421785e+02), 
             std::complex<double> (3.423629685263381e+01 , -6.959294146243286e+02), 
             std::complex<double> (1.783994114605263e+01 , -5.025871221789065e+02), 
             std::complex<double> (9.294217762167191e+00 , -3.628453650249633e+02), 
             std::complex<double> (4.841264833498431e+00 , -2.619152047595795e+02), 
             std::complex<double> (2.521291169347260e+00 , -1.890434270159062e+02), 
             std::complex<double> (1.312735026613885e+00 , -1.364394394263236e+02), 
             std::complex<double> (6.832947811272498e-01 , -9.846937060320541e+01), 
             std::complex<double> (3.555923785838337e-01 , -7.106279792856250e+01), 
             std::complex<double> (1.850600723880460e-01 , -5.128034405528786e+01), 
             std::complex<double> (9.636167974674677e-02 , -3.699987287052112e+01), 
             std::complex<double> (5.024118718300388e-02 , -2.668928681309445e+01), 
             std::complex<double> (2.624658124979181e-02 , -1.924230801677728e+01), 
             std::complex<double> (1.374180896525411e-02 , -1.385980443983916e+01), 
             std::complex<double> (7.210358199452321e-03 , -9.963940963742065e+00), 
             std::complex<double> (3.793251273097910e-03 , -7.135934354100351e+00), 
             std::complex<double> (2.004695599067881e-03 , -5.071173978093083e+00), 
             std::complex<double> (1.058010046511066e-03 , -3.554451988231489e+00), 
             std::complex<double> (4.736175536187006e-04 , -2.510347587968260e+00), 
             std::complex<double> (8.787189316466357e-05 , -2.061278441304192e+00), 
             std::complex<double> (2.868426589209481e-06 , -2.001467651128709e+00), 
             std::complex<double> (1.153600932097999e-08 , -2.000004878740090e+00), 
             std::complex<double> (2.795693697430165e-12 , -2.000000000944322e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 8.134736423706821e+11;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.440809643864922e+12 , 5.008741310420167e+11), 
             std::complex<double> (-6.515196321257891e+11 , 7.922632408784691e+11), 
             std::complex<double> (-1.960801373321361e+11 , 5.274546032541372e+11), 
             std::complex<double> (-5.161255926383042e+10 , 2.840542831706465e+11), 
             std::complex<double> (-1.309981538565485e+10 , 1.448574651533858e+11), 
             std::complex<double> (-3.294018128354926e+09 , 7.286123472575073e+10), 
             std::complex<double> (-8.263543820442933e+08 , 3.652152745333930e+10), 
             std::complex<double> (-2.071812928383177e+08 , 1.829044690030399e+10), 
             std::complex<double> (-5.193624110930210e+07 , 9.158095764230179e+09), 
             std::complex<double> (-1.301890294393754e+07 , 4.585242931673070e+09), 
             std::complex<double> (-3.263429464993271e+06 , 2.295691785310366e+09), 
             std::complex<double> (-8.180372188950301e+05 , 1.149379168079005e+09), 
             std::complex<double> (-2.050556062820008e+05 , 5.754567513124533e+08), 
             std::complex<double> (-5.140083486249537e+04 , 2.881124092881241e+08), 
             std::complex<double> (-1.288453349035496e+04 , 1.442484699582121e+08), 
             std::complex<double> (-3.229737823183041e+03 , 7.222049496261081e+07), 
             std::complex<double> (-8.095915196311006e+02 , 3.615844158660392e+07), 
             std::complex<double> (-2.029386866901165e+02 , 1.810334999136597e+07), 
             std::complex<double> (-5.087028352358930e+01 , 9.063755694635114e+06), 
             std::complex<double> (-1.275158900550609e+01 , 4.537926258220419e+06), 
             std::complex<double> (-3.196435182078091e+00 , 2.271991370749227e+06), 
             std::complex<double> (-8.012534376710599e-01 , 1.137511826999906e+06), 
             std::complex<double> (-2.008530083475540e-01 , 5.695149960757114e+05), 
             std::complex<double> (-5.034941129638724e-02 , 2.851375459029769e+05), 
             std::complex<double> (-1.262183255250669e-02 , 1.427590505521357e+05), 
             std::complex<double> (-3.164206846249448e-03 , 7.147479116274066e+04), 
             std::complex<double> (-7.932941996121635e-04 , 3.578509231837100e+04), 
             std::complex<double> (-1.989274626938337e-04 , 1.791642679792211e+04), 
             std::complex<double> (-4.991617688647164e-05 , 8.970170535750925e+03), 
             std::complex<double> (-1.254209880189170e-05 , 4.491073575979819e+03), 
             std::complex<double> (-3.157576840181941e-06 , 2.248538454431801e+03), 
             std::complex<double> (-7.959306737183112e-07 , 1.125779086624337e+03), 
             std::complex<double> (-2.002430377521418e-07 , 5.636595075036817e+02), 
             std::complex<double> (-5.029126134722351e-08 , 2.822432614075452e+02), 
             std::complex<double> (-1.270799281646027e-08 , 1.413845708853201e+02), 
             std::complex<double> (-3.228987431756206e-09 , 7.093527042970345e+01), 
             std::complex<double> (-7.957795172405493e-10 , 3.580815253095030e+01), 
             std::complex<double> (-1.633323393295353e-10 , 1.845658665182304e+01), 
             std::complex<double> (-1.340243786813576e-11 , 9.584244111693033e+00), 
             std::complex<double> (-5.701085395132057e-14 , 3.142091104567096e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (2.056207481535771e+11 , 2.600029590914943e+11), 
             std::complex<double> (2.203134529581333e+11 , -4.336552182703562e+10), 
             std::complex<double> (8.080151689298392e+10 , -9.365894433132759e+10), 
             std::complex<double> (2.235463725189983e+10 , -5.948444658785091e+10), 
             std::complex<double> (5.745317216740170e+09 , -3.150601300483326e+10), 
             std::complex<double> (1.449234738065675e+09 , -1.599523850540641e+10), 
             std::complex<double> (3.638488843256133e+08 , -8.036209551752358e+09), 
             std::complex<double> (9.124119725197355e+07 , -4.026975908568357e+09), 
             std::complex<double> (2.287349027494256e+07 , -2.016615694224953e+09), 
             std::complex<double> (5.733788895445901e+06 , -1.009708855536060e+09), 
             std::complex<double> (1.437284983727629e+06 , -5.055351633723626e+08), 
             std::complex<double> (3.602815783340446e+05 , -2.531058053673864e+08), 
             std::complex<double> (9.031101481455479e+04 , -1.267219141238533e+08), 
             std::complex<double> (2.263806297092661e+04 , -6.344553464017255e+07), 
             std::complex<double> (5.674632929256124e+03 , -3.176510741323538e+07), 
             std::complex<double> (1.422447651755881e+03 , -1.590375136294847e+07), 
             std::complex<double> (3.565618397870513e+02 , -7.962488574325078e+06), 
             std::complex<double> (8.937860016321510e+01 , -3.986557805603319e+06), 
             std::complex<double> (2.240435461492090e+01 , -1.995939206553575e+06), 
             std::complex<double> (5.616060491211885e+00 , -9.993015302610622e+05), 
             std::complex<double> (1.407771168802800e+00 , -5.003176174197889e+05), 
             std::complex<double> (3.528854391884821e-01 , -2.504926798531449e+05), 
             std::complex<double> (8.845817952176574e-02 , -1.254134982935338e+05), 
             std::complex<double> (2.217414904179711e-02 , -6.279043987205452e+04), 
             std::complex<double> (5.558597554302266e-03 , -3.143712113477009e+04), 
             std::complex<double> (1.393463586730668e-03 , -1.573953910379059e+04), 
             std::complex<double> (3.493274824860981e-04 , -7.880272770808787e+03), 
             std::complex<double> (8.757428434491133e-05 , -3.945394906678375e+03), 
             std::complex<double> (2.195897676440100e-05 , -1.975330028990690e+03), 
             std::complex<double> (5.510747522603683e-06 , -9.889826714095484e+02), 
             std::complex<double> (1.385290938447173e-06 , -4.951502691389830e+02), 
             std::complex<double> (3.492409888766843e-07 , -2.479034944790983e+02), 
             std::complex<double> (8.806072942349475e-08 , -1.241130577767191e+02), 
             std::complex<double> (2.209550055426092e-08 , -6.213112609130783e+01), 
             std::complex<double> (5.543483469188144e-09 , -3.109062229807737e+01), 
             std::complex<double> (1.413668111604913e-09 , -1.553354716842355e+01), 
             std::complex<double> (3.657274168660119e-10 , -7.715176851404666e+00), 
             std::complex<double> (9.523575611016325e-11 , -3.790763970499301e+00), 
             std::complex<double> (1.490208229386450e-11 , -2.198773295701914e+00), 
             std::complex<double> (1.487538366356832e-13 , -2.001307672244305e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 6.819410646545146e+10;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.233063093947589e+11 , 3.886593358749860e+10), 
             std::complex<double> (-6.295005529345710e+10 , 6.760619481505519e+10), 
             std::complex<double> (-2.202929806979887e+10 , 5.000930736702473e+10), 
             std::complex<double> (-6.684421486831000e+09 , 2.935025022025160e+10), 
             std::complex<double> (-1.939656294074093e+09 , 1.609879816721802e+10), 
             std::complex<double> (-5.555174167453873e+08 , 8.659997388731089e+09), 
             std::complex<double> (-1.585026682299974e+08 , 4.632600656020922e+09), 
             std::complex<double> (-4.517611163502118e+07 , 2.474244838119772e+09), 
             std::complex<double> (-1.287206096674998e+07 , 1.320882090064760e+09), 
             std::complex<double> (-3.667325188798876e+06 , 7.050654991503648e+08), 
             std::complex<double> (-1.044816336692254e+06 , 3.763387441777830e+08), 
             std::complex<double> (-2.976647517980348e+05 , 2.008740593667665e+08), 
             std::complex<double> (-8.480354131036608e+04 , 1.072179468827314e+08), 
             std::complex<double> (-2.416018814143966e+04 , 5.722828726793589e+07), 
             std::complex<double> (-6.883139240223544e+03 , 3.054596890449500e+07), 
             std::complex<double> (-1.960978209527619e+03 , 1.630410757626035e+07), 
             std::complex<double> (-5.586746450356742e+02 , 8.702422223743515e+06), 
             std::complex<double> (-1.591641011210570e+02 , 4.644973776248476e+06), 
             std::complex<double> (-4.534519000027903e+01 , 2.479284593753169e+06), 
             std::complex<double> (-1.291865277349106e+01 , 1.323334078974005e+06), 
             std::complex<double> (-3.680466919741806e+00 , 7.063380657425054e+05), 
             std::complex<double> (-1.048547944916449e+00 , 3.770124801100980e+05), 
             std::complex<double> (-2.987266178599685e-01 , 2.012328332308975e+05), 
             std::complex<double> (-8.510620527872556e-02 , 1.074093175638488e+05), 
             std::complex<double> (-2.424663762151651e-02 , 5.733041346228250e+04), 
             std::complex<double> (-6.907856556376299e-03 , 3.060047688687298e+04), 
             std::complex<double> (-1.968001888140749e-03 , 1.633320188958294e+04), 
             std::complex<double> (-5.606260664002580e-04 , 8.717952884880904e+03), 
             std::complex<double> (-1.596933206252242e-04 , 4.653266077629039e+03), 
             std::complex<double> (-4.548962854141386e-05 , 2.483715760318971e+03), 
             std::complex<double> (-1.296021053810168e-05 , 1.325708795424703e+03), 
             std::complex<double> (-3.694261449336069e-06 , 7.076234834857372e+02), 
             std::complex<double> (-1.053615766760628e-06 , 3.777321101946091e+02), 
             std::complex<double> (-3.006055470841920e-07 , 2.016797505912135e+02), 
             std::complex<double> (-8.570477437377466e-08 , 1.077654749773149e+02), 
             std::complex<double> (-2.421283676088920e-08 , 5.774033972822340e+01), 
             std::complex<double> (-6.409753640296327e-09 , 3.122497845132783e+01), 
             std::complex<double> (-1.245535805161496e-09 , 1.732451124674103e+01), 
             std::complex<double> (-7.404040085441005e-11 , 9.494141891652125e+00), 
             std::complex<double> (-2.050132095128929e-13 , 3.141735081105737e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.465118547022579e+10 , 2.093220829935340e+10), 
             std::complex<double> (1.830958305082091e+10 , -1.307643934371990e+09), 
             std::complex<double> (8.041663491442461e+09 , -7.356606943762166e+09), 
             std::complex<double> (2.603368892097106e+09 , -5.419038711361868e+09), 
             std::complex<double> (7.695417420646356e+08 , -3.147170042352588e+09), 
             std::complex<double> (2.215627186385823e+08 , -1.719871514943714e+09), 
             std::complex<double> (6.331238179280723e+07 , -9.241410089746137e+08), 
             std::complex<double> (1.805290078169994e+07 , -4.942037286669359e+08), 
             std::complex<double> (5.144453523510835e+06 , -2.639270909474180e+08), 
             std::complex<double> (1.465735762659553e+06 , -1.408944810820993e+08), 
             std::complex<double> (4.175904139532322e+05 , -7.520663452194430e+07), 
             std::complex<double> (1.189704769722787e+05 , -4.014252614512323e+07), 
             std::complex<double> (3.389425806841782e+04 , -2.142640751692654e+07), 
             std::complex<double> (9.656340261115569e+03 , -1.143649425888393e+07), 
             std::complex<double> (2.751052260642191e+03 , -6.104304042796352e+06), 
             std::complex<double> (7.837635623842881e+02 , -3.258211764940087e+06), 
             std::complex<double> (2.232910321480046e+02 , -1.739091500941681e+06), 
             std::complex<double> (6.361469922974091e+01 , -9.282512666781858e+05), 
             std::complex<double> (1.812356506508652e+01 , -4.954600775780139e+05), 
             std::complex<double> (5.163328796549208e+00 , -2.644549993620467e+05), 
             std::complex<double> (1.471010540297020e+00 , -1.411545547746814e+05), 
             std::complex<double> (4.190841937544170e-01 , -7.534215037820653e+04), 
             std::complex<double> (1.193950308609848e-01 , -4.021435674664161e+04), 
             std::complex<double> (3.401508311031978e-02 , -2.146467122491296e+04), 
             std::complex<double> (9.690795370584563e-03 , -1.145690611003498e+04), 
             std::complex<double> (2.760910937221205e-03 , -6.115197152070576e+03), 
             std::complex<double> (7.865951615107165e-04 , -3.264025608666328e+03), 
             std::complex<double> (2.240959842166102e-04 , -1.742194346241094e+03), 
             std::complex<double> (6.383497488331847e-05 , -9.299068801746106e+02), 
             std::complex<double> (1.818260672086020e-05 , -4.963427519773529e+02), 
             std::complex<double> (5.179019680093155e-06 , -2.649242232774399e+02), 
             std::complex<double> (1.475538253956464e-06 , -1.414014294834601e+02), 
             std::complex<double> (4.206723634968911e-07 , -7.546722142743826e+01), 
             std::complex<double> (1.200081716501330e-07 , -4.026856732529073e+01), 
             std::complex<double> (3.427478319649621e-08 , -2.147014018045991e+01), 
             std::complex<double> (9.860467176147622e-09 , -1.141623593098393e+01), 
             std::complex<double> (2.879767906246125e-09 , -6.017089153129619e+00), 
             std::complex<double> (7.941984902259705e-10 , -3.158204500259297e+00), 
             std::complex<double> (9.333783894965729e-11 , -2.095194761873083e+00), 
             std::complex<double> (5.467232781849662e-13 , -2.000392276848583e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 8.767608340531406e+09;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.610661672721688e+10 , 4.643816589270030e+09), 
             std::complex<double> (-9.053164371726076e+09 , 8.716494696091959e+09), 
             std::complex<double> (-3.591313491142742e+09 , 7.053705479160174e+09), 
             std::complex<double> (-1.229272700672112e+09 , 4.464770676088959e+09), 
             std::complex<double> (-3.994211268415112e+08 , 2.609323249948285e+09), 
             std::complex<double> (-1.275801907840611e+08 , 1.486413376230754e+09), 
             std::complex<double> (-4.052780857353461e+07 , 8.398724705676417e+08), 
             std::complex<double> (-1.285183950712821e+07 , 4.733304543902467e+08), 
             std::complex<double> (-4.073212774435883e+06 , 2.665383302636197e+08), 
             std::complex<double> (-1.290721991356240e+06 , 1.500521143525332e+08), 
             std::complex<double> (-4.089819839943156e+05 , 8.446734237943131e+07), 
             std::complex<double> (-1.295889574142116e+05 , 4.754712014277605e+07), 
             std::complex<double> (-4.106098815384283e+04 , 2.676430788050086e+07), 
             std::complex<double> (-1.301038180154414e+04 , 1.506560858066691e+07), 
             std::complex<double> (-4.122403062976586e+03 , 8.480412595076224e+06), 
             std::complex<double> (-1.306203329104089e+03 , 4.773612595763293e+06), 
             std::complex<double> (-4.138768101624122e+02 , 2.687059746438675e+06), 
             std::complex<double> (-1.311388484187883e+02 , 1.512542072844161e+06), 
             std::complex<double> (-4.155196172073543e+01 , 8.514077538031481e+05), 
             std::complex<double> (-1.316593020099958e+01 , 4.792561968152644e+05), 
             std::complex<double> (-4.171685799708315e+00 , 2.697726218716812e+05), 
             std::complex<double> (-1.321819325120319e+00 , 1.518546197359088e+05), 
             std::complex<double> (-4.188259805506248e-01 , 8.547874652668991e+04), 
             std::complex<double> (-1.327078198553247e-01 , 4.811586339718426e+04), 
             std::complex<double> (-4.204951108709425e-02 , 2.708435063051265e+04), 
             std::complex<double> (-1.332373286169278e-02 , 1.524574286399132e+04), 
             std::complex<double> (-4.221693620706279e-03 , 8.581808388287651e+03), 
             std::complex<double> (-1.337622638041483e-03 , 4.830690628082541e+03), 
             std::complex<double> (-4.237972052053582e-04 , 2.719194235821104e+03), 
             std::complex<double> (-1.342577676804125e-04 , 1.530640196045328e+03), 
             std::complex<double> (-4.251824176118857e-05 , 8.616123617127558e+02), 
             std::complex<double> (-1.345292221624769e-05 , 4.850309180194964e+02), 
             std::complex<double> (-4.247215505066930e-06 , 2.730774987242593e+02), 
             std::complex<double> (-1.333944594564809e-06 , 1.538113817636532e+02), 
             std::complex<double> (-4.140617757977241e-07 , 8.675151323219090e+01), 
             std::complex<double> (-1.245937376871465e-07 , 4.913607539941756e+01), 
             std::complex<double> (-3.359434021286551e-08 , 2.818825776567928e+01), 
             std::complex<double> (-5.849129029473393e-09 , 1.663323083835730e+01), 
             std::complex<double> (-2.410291864794301e-10 , 9.453485964597251e+00), 
             std::complex<double> (-3.727341490860978e-13 , 3.141632582666789e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.616973733100822e+09 , 2.571055649401248e+09), 
             std::complex<double> (2.284977120198832e+09 , 8.661491429684791e+07), 
             std::complex<double> (1.168363898460488e+09 , -8.499632885645322e+08), 
             std::complex<double> (4.332690546151149e+08 , -7.271804345050156e+08), 
             std::complex<double> (1.444070136134930e+08 , -4.606358184700190e+08), 
             std::complex<double> (4.649915091918428e+07 , -2.688810025149615e+08), 
             std::complex<double> (1.480900101320543e+07 , -1.530888462639134e+08), 
             std::complex<double> (4.699914577101156e+06 , -8.648460985240068e+07), 
             std::complex<double> (1.489955560392514e+06 , -4.873763581290583e+07), 
             std::complex<double> (4.721763958162231e+05 , -2.744426156481187e+07), 
             std::complex<double> (1.496190661809330e+05 , -1.545010421633239e+07), 
             std::complex<double> (4.740829143913260e+04 , -8.697156896988722e+06), 
             std::complex<double> (1.502162169454108e+04 , -4.895673346130614e+06), 
             std::complex<double> (4.759680604375374e+03 , -2.755777500460585e+06), 
             std::complex<double> (1.508128441344631e+03 , -1.551224969096821e+06), 
             std::complex<double> (4.778578243255149e+02 , -8.731826162070866e+05), 
             std::complex<double> (1.514115571567000e+02 , -4.915132917445675e+05), 
             std::complex<double> (4.797548258273670e+01 , -2.766721329085908e+05), 
             std::complex<double> (1.520126210499702e+01 , -1.557383462331809e+05), 
             std::complex<double> (4.816590909208264e+00 , -8.766489071013530e+04), 
             std::complex<double> (1.526158426023936e+00 , -4.934644054491052e+04), 
             std::complex<double> (4.835701513536483e-01 , -2.777704017799139e+04), 
             std::complex<double> (1.532215720993542e-01 , -1.563565580936367e+04), 
             std::complex<double> (4.854913717611598e-02 , -8.801288039593481e+03), 
             std::complex<double> (1.538313018094348e-02 , -4.954232238658193e+03), 
             std::complex<double> (4.874280665100051e-03 , -2.788730001182036e+03), 
             std::complex<double> (1.544466933006866e-03 , -1.569771778059184e+03), 
             std::complex<double> (4.893754466976076e-04 , -8.836217076683749e+02), 
             std::complex<double> (1.550562516061321e-04 , -4.973883872116578e+02), 
             std::complex<double> (4.912876544759819e-05 , -2.799774341346834e+02), 
             std::complex<double> (1.556607337791191e-05 , -1.575957468364490e+02), 
             std::complex<double> (4.931058801125233e-06 , -8.870482806430422e+01), 
             std::complex<double> (1.561223857528711e-06 , -4.992188851476876e+01), 
             std::complex<double> (4.935854222620059e-07 , -2.808331872815572e+01), 
             std::complex<double> (1.555119962643026e-07 , -1.577675284693110e+01), 
             std::complex<double> (4.888966480993283e-08 , -8.825588831853372e+00), 
             std::complex<double> (1.533442783581276e-08 , -4.879445917161496e+00), 
             std::complex<double> (4.147180446793749e-09 , -2.733110063583577e+00), 
             std::complex<double> (3.420152555682672e-10 , -2.042980526121471e+00), 
             std::complex<double> (1.082387178464305e-12 , -2.000114689334851e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 1.563082652204966e+09;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.907347511144252e+09 , 7.724039573343986e+08), 
             std::complex<double> (-1.765747824695564e+09 , 1.542115710523208e+09), 
             std::complex<double> (-7.782069398337824e+08 , 1.347815058517721e+09), 
             std::complex<double> (-2.954009947906344e+08 , 9.122444485838728e+08), 
             std::complex<double> (-1.057749158935667e+08 , 5.639519953938454e+08), 
             std::complex<double> (-3.708246196328126e+07 , 3.377070823419682e+08), 
             std::complex<double> (-1.290387685867627e+07 , 1.999937736846150e+08), 
             std::complex<double> (-4.478624662514562e+06 , 1.179826476320135e+08), 
             std::complex<double> (-1.553022292221680e+06 , 6.950866578235863e+07), 
             std::complex<double> (-5.383626840240620e+05 , 4.093156927606445e+07), 
             std::complex<double> (-1.866058016054107e+05 , 2.409950029619724e+07), 
             std::complex<double> (-6.467836273928163e+04 , 1.418840204411208e+07), 
             std::complex<double> (-2.241750293020049e+04 , 8.353155175602375e+06), 
             std::complex<double> (-7.769865233301661e+03 , 4.917730285754640e+06), 
             std::complex<double> (-2.693016676141433e+03 , 2.895195235813679e+06), 
             std::complex<double> (-9.333925532483522e+02 , 1.704475078466238e+06), 
             std::complex<double> (-3.235113459850519e+02 , 1.003467554684307e+06), 
             std::complex<double> (-1.121281675092420e+02 , 5.907666523194790e+05), 
             std::complex<double> (-3.886332883033633e+01 , 3.477992130311116e+05), 
             std::complex<double> (-1.346992282741687e+01 , 2.047581574452301e+05), 
             std::complex<double> (-4.668632028598083e+00 , 1.205462845646340e+05), 
             std::complex<double> (-1.618131309090977e+00 , 7.096863400428811e+04), 
             std::complex<double> (-5.608396712711295e-01 , 4.178102262031913e+04), 
             std::complex<double> (-1.943862624114411e-01 , 2.459754097613390e+04), 
             std::complex<double> (-6.737411332605869e-02 , 1.448119304279567e+04), 
             std::complex<double> (-2.335181629214159e-02 , 8.525445051947156e+03), 
             std::complex<double> (-8.093880894964935e-03 , 5.019147606715412e+03), 
             std::complex<double> (-2.805530386959645e-03 , 2.954904361673518e+03), 
             std::complex<double> (-9.725074005272411e-04 , 1.739636397084562e+03), 
             std::complex<double> (-3.371090402953826e-04 , 1.024184371617434e+03), 
             std::complex<double> (-1.168360683456084e-04 , 6.029913504357919e+02), 
             std::complex<double> (-4.046239368609373e-05 , 3.550441244007910e+02), 
             std::complex<double> (-1.397941032449423e-05 , 2.091048387833814e+02), 
             std::complex<double> (-4.800772076943351e-06 , 1.232436196415769e+02), 
             std::complex<double> (-1.624406049955491e-06 , 7.279170935439704e+01), 
             std::complex<double> (-5.262502423559682e-07 , 4.325330555897848e+01), 
             std::complex<double> (-1.465199993035186e-07 , 2.612498896184864e+01), 
             std::complex<double> (-2.303827143522598e-08 , 1.621892036988437e+01), 
             std::complex<double> (-6.970580647162130e-10 , 9.436111337656246e+00), 
             std::complex<double> (-8.193513926457782e-13 , 3.141603672944739e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (2.497337894786215e+08 , 4.368282941921436e+08), 
             std::complex<double> (3.898206589092215e+08 , 5.295061029987082e+07), 
             std::complex<double> (2.268112309483741e+08 , -1.309341943075704e+08), 
             std::complex<double> (9.470866342345458e+07 , -1.309034027352926e+08), 
             std::complex<double> (3.505384241105291e+07 , -9.015958773598580e+07), 
             std::complex<double> (1.243110475743731e+07 , -5.592204097099373e+07), 
             std::complex<double> (4.343012809978185e+06 , -3.351549463081950e+07), 
             std::complex<double> (1.509437704910050e+06 , -1.985328670118356e+07), 
             std::complex<double> (5.236680280354188e+05 , -1.171306272258705e+07), 
             std::complex<double> (1.815621681942994e+05 , -6.900867186258990e+06), 
             std::complex<double> (6.293620093860045e+04 , -4.063753697811282e+06), 
             std::complex<double> (2.181438914464249e+04 , -2.392646258744826e+06), 
             std::complex<double> (7.560913645580136e+03 , -1.408654394688695e+06), 
             std::complex<double> (2.620605071047464e+03 , -8.293191514502283e+05), 
             std::complex<double> (9.082962391113452e+02 , -4.882428732753578e+05), 
             std::complex<double> (3.148132225519375e+02 , -2.874412437632169e+05), 
             std::complex<double> (1.091134247629300e+02 , -1.692239745221116e+05), 
             std::complex<double> (3.781841513749194e+01 , -9.962643106448060e+04), 
             std::complex<double> (1.310776335399662e+01 , -5.865259220404935e+04), 
             std::complex<double> (4.543118410465391e+00 , -3.453025883110186e+04), 
             std::complex<double> (1.574632584801455e+00 , -2.032883312959338e+04), 
             std::complex<double> (5.457619972253657e-01 , -1.196809603519004e+04), 
             std::complex<double> (1.891588768114733e-01 , -7.045919492437385e+03), 
             std::complex<double> (6.556196983586285e-02 , -4.148110143203874e+03), 
             std::complex<double> (2.272374288353910e-02 , -2.442096714678605e+03), 
             std::complex<double> (7.875992491065385e-03 , -1.437723504272329e+03), 
             std::complex<double> (2.729770835270368e-03 , -8.464234647792599e+02), 
             std::complex<double> (9.461524819415148e-04 , -4.983098401341786e+02), 
             std::complex<double> (3.279649621861561e-04 , -2.933659036596868e+02), 
             std::complex<double> (1.136902389785451e-04 , -1.727090945383379e+02), 
             std::complex<double> (3.941614129449380e-05 , -1.016734350602007e+02), 
             std::complex<double> (1.366775124773620e-05 , -5.984962817116728e+01), 
             std::complex<double> (4.738818693738254e-06 , -3.522124875771384e+01), 
             std::complex<double> (1.641113023885773e-06 , -2.071229125324714e+01), 
             std::complex<double> (5.674218636329427e-07 , -1.215416621367042e+01), 
             std::complex<double> (1.962545066726201e-07 , -7.088712693028019e+00), 
             std::complex<double> (6.746632821637535e-08 , -4.077236899881440e+00), 
             std::complex<double> (1.787340333528884e-08 , -2.448381223000574e+00), 
             std::complex<double> (1.069635178400833e-09 , -2.018319335566692e+00), 
             std::complex<double> (2.402486204385628e-12 , -2.000032837342133e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 3.597493655870090e+08;
             std::complex<double> zvec1[] = {
             std::complex<double> (-6.757929558574898e+08 , 1.664713583363300e+08), 
             std::complex<double> (-4.372002913964674e+08 , 3.496303185812742e+08), 
             std::complex<double> (-2.106453818015416e+08 , 3.265260638143334e+08), 
             std::complex<double> (-8.745943473669402e+07 , 2.346099468608779e+08), 
             std::complex<double> (-3.407926908881292e+07 , 1.525340805277402e+08), 
             std::complex<double> (-1.295083008629596e+07 , 9.547481054111528e+07), 
             std::complex<double> (-4.874737200563508e+06 , 5.891054016331727e+07), 
             std::complex<double> (-1.828261639255609e+06 , 3.615462071028834e+07), 
             std::complex<double> (-6.847587043800608e+05 , 2.214422265356991e+07), 
             std::complex<double> (-2.563400778194122e+05 , 1.355282368658760e+07), 
             std::complex<double> (-9.594293369958854e+04 , 8.292331171473768e+06), 
             std::complex<double> (-3.590695713795506e+04 , 5.073149276407419e+06), 
             std::complex<double> (-1.343793858490837e+04 , 3.103569584341675e+06), 
             std::complex<double> (-5.029008526599193e+03 , 1.898623759716890e+06), 
             std::complex<double> (-1.882047047445277e+03 , 1.161485875297074e+06), 
             std::complex<double> (-7.043328998930591e+02 , 7.105392198459121e+05), 
             std::complex<double> (-2.635877568499291e+02 , 4.346721750592028e+05), 
             std::complex<double> (-9.864439477907760e+01 , 2.659105081487501e+05), 
             std::complex<double> (-3.691642612809024e+01 , 1.626706165169031e+05), 
             std::complex<double> (-1.381550956016815e+01 , 9.951366200834389e+04), 
             std::complex<double> (-5.170274872813414e+00 , 6.087742838242603e+04), 
             std::complex<double> (-1.934900695949033e+00 , 3.724173367350674e+04), 
             std::complex<double> (-7.241028989448756e-01 , 2.278261078161663e+04), 
             std::complex<double> (-2.709794933849554e-01 , 1.393725088986471e+04), 
             std::complex<double> (-1.014067814475907e-01 , 8.526107754745986e+03), 
             std::complex<double> (-3.794857907721268e-02 , 5.215845374286750e+03), 
             std::complex<double> (-1.420128602616471e-02 , 3.190796539536440e+03), 
             std::complex<double> (-5.314641373456001e-03 , 1.951977863633049e+03), 
             std::complex<double> (-1.989073462021171e-03 , 1.194137478514606e+03), 
             std::complex<double> (-7.444727058733413e-04 , 7.305392641779973e+02), 
             std::complex<double> (-2.785935492041724e-04 , 4.469500227108133e+02), 
             std::complex<double> (-1.041697874742726e-04 , 2.734917816578494e+02), 
             std::complex<double> (-3.885029775980553e-05 , 1.674234863546034e+02), 
             std::complex<double> (-1.438085862810561e-05 , 1.026094281347010e+02), 
             std::complex<double> (-5.210915924018490e-06 , 6.307942476532812e+01), 
             std::complex<double> (-1.771743570475432e-06 , 3.909277798924850e+01), 
             std::complex<double> (-4.871556936491854e-07 , 2.470837551989106e+01), 
             std::complex<double> (-6.451344870642707e-08 , 1.597850211276049e+01), 
             std::complex<double> (-1.364027307617252e-09 , 9.429065404812917e+00), 
             std::complex<double> (-1.106422446541182e-12 , 3.141595654835673e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (5.021229629899424e+07 , 9.573426248688118e+07), 
             std::complex<double> (8.503658702053596e+07 , 1.916561491759082e+07), 
             std::complex<double> (5.525644276706799e+07 , -2.500386205323650e+07), 
             std::complex<double> (2.561522955845227e+07 , -2.958192196709283e+07), 
             std::complex<double> (1.040165584109257e+07 , -2.211621171347208e+07), 
             std::complex<double> (4.014428794357887e+06 , -1.452511290783674e+07), 
             std::complex<double> (1.519817102089311e+06 , -9.120511098008595e+06), 
             std::complex<double> (5.712411633161930e+05 , -5.633820421544250e+06), 
             std::complex<double> (2.141267423762998e+05 , -3.458982067629103e+06), 
             std::complex<double> (8.018288424430869e+04 , -2.118895560447226e+06), 
             std::complex<double> (3.001424968198962e+04 , -1.296889525250297e+06), 
             std::complex<double> (1.123340896130684e+04 , -7.935217430994142e+05), 
             std::complex<double> (4.204095139347477e+03 , -4.854708974808910e+05), 
             std::complex<double> (1.573348184128247e+03 , -2.969944305593917e+05), 
             std::complex<double> (5.888082824190428e+02 , -1.816879818254871e+05), 
             std::complex<double> (2.203544256809051e+02 , -1.111479350313492e+05), 
             std::complex<double> (8.246490999967402e+01 , -6.799478212729240e+04), 
             std::complex<double> (3.086145346157888e+01 , -4.159579187759199e+04), 
             std::complex<double> (1.154950870899146e+01 , -2.544620754551074e+04), 
             std::complex<double> (4.322260374810132e+00 , -1.556670445547073e+04), 
             std::complex<double> (1.617553585390565e+00 , -9.522923024249372e+03), 
             std::complex<double> (6.053495386271041e-01 , -5.825642854377887e+03), 
             std::complex<double> (2.265439252257474e-01 , -3.563833664786646e+03), 
             std::complex<double> (8.478026332952812e-02 , -2.180173054066974e+03), 
             std::complex<double> (3.172705425842923e-02 , -1.333719304104711e+03), 
             std::complex<double> (1.187292562038440e-02 , -8.159014677217300e+02), 
             std::complex<double> (4.443061545970862e-03 , -4.991263273560290e+02), 
             std::complex<double> (1.662668710382170e-03 , -3.053387126875481e+02), 
             std::complex<double> (6.222243967626335e-04 , -1.867882689298264e+02), 
             std::complex<double> (2.328909259342514e-04 , -1.142635017832522e+02), 
             std::complex<double> (8.718629119580242e-05 , -6.989390121133002e+01), 
             std::complex<double> (3.264631534860909e-05 , -4.274654572583646e+01), 
             std::complex<double> (1.222793440898706e-05 , -2.613216178496151e+01), 
             std::complex<double> (4.581577656949828e-06 , -1.595682388476035e+01), 
             std::complex<double> (1.717608884935356e-06 , -9.713170677780306e+00), 
             std::complex<double> (6.443895273476331e-07 , -5.863918036472604e+00), 
             std::complex<double> (2.355597460701598e-07 , -3.490959616758682e+00), 
             std::complex<double> (5.581060396181204e-08 , -2.262430565379247e+00), 
             std::complex<double> (2.262389730855300e-09 , -2.007407589661226e+00), 
             std::complex<double> (3.375789190314810e-12 , -2.000009240487879e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 1.014613602735905e+08;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.921298399391937e+08 , 4.411810006779376e+07), 
             std::complex<double> (-1.309494937489181e+08 , 9.663406920712876e+07), 
             std::complex<double> (-6.807544958308402e+07 , 9.557567298880142e+07), 
             std::complex<double> (-3.056719370227870e+07 , 7.244434603880116e+07), 
             std::complex<double> (-1.282991474460506e+07 , 4.929886967980125e+07), 
             std::complex<double> (-5.232896374310734e+06 , 3.210958024437876e+07), 
             std::complex<double> (-2.109381616336855e+06 , 2.054738737514637e+07), 
             std::complex<double> (-8.462635041149040e+05 , 1.305563208140372e+07), 
             std::complex<double> (-3.388658761826573e+05 , 8.271901465431931e+06), 
             std::complex<double> (-1.355870790070452e+05 , 5.235028754150138e+06), 
             std::complex<double> (-5.423454780303057e+04 , 3.311580310229713e+06), 
             std::complex<double> (-2.169105448313469e+04 , 2.094462173327443e+06), 
             std::complex<double> (-8.674892818394434e+03 , 1.324579897681296e+06), 
             std::complex<double> (-3.469278272680851e+03 , 8.376665324502613e+05), 
             std::complex<double> (-1.387429096730825e+03 , 5.297355469988022e+05), 
             std::complex<double> (-5.548569922154544e+02 , 3.350001862007834e+05), 
             std::complex<double> (-2.218965616189625e+02 , 2.118508313338263e+05), 
             std::complex<double> (-8.874000063767203e+01 , 1.339722878806102e+05), 
             std::complex<double> (-3.548850276807578e+01 , 8.472267256430519e+04), 
             std::complex<double> (-1.419239245410664e+01 , 5.357772519806287e+04), 
             std::complex<double> (-5.675754850035216e+00 , 3.388198749163363e+04), 
             std::complex<double> (-2.269822844993397e+00 , 2.142661118979817e+04), 
             std::complex<double> (-9.077373037127995e-01 , 1.354996355298606e+04), 
             std::complex<double> (-3.630169283001943e-01 , 8.568855909327011e+03), 
             std::complex<double> (-1.451744119223127e-01 , 5.418857655487433e+03), 
             std::complex<double> (-5.805597115309982e-02 , 3.426834435421355e+03), 
             std::complex<double> (-2.321618407667637e-02 , 2.167103697694515e+03), 
             std::complex<double> (-9.283303310394741e-03 , 1.370469083810355e+03), 
             std::complex<double> (-3.711453141871698e-03 , 8.666949000740808e+02), 
             std::complex<double> (-1.483330493985827e-03 , 5.481278448477266e+02), 
             std::complex<double> (-5.923672491767484e-04 , 3.466921929181578e+02), 
             std::complex<double> (-2.361181712478414e-04 , 2.193424625752881e+02), 
             std::complex<double> (-9.369480092950490e-05 , 1.388649185253084e+02), 
             std::complex<double> (-3.678071268751481e-05 , 8.806232953845476e+01), 
             std::complex<double> (-1.403861979031402e-05 , 5.607940283582973e+01), 
             std::complex<double> (-4.926490814993237e-06 , 3.608119254060830e+01), 
             std::complex<double> (-1.300310160202216e-06 , 2.373799429271499e+01), 
             std::complex<double> (-1.387999663190090e-07 , 1.584477852379612e+01), 
             std::complex<double> (-2.025137578642771e-09 , 9.426340415696112e+00), 
             std::complex<double> (-1.159772172168061e-12 , 3.141593461908113e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.246639116695216e+07 , 2.571362615571233e+07), 
             std::complex<double> (2.258771346510689e+07 , 6.970091021206074e+06), 
             std::complex<double> (1.614048604932548e+07 , -5.582181370644691e+06), 
             std::complex<double> (8.207896193435946e+06 , -7.994748690941910e+06), 
             std::complex<double> (3.620733203720899e+06 , -6.485179102127224e+06), 
             std::complex<double> (1.506466749955960e+06 , -4.499161663663012e+06), 
             std::complex<double> (6.121131110691743e+05 , -2.949863042766987e+06), 
             std::complex<double> (2.463575349610918e+05 , -1.892341784543781e+06), 
             std::complex<double> (9.877374058162920e+04 , -1.203537403010204e+06), 
             std::complex<double> (3.954151368307140e+04 , -7.628386133020810e+05), 
             std::complex<double> (1.581974580395497e+04 , -4.828501078634487e+05), 
             std::complex<double> (6.327606980469617e+03 , -3.054603742638551e+05), 
             std::complex<double> (2.530678876673760e+03 , -1.931980016907821e+05), 
             std::complex<double> (1.012086559791048e+03 , -1.221834755304546e+05), 
             std::complex<double> (4.047543897120896e+02 , -7.726933251349335e+04), 
             std::complex<double> (1.618686916350659e+02 , -4.886476169363668e+04), 
             std::complex<double> (6.473409570981391e+01 , -3.090167250891183e+04), 
             std::complex<double> (2.588824778407673e+01 , -1.954191934360788e+04), 
             std::complex<double> (1.035312496177912e+01 , -1.235811007077327e+04), 
             std::complex<double> (4.140373386908720e+00 , -7.815139687146089e+03), 
             std::complex<double> (1.655797620187784e+00 , -4.942211922368764e+03), 
             std::complex<double> (6.621788179877617e-01 , -3.125402433678355e+03), 
             std::complex<double> (2.648160280637075e-01 , -1.976471185508419e+03), 
             std::complex<double> (1.059043000700572e-01 , -1.249899078839237e+03), 
             std::complex<double> (4.235278368110191e-02 , -7.904223581419200e+02), 
             std::complex<double> (1.693744796849116e-02 , -4.998538111026916e+02), 
             std::complex<double> (6.773502714409515e-03 , -3.161008094278062e+02), 
             std::complex<double> (2.708801934770430e-03 , -1.998965180502970e+02), 
             std::complex<double> (1.083242235716243e-03 , -1.264088268013140e+02), 
             std::complex<double> (4.331659849709797e-04 , -7.993388798257237e+01), 
             std::complex<double> (1.732057476218584e-04 , -5.054030369263095e+01), 
             std::complex<double> (6.925184332833776e-05 , -3.194684470569799e+01), 
             std::complex<double> (2.768038637046936e-05 , -2.018017324127766e+01), 
             std::complex<double> (1.105830188665991e-05 , -1.272570840892989e+01), 
             std::complex<double> (4.422861272282279e-06 , -7.990248121945951e+00), 
             std::complex<double> (1.776432828792962e-06 , -4.964258008290898e+00), 
             std::complex<double> (6.753436106737522e-07 , -3.053743004207058e+00), 
             std::complex<double> (1.341813225842820e-07 , -2.145995625494450e+00), 
             std::complex<double> (3.604001177576182e-09 , -2.002860125373691e+00), 
             std::complex<double> (3.645919624456548e-12 , -2.000002562708367e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 3.372849818068297e+07;
             std::complex<double> zvec1[] = {
             std::complex<double> (-6.429063117899310e+07 , 1.382536283779675e+07), 
             std::complex<double> (-4.576585974019063e+07 , 3.136665907857788e+07), 
             std::complex<double> (-2.539679001577676e+07 , 3.260805476536937e+07), 
             std::complex<double> (-1.221618159512390e+07 , 2.593203722131532e+07), 
             std::complex<double> (-5.476874336357435e+06 , 1.839561696629719e+07), 
             std::complex<double> (-2.378268155268686e+06 , 1.242235812413723e+07), 
             std::complex<double> (-1.018431273742565e+06 , 8.213785221815328e+06), 
             std::complex<double> (-4.335131104676636e+05 , 5.382543119647412e+06), 
             std::complex<double> (-1.840622034209667e+05 , 3.513805748568872e+06), 
             std::complex<double> (-7.806498146729966e+04 , 2.290163401848242e+06), 
             std::complex<double> (-3.309391739017106e+04 , 1.491618615944708e+06), 
             std::complex<double> (-1.402669570518078e+04 , 9.712321916524942e+05), 
             std::complex<double> (-5.944654458337495e+03 , 6.323170914697213e+05), 
             std::complex<double> (-2.519315610627144e+03 , 4.116462213842663e+05), 
             std::complex<double> (-1.067657749108736e+03 , 2.679808442277437e+05), 
             std::complex<double> (-4.524586361738204e+02 , 1.744533454528707e+05), 
             std::complex<double> (-1.917453332052945e+02 , 1.135672552748869e+05), 
             std::complex<double> (-8.125879601521375e+01 , 7.393093982221054e+04), 
             std::complex<double> (-3.443623343777325e+01 , 4.812813325335967e+04), 
             std::complex<double> (-1.459356079119954e+01 , 3.133081424791512e+04), 
             std::complex<double> (-6.184561271340652e+00 , 2.039596719751718e+04), 
             std::complex<double> (-2.620954485901897e+00 , 1.327751917689954e+04), 
             std::complex<double> (-1.110737654189575e+00 , 8.643500052475230e+03), 
             std::complex<double> (-4.707177620204561e-01 , 5.626813540685920e+03), 
             std::complex<double> (-1.994807274615887e-01 , 3.662990984521928e+03), 
             std::complex<double> (-8.453348245830929e-02 , 2.384570786916154e+03), 
             std::complex<double> (-3.582158522626007e-02 , 1.552340727575078e+03), 
             std::complex<double> (-1.517896501240274e-02 , 1.010577670156192e+03), 
             std::complex<double> (-6.430981294561975e-03 , 6.579093549166789e+02), 
             std::complex<double> (-2.723490230949166e-03 , 4.283460804032580e+02), 
             std::complex<double> (-1.152076757011837e-03 , 2.789330163420434e+02), 
             std::complex<double> (-4.860211591947028e-04 , 1.817127186239851e+02), 
             std::complex<double> (-2.037330620376564e-04 , 1.184939572294430e+02), 
             std::complex<double> (-8.407337628012382e-05 , 7.744816578929685e+01), 
             std::complex<double> (-3.331415810425147e-05 , 5.089674660136445e+01), 
             std::complex<double> (-1.176643899885547e-05 , 3.386996437350400e+01), 
             std::complex<double> (-2.864688946140982e-06 , 2.308235712690259e+01), 
             std::complex<double> (-2.388707580183042e-07 , 1.577395525089171e+01), 
             std::complex<double> (-2.496845245240810e-09 , 9.425329236788178e+00), 
             std::complex<double> (-1.085884088875372e-12 , 3.141592869202212e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (3.673504009246490e+06 , 8.146271161489838e+06), 
             std::complex<double> (7.043023862767386e+06 , 2.724552630469848e+06), 
             std::complex<double> (5.463745970945064e+06 , -1.379851799255619e+06), 
             std::complex<double> (3.016952272996435e+06 , -2.491514822179004e+06), 
             std::complex<double> (1.433817366155757e+06 , -2.194495400715092e+06), 
             std::complex<double> (6.382075633756730e+05 , -1.605675377863806e+06), 
             std::complex<double> (2.761766029068796e+05 , -1.096579009107399e+06), 
             std::complex<double> (1.180831408468410e+05 , -7.283108575517322e+05), 
             std::complex<double> (5.023061834893578e+04 , -4.781435209079718e+05), 
             std::complex<double> (2.132097134070465e+04 , -3.123792417990398e+05), 
             std::complex<double> (9.041613355595406e+03 , -2.036627156244318e+05), 
             std::complex<double> (3.832794113550935e+03 , -1.326668206383934e+05), 
             std::complex<double> (1.624475676918181e+03 , -8.638787111657622e+04), 
             std::complex<double> (6.884627251406737e+02 , -5.624388207442427e+04), 
             std::complex<double> (2.917659804910053e+02 , -3.661584160057628e+04), 
             std::complex<double> (1.236469138146348e+02 , -2.383694279212383e+04), 
             std::complex<double> (5.239981247836497e+01 , -1.551768449569411e+04), 
             std::complex<double> (2.220627437836546e+01 , -1.010185399194996e+04), 
             std::complex<double> (9.410683548187754e+00 , -6.576189564935177e+03), 
             std::complex<double> (3.988095061160299e+00 , -4.281019054152628e+03), 
             std::complex<double> (1.690090259784373e+00 , -2.786889937655153e+03), 
             std::complex<double> (7.162379626257696e-01 , -1.814230104398216e+03), 
             std::complex<double> (3.035356743872196e-01 , -1.181040552752193e+03), 
             std::complex<double> (1.286369037874478e-01 , -7.688419116554246e+02), 
             std::complex<double> (5.451539715526936e-02 , -5.005054858092889e+02), 
             std::complex<double> (2.310264818588378e-02 , -3.258213871250342e+02), 
             std::complex<double> (9.790199154645674e-03 , -2.121035172713848e+02), 
             std::complex<double> (4.148831560431325e-03 , -1.380735011578393e+02), 
             std::complex<double> (1.758261007413990e-03 , -8.987917550844391e+01), 
             std::complex<double> (7.452078377491277e-04 , -5.850263379504320e+01), 
             std::complex<double> (3.158664791863635e-04 , -3.807283265103564e+01), 
             std::complex<double> (1.338757334037931e-04 , -2.476701737085556e+01), 
             std::complex<double> (5.673970920817093e-05 , -1.609538732124914e+01), 
             std::complex<double> (2.405950423860381e-05 , -1.043511266833778e+01), 
             std::complex<double> (1.021201855768473e-05 , -6.726705277328384e+00), 
             std::complex<double> (4.327005335681858e-06 , -4.281318911346707e+00), 
             std::complex<double> (1.630788644829155e-06 , -2.726418689856930e+00), 
             std::complex<double> (2.556501286869910e-07 , -2.076920187151913e+00), 
             std::complex<double> (4.688726700939065e-09 , -2.001061330312164e+00), 
             std::complex<double> (3.501036610277353e-12 , -2.000000701955252e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 1.283029586327183e+07;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.459007912784178e+07 , 4.972291833663177e+06), 
             std::complex<double> (-1.815748505686428e+07 , 1.162135362909335e+07), 
             std::complex<double> (-1.066073151038934e+07 , 1.261806928501339e+07), 
             std::complex<double> (-5.449447679683310e+06 , 1.047782943559382e+07), 
             std::complex<double> (-2.591250137642941e+06 , 7.720737387186867e+06), 
             std::complex<double> (-1.190014195496150e+06 , 5.389144786123645e+06), 
             std::complex<double> (-5.377897566460994e+05 , 3.670937484082993e+06), 
             std::complex<double> (-2.412729271159523e+05 , 2.473315429862499e+06), 
             std::complex<double> (-1.078904988188915e+05 , 1.658273558144286e+06), 
             std::complex<double> (-4.817500819771435e+04 , 1.109388473927277e+06), 
             std::complex<double> (-2.149692561599358e+04 , 7.414597361508097e+05), 
             std::complex<double> (-9.589682153767131e+03 , 4.953390317278366e+05), 
             std::complex<double> (-4.277358223565303e+03 , 3.308516096796244e+05), 
             std::complex<double> (-1.907751975625013e+03 , 2.209664625716690e+05), 
             std::complex<double> (-8.508578991268816e+02 , 1.475715739028376e+05), 
             std::complex<double> (-3.794785178178740e+02 , 9.855339265817212e+04), 
             std::complex<double> (-1.692447345084124e+02 , 6.581685364627738e+04), 
             std::complex<double> (-7.548178769318051e+01 , 4.395428011011957e+04), 
             std::complex<double> (-3.366421561355417e+01 , 2.935381611835436e+04), 
             std::complex<double> (-1.501391087611332e+01 , 1.960323231200959e+04), 
             std::complex<double> (-6.696036319041447e+00 , 1.309153910632633e+04), 
             std::complex<double> (-2.986338900259917e+00 , 8.742864564537293e+03), 
             std::complex<double> (-1.331849398916066e+00 , 5.838710649438138e+03), 
             std::complex<double> (-5.939705045141928e-01 , 3.899245316699059e+03), 
             std::complex<double> (-2.648947021705502e-01 , 2.604024617952774e+03), 
             std::complex<double> (-1.181392386147316e-01 , 1.739048440096345e+03), 
             std::complex<double> (-5.269148176379781e-02 , 1.161402903204316e+03), 
             std::complex<double> (-2.350180750512550e-02 , 7.756477359771625e+02), 
             std::complex<double> (-1.048130427736623e-02 , 5.180473558027376e+02), 
             std::complex<double> (-4.672222488715252e-03 , 3.460402923009413e+02), 
             std::complex<double> (-2.080102598231416e-03 , 2.312070734654582e+02), 
             std::complex<double> (-9.231664120304972e-04 , 1.545747381390437e+02), 
             std::complex<double> (-4.063669606927029e-04 , 1.034822572748718e+02), 
             std::complex<double> (-1.751682455950118e-04 , 6.948917145018561e+01), 
             std::complex<double> (-7.156160098429197e-05 , 4.698182129914718e+01), 
             std::complex<double> (-2.512004300949283e-05 , 3.223630121911779e+01), 
             std::complex<double> (-5.470540831325066e-06 , 2.264983957612602e+01), 
             std::complex<double> (-3.458248330082421e-07 , 1.573834973377633e+01), 
             std::complex<double> (-2.589926405543323e-09 , 9.424967102279309e+00), 
             std::complex<double> (-8.348337637419898e-13 , 3.141592710621604e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (1.246531385938415e+06 , 2.956282374766850e+06), 
             std::complex<double> (2.506733065592390e+06 , 1.156096636249635e+06), 
             std::complex<double> (2.088609109973678e+06 , -3.537315818648565e+05), 
             std::complex<double> (1.241272434514447e+06 , -8.705288026002197e+05), 
             std::complex<double> (6.310717591156868e+05 , -8.342515372468245e+05), 
             std::complex<double> (2.985943629114640e+05 , -6.431457804192142e+05), 
             std::complex<double> (1.367492516623627e+05 , -4.567064447380460e+05), 
             std::complex<double> (6.171645828288417e+04 , -3.133208943835331e+05), 
             std::complex<double> (2.767110144842042e+04 , -2.117518006956108e+05), 
             std::complex<double> (1.237024616501191e+04 , -1.421642631055339e+05), 
             std::complex<double> (5.522831279230314e+03 , -9.516512203510279e+04), 
             std::complex<double> (2.464288902873484e+03 , -6.362053236621246e+04), 
             std::complex<double> (1.099280272202903e+03 , -4.250732275092662e+04), 
             std::complex<double> (4.903146987665073e+02 , -2.839339982309546e+04), 
             std::complex<double> (2.186850046483490e+02 , -1.896359789719047e+04), 
             std::complex<double> (9.753335246203429e+01 , -1.266489642136532e+04), 
             std::complex<double> (4.349935822891281e+01 , -8.458094667430363e+03), 
             std::complex<double> (1.940040088961305e+01 , -5.648576033440358e+03), 
             std::complex<double> (8.652428552847550e+00 , -3.772275653281644e+03), 
             std::complex<double> (3.858911175449535e+00 , -2.519225004627943e+03), 
             std::complex<double> (1.721039751182882e+00 , -1.682403108226222e+03), 
             std::complex<double> (7.675682990098244e-01 , -1.123551291161626e+03), 
             std::complex<double> (3.423268282506078e-01 , -7.503354719568061e+02), 
             std::complex<double> (1.526712132514274e-01 , -5.010921601438338e+02), 
             std::complex<double> (6.808639022063472e-02 , -3.346407359195942e+02), 
             std::complex<double> (3.036374579458956e-02 , -2.234796271338597e+02), 
             std::complex<double> (1.354158902757430e-02 , -1.492424717676142e+02), 
             std::complex<double> (6.039988646576623e-03 , -9.966358958671667e+01), 
             std::complex<double> (2.694513074849377e-03 , -6.655141240510318e+01), 
             std::complex<double> (1.202270684309551e-03 , -4.443504688348371e+01), 
             std::complex<double> (5.365218146476245e-04 , -2.966035001247701e+01), 
             std::complex<double> (2.395448754963829e-04 , -1.978615789013239e+01), 
             std::complex<double> (1.070874682733313e-04 , -1.318088343454806e+01), 
             std::complex<double> (4.795309788625994e-05 , -8.752786941800668e+00), 
             std::complex<double> (2.151030805034829e-05 , -5.769840862467224e+00), 
             std::complex<double> (9.578239609139342e-06 , -3.749469221810053e+00), 
             std::complex<double> (3.453345215971069e-06 , -2.484392132570860e+00), 
             std::complex<double> (4.072444940948181e-07 , -2.038362597345107e+00), 
             std::complex<double> (5.123229727783894e-09 , -2.000380685759731e+00), 
             std::complex<double> (2.760544492680041e-12 , -2.000000190222195e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 5.458426252894914e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.050950882538181e+07 , 2.005406502269755e+06), 
             std::complex<double> (-8.005420745274840e+06 , 4.807254759041592e+06), 
             std::complex<double> (-4.936254557387851e+06 , 5.422172021946796e+06), 
             std::complex<double> (-2.663175432027383e+06 , 4.681375369207916e+06), 
             std::complex<double> (-1.335063274229627e+06 , 3.572023537257050e+06), 
             std::complex<double> (-6.448038119832833e+05 , 2.570514667934843e+06), 
             std::complex<double> (-3.058391087485874e+05 , 1.799365472380248e+06), 
             std::complex<double> (-1.438199255488981e+05 , 1.243313053870234e+06), 
             std::complex<double> (-6.735725605682332e+04 , 8.538905519435105e+05), 
             std::complex<double> (-3.148654309104434e+04 , 5.847770635030086e+05), 
             std::complex<double> (-1.470552286419278e+04 , 3.999469108944969e+05), 
             std::complex<double> (-6.865252003293586e+03 , 2.733665768552077e+05), 
             std::complex<double> (-3.204415165340534e+03 , 1.867940129522092e+05), 
             std::complex<double> (-1.495552756769021e+03 , 1.276209314940099e+05), 
             std::complex<double> (-6.979689025020568e+02 , 8.718736166367456e+04), 
             std::complex<double> (-3.257331653345278e+02 , 5.956243020759933e+04), 
             std::complex<double> (-1.520144631556065e+02 , 4.068977902915746e+04), 
             std::complex<double> (-7.094262808151849e+01 , 2.779684341628767e+04), 
             std::complex<double> (-3.310764475927235e+01 , 1.898909872540562e+04), 
             std::complex<double> (-1.545055192872436e+01 , 1.297217104912800e+04), 
             std::complex<double> (-7.210215250877744e+00 , 8.861776571039685e+03), 
             std::complex<double> (-3.364597549280587e+00 , 6.053812791942889e+03), 
             std::complex<double> (-1.569960040455586e+00 , 4.135590284920948e+03), 
             std::complex<double> (-7.324865793861162e-01 , 2.825184511724380e+03), 
             std::complex<double> (-3.416945563722042e-01 , 1.930002467565401e+03), 
             std::complex<double> (-1.593466837912913e-01 , 1.318477327204231e+03), 
             std::complex<double> (-7.426773456016582e-02 , 9.007318864942993e+02), 
             std::complex<double> (-3.458198687465494e-02 , 6.153691420655362e+02), 
             std::complex<double> (-1.607983186921572e-02 , 4.204486648841556e+02), 
             std::complex<double> (-7.459731728351108e-03 , 2.873225506878231e+02), 
             std::complex<double> (-3.447373385986109e-03 , 1.964250119542468e+02), 
             std::complex<double> (-1.582531414366485e-03 , 1.343968557071215e+02), 
             std::complex<double> (-7.174781164518479e-04 , 9.212236643358214e+01), 
             std::complex<double> (-3.167681984728014e-04 , 6.339035462458018e+01), 
             std::complex<double> (-1.308773427978208e-04 , 4.398198389333952e+01), 
             std::complex<double> (-4.453472286349769e-05 , 3.103115269401535e+01), 
             std::complex<double> (-8.366223443216935e-06 , 2.237367721577699e+01), 
             std::complex<double> (-3.958304402333969e-07 , 1.572135584013273e+01), 
             std::complex<double> (-2.213725251461898e-09 , 9.424841298880422e+00), 
             std::complex<double> (-5.194500495907792e-13 , 3.141592668564982e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (4.757858459331138e+05 , 1.201304323932755e+06), 
             std::complex<double> (9.964623402468410e+05 , 5.305061412770180e+05), 
             std::complex<double> (8.836873918714785e+05 , -8.309101091543230e+04), 
             std::complex<double> (5.609317465831924e+05 , -3.334541928323549e+05), 
             std::complex<double> (3.032083367062090e+05 , -3.489610656505206e+05), 
             std::complex<double> (1.516819854666219e+05 , -2.833162001470454e+05), 
             std::complex<double> (7.313528024965849e+04 , -2.089256799363761e+05), 
             std::complex<double> (3.465612628387137e+04 , -1.477955826496697e+05), 
             std::complex<double> (1.628912651422027e+04 , -1.026067143666145e+05), 
             std::complex<double> (7.627149464015731e+03 , -7.062178822430955e+04), 
             std::complex<double> (3.564961875371606e+03 , -4.841304736960518e+04), 
             std::complex<double> (1.664899007210490e+03 , -3.312660429421671e+04), 
             std::complex<double> (7.772372518024200e+02 , -2.264718871064117e+04), 
             std::complex<double> (3.627782836840460e+02 , -1.547660858645463e+04), 
             std::complex<double> (1.693136338353727e+02 , -1.057438839966940e+04), 
             std::complex<double> (7.901773837936850e+01 , -7.224311267623453e+03), 
             std::complex<double> (3.687645822131733e+01 , -4.935370772249512e+03), 
             std::complex<double> (1.720965173932001e+01 , -3.371590100640262e+03), 
             std::complex<double> (8.031503359451282e+00 , -2.303275196884474e+03), 
             std::complex<double> (3.748203767107624e+00 , -1.573457533013899e+03), 
             std::complex<double> (1.749231366349186e+00 , -1.074888202991120e+03), 
             std::complex<double> (8.163257877457771e-01 , -7.342957535891645e+02), 
             std::complex<double> (3.809480597010548e-01 , -5.016238760882882e+02), 
             std::complex<double> (1.777662118030407e-01 , -3.426766377143156e+02), 
             std::complex<double> (8.294930280633252e-02 , -2.340933069973269e+02), 
             std::complex<double> (3.870393210137684e-02 , -1.599151804894827e+02), 
             std::complex<double> (1.805652566538811e-02 , -1.092401483936102e+02), 
             std::complex<double> (8.420473299673142e-03 , -7.462039478868475e+01), 
             std::complex<double> (3.924328566548329e-03 , -5.096777439599738e+01), 
             std::complex<double> (1.827505499097473e-03 , -3.480599683173383e+01), 
             std::complex<double> (8.500197486023599e-04 , -2.375969935298117e+01), 
             std::complex<double> (3.946455976746133e-04 , -1.620529705270146e+01), 
             std::complex<double> (1.828759168156255e-04 , -1.103225677647784e+01), 
             std::complex<double> (8.472003119270014e-05 , -7.479680090904045e+00), 
             std::complex<double> (3.940326231905205e-05 , -5.025260861581506e+00), 
             std::complex<double> (1.803614519869294e-05 , -3.327887497132883e+00), 
             std::complex<double> (5.890845876225736e-06 , -2.310264685358018e+00), 
             std::complex<double> (5.063827361806577e-07 , -2.018152736228967e+00), 
             std::complex<double> (4.561826352117045e-09 , -2.000132625052268e+00), 
             std::complex<double> (1.820055215917564e-12 , -2.000000051069919e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 2.550702517738841e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-4.930174595007643e+06 , 8.906123945983183e+05), 
             std::complex<double> (-3.856786835171367e+06 , 2.181816114949857e+06), 
             std::complex<double> (-2.482148829939292e+06 , 2.544738896772400e+06), 
             std::complex<double> (-1.405115674755809e+06 , 2.275778593632792e+06), 
             std::complex<double> (-7.387663093551483e+05 , 1.793131224091675e+06), 
             std::complex<double> (-3.734605338130355e+05 , 1.327339812411121e+06), 
             std::complex<double> (-1.850575188533130e+05 , 9.528317576234050e+05), 
             std::complex<double> (-9.079378741949148e+04 , 6.737880814234124e+05), 
             std::complex<double> (-4.432882114956781e+04 , 4.729842074494265e+05), 
             std::complex<double> (-2.159138740213971e+04 , 3.308413813921391e+05), 
             std::complex<double> (-1.050437545701187e+04 , 2.310143674796951e+05), 
             std::complex<double> (-5.107570502079680e+03 , 1.611727248228040e+05), 
             std::complex<double> (-2.482785844932301e+03 , 1.123998936635827e+05), 
             std::complex<double> (-1.206721733980071e+03 , 7.837068482543784e+04), 
             std::complex<double> (-5.864755628960727e+02 , 5.463856613034698e+04), 
             std::complex<double> (-2.850266897319742e+02 , 3.809118596369446e+04), 
             std::complex<double> (-1.385241020326192e+02 , 2.655459986136343e+04), 
             std::complex<double> (-6.732505605045868e+01 , 1.851186670401653e+04), 
             std::complex<double> (-3.272226359482976e+01 , 1.290501053080275e+04), 
             std::complex<double> (-1.590471081285379e+01 , 8.996331796926745e+03), 
             std::complex<double> (-7.730727388922802e+00 , 6.271511018285853e+03), 
             std::complex<double> (-3.757635405299777e+00 , 4.371988368460032e+03), 
             std::complex<double> (-1.826356354800551e+00 , 3.047799586724509e+03), 
             std::complex<double> (-8.875570327597455e-01 , 2.124688530802362e+03), 
             std::complex<double> (-4.312001441356908e-01 , 1.481177996749026e+03), 
             std::complex<double> (-2.093635235643192e-01 , 1.032584631985095e+03), 
             std::complex<double> (-1.015398924312514e-01 , 7.198752872527600e+02), 
             std::complex<double> (-4.915494178042067e-02 , 5.018986144612410e+02), 
             std::complex<double> (-2.373052610205995e-02 , 3.499697563427933e+02), 
             std::complex<double> (-1.141430363810807e-02 , 2.440955463542509e+02), 
             std::complex<double> (-5.464920162706266e-03 , 1.703435249312997e+02), 
             std::complex<double> (-2.600067462292022e-03 , 1.190086518281665e+02), 
             std::complex<double> (-1.222464485916535e-03 , 8.333665338763376e+01), 
             std::complex<double> (-5.581337382325661e-04 , 5.863631192932306e+01), 
             std::complex<double> (-2.349494487813428e-04 , 4.166212288857597e+01), 
             std::complex<double> (-7.721358338537288e-05 , 3.015006071783224e+01), 
             std::complex<double> (-1.231637423407742e-05 , 2.220420220339163e+01), 
             std::complex<double> (-4.391306357772723e-07 , 1.571363395468869e+01), 
             std::complex<double> (-1.865670154700320e-09 , 9.424798723514536e+00), 
             std::complex<double> (-3.589819255287683e-13 , 3.141592657496090e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (2.005193961994937e+05 , 5.368967287834879e+05), 
             std::complex<double> (4.348508784911602e+05 , 2.613421677037183e+05), 
             std::complex<double> (4.073271579442536e+05 , -1.014528678277312e+04), 
             std::complex<double> (2.743122433151479e+05 , -1.374604130298962e+05), 
             std::complex<double> (1.568032983927115e+05 , -1.579952321612911e+05), 
             std::complex<double> (8.255061799191230e+04 , -1.350862380378558e+05), 
             std::complex<double> (4.172060763041765e+04 , -1.033551459827106e+05), 
             std::complex<double> (2.066641886699318e+04 , -7.529131609292091e+04), 
             std::complex<double> (1.013724643425412e+04 , -5.360674414854537e+04), 
             std::complex<double> (4.948784662546389e+03 , -3.775325697210731e+04), 
             std::complex<double> (2.410274668871462e+03 , -2.644883673064380e+04), 
             std::complex<double> (1.172581561484868e+03 , -1.848221639779480e+04), 
             std::complex<double> (5.701384432699806e+02 , -1.289928872520819e+04), 
             std::complex<double> (2.771408365731337e+02 , -8.997406917602681e+03), 
             std::complex<double> (1.346987559495062e+02 , -6.273972563569732e+03), 
             std::complex<double> (6.546364866130867e+01 , -4.374279045122077e+03), 
             std::complex<double> (3.181472865977212e+01 , -3.049583449110763e+03), 
             std::complex<double> (1.546180611535955e+01 , -2.125984214483195e+03), 
             std::complex<double> (7.514579184335930e+00 , -1.482082724595209e+03), 
             std::complex<double> (3.652298515360992e+00 , -1.033192751513496e+03), 
             std::complex<double> (1.775211712738303e+00 , -7.202585951329638e+02), 
             std::complex<double> (8.628959985261258e-01 , -5.021048281924063e+02), 
             std::complex<double> (4.194526383070269e-01 , -3.500251008453313e+02), 
             std::complex<double> (2.038982737415696e-01 , -2.440069713726977e+02), 
             std::complex<double> (9.911655548984351e-02 , -1.700991182751359e+02), 
             std::complex<double> (4.817947036381683e-02 , -1.185756387351027e+02), 
             std::complex<double> (2.341375779258309e-02 , -8.265625933722868e+01), 
             std::complex<double> (1.137013204115651e-02 , -5.761412185658138e+01), 
             std::complex<double> (5.513422775337461e-03 , -4.015376253881021e+01), 
             std::complex<double> (2.666972343186617e-03 , -2.797744576939514e+01), 
             std::complex<double> (1.285914200590372e-03 , -1.948276915915976e+01), 
             std::complex<double> (6.184050788925227e-04 , -1.355172852903723e+01), 
             std::complex<double> (2.975845442855294e-04 , -9.403442074032634e+00), 
             std::complex<double> (1.440608827015309e-04 , -6.491176712996634e+00), 
             std::complex<double> (7.056435232467201e-05 , -4.432329864276197e+00), 
             std::complex<double> (3.341080809385670e-05 , -2.990823380388287e+00), 
             std::complex<double> (9.616118601422761e-06 , -2.189864930972624e+00), 
             std::complex<double> (6.028777387097391e-07 , -2.008183316055813e+00), 
             std::complex<double> (3.999960502862214e-09 , -2.000045049887408e+00), 
             std::complex<double> (1.243766270842621e-12 , -2.000000013597790e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 1.290546341027422e+06;
             std::complex<double> zvec1[] = {
             std::complex<double> (-2.502765248999940e+06 , 4.292332843960182e+05), 
             std::complex<double> (-2.003326421972529e+06 , 1.071454020727822e+06), 
             std::complex<double> (-1.338659258322089e+06 , 1.287200834049083e+06), 
             std::complex<double> (-7.910912521178401e+05 , 1.188412953890208e+06), 
             std::complex<double> (-4.342649931655817e+05 , 9.644918009649816e+05), 
             std::complex<double> (-2.288324369304984e+05 , 7.329354018854547e+05), 
             std::complex<double> (-1.179944225906926e+05 , 5.385868176037319e+05), 
             std::complex<double> (-6.016333578097559e+04 , 3.890816859074764e+05), 
             std::complex<double> (-3.050089300979179e+04 , 2.786618189732523e+05), 
             std::complex<double> (-1.541806340674090e+04 , 1.987104100739138e+05), 
             std::complex<double> (-7.782317924282151e+03 , 1.413866083152952e+05), 
             std::complex<double> (-3.925248858357621e+03 , 1.004879595754572e+05), 
             std::complex<double> (-1.979086822803152e+03 , 7.138004178477790e+04), 
             std::complex<double> (-9.976616590124471e+02 , 5.068939212879608e+04), 
             std::complex<double> (-5.028805454733488e+02 , 3.599114187100624e+04), 
             std::complex<double> (-2.534727979070896e+02 , 2.555306708069782e+04), 
             std::complex<double> (-1.277589638796787e+02 , 1.814156563984228e+04), 
             std::complex<double> (-6.439381958124774e+01 , 1.287948955138966e+04), 
             std::complex<double> (-3.245529386317756e+01 , 9.143631015748384e+03), 
             std::complex<double> (-1.635741123564978e+01 , 6.491378530304735e+03), 
             std::complex<double> (-8.243956194066978e+00 , 4.608445687693074e+03), 
             std::complex<double> (-4.154832445243196e+00 , 3.271689408029632e+03), 
             std::complex<double> (-2.093942435388286e+00 , 2.322686574716394e+03), 
             std::complex<double> (-1.055246730662401e+00 , 1.648965485093894e+03), 
             std::complex<double> (-5.317488508679939e-01 , 1.170678422233750e+03), 
             std::complex<double> (-2.679389925487527e-01 , 8.311394422823523e+02), 
             std::complex<double> (-1.350208595614516e-01 , 5.901066656600286e+02), 
             std::complex<double> (-6.807054249850865e-02 , 4.190130531100916e+02), 
             std::complex<double> (-3.436356661500544e-02 , 2.975807024664036e+02), 
             std::complex<double> (-1.739197353339325e-02 , 2.114175865143284e+02), 
             std::complex<double> (-8.827640522209079e-03 , 1.503119386905098e+02), 
             std::complex<double> (-4.481563739332731e-03 , 1.070222323054811e+02), 
             std::complex<double> (-2.253748345989574e-03 , 7.641963076905652e+01), 
             std::complex<double> (-1.093603267092463e-03 , 5.488155783339397e+01), 
             std::complex<double> (-4.768397554516450e-04 , 3.986045466394683e+01), 
             std::complex<double> (-1.506747084130995e-04 , 2.951642374143406e+01), 
             std::complex<double> (-2.013620258788699e-05 , 2.210476653929401e+01), 
             std::complex<double> (-5.552313613821935e-07 , 1.571027942234789e+01), 
             std::complex<double> (-1.926323040238913e-09 , 9.424784640021977e+00), 
             std::complex<double> (-3.009536009665105e-13 , 3.141592654602796e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (9.194032855037763e+04 , 2.601578403340305e+05), 
             std::complex<double> (2.054760529120878e+05 , 1.371438373866472e+05), 
             std::complex<double> (2.019780009585420e+05 , 7.919171480510266e+03), 
             std::complex<double> (1.434635054299100e+05 , -6.000910599486931e+04), 
             std::complex<double> (8.630916969145439e+04 , -7.641489229142034e+04), 
             std::complex<double> (4.762473385929766e+04 , -6.883488486112833e+04), 
             std::complex<double> (2.513411361413908e+04 , -5.460905206184494e+04), 
             std::complex<double> (1.296642178456633e+04 , -4.092467458172722e+04), 
             std::complex<double> (6.612492775672838e+03 , -2.984423565344748e+04), 
             std::complex<double> (3.352544120517270e+03 , -2.147372364832832e+04), 
             std::complex<double> (1.694740554350311e+03 , -1.534797098873945e+04), 
             std::complex<double> (8.554348434150020e+02 , -1.093300913774342e+04), 
             std::complex<double> (4.314660545210022e+02 , -7.774942791214966e+03), 
             std::complex<double> (2.175421210943460e+02 , -5.524419367812604e+03), 
             std::complex<double> (1.096627197460616e+02 , -3.923654369504035e+03), 
             std::complex<double> (5.527624014365677e+01 , -2.786130441099687e+03), 
             std::complex<double> (2.786165047015301e+01 , -1.978176201054154e+03), 
             std::complex<double> (1.404343366449573e+01 , -1.404444668016473e+03), 
             std::complex<double> (7.078373224684772e+00 , -9.970851771857878e+02), 
             std::complex<double> (3.567616903101862e+00 , -7.078705806875205e+02), 
             std::complex<double> (1.798069886400765e+00 , -5.025419344202305e+02), 
             std::complex<double> (9.062097114755603e-01 , -3.567702274479664e+02), 
             std::complex<double> (4.567286625107881e-01 , -2.532810169071370e+02), 
             std::complex<double> (2.301938531609396e-01 , -1.798098374108593e+02), 
             std::complex<double> (1.160111695895787e-01 , -1.276493995298795e+02), 
             std::complex<double> (5.845649695944904e-02 , -9.061784493274168e+01), 
             std::complex<double> (2.944895517478061e-02 , -6.432625133024436e+01), 
             std::complex<double> (1.482962918015422e-02 , -4.565856346717736e+01), 
             std::complex<double> (7.464831113809746e-03 , -3.240229434032441e+01), 
             std::complex<double> (3.761418188112178e-03 , -2.298628697264040e+01), 
             std::complex<double> (1.902835746509704e-03 , -1.629447729925296e+01), 
             std::complex<double> (9.702052938189658e-04 , -1.153353562848732e+01), 
             std::complex<double> (5.012052927153164e-04 , -8.138637262489231e+00), 
             std::complex<double> (2.638929596206073e-04 , -5.706366558981996e+00), 
             std::complex<double> (1.417189155002111e-04 , -3.951004496819755e+00), 
             std::complex<double> (7.010775597111909e-05 , -2.721903857215691e+00), 
             std::complex<double> (1.727579999046553e-05 , -2.110614293887127e+00), 
             std::complex<double> (8.049919651054121e-07 , -2.003532197178765e+00), 
             std::complex<double> (4.230943854817009e-09 , -2.000014967475692e+00), 
             std::complex<double> (1.094499551766199e-12 , -2.000000003594051e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

        { 
            int np = 40;
            int m = 2;
            double b = 6.988380680095612e+05;
             std::complex<double> zvec1[] = {
             std::complex<double> (-1.359172146074914e+06 , 2.218902434034566e+05), 
             std::complex<double> (-1.109828775325050e+06 , 5.629714892325068e+05), 
             std::complex<double> (-7.665856554752493e+05 , 6.942945071882324e+05), 
             std::complex<double> (-4.708340500485759e+05 , 6.597889805515391e+05), 
             std::complex<double> (-2.687890174238547e+05 , 5.502936302554040e+05), 
             std::complex<double> (-1.471107180914611e+05 , 4.285292880669166e+05), 
             std::complex<double> (-7.866763494098514e+04 , 3.218442880336218e+05), 
             std::complex<double> (-4.154701737698323e+04 , 2.371657967194124e+05), 
             std::complex<double> (-2.179834245361820e+04 , 1.730353515984738e+05), 
             std::complex<double> (-1.139739780817430e+04 , 1.255912705923701e+05), 
             std::complex<double> (-5.948411704568179e+03 , 9.090883846097489e+04), 
             std::complex<double> (-3.101558818577540e+03 , 6.571110150340381e+04), 
             std::complex<double> (-1.616336520195704e+03 , 4.746259906208562e+04), 
             std::complex<double> (-8.420747888118427e+02 , 3.426870863800703e+04), 
             std::complex<double> (-4.386133876776503e+02 , 2.473758319865785e+04), 
             std::complex<double> (-2.284258426830214e+02 , 1.785548516887275e+04), 
             std::complex<double> (-1.189451942485062e+02 , 1.288731690608702e+04), 
             std::complex<double> (-6.192782050055680e+01 , 9.301245435648771e+03), 
             std::complex<double> (-3.223786217893775e+01 , 6.712950532610028e+03), 
             std::complex<double> (-1.678095837619314e+01 , 4.844877806094157e+03), 
             std::complex<double> (-8.735467064786857e+00 , 3.496641695320403e+03), 
             std::complex<double> (-4.548046012081870e+00 , 2.523594992618343e+03), 
             std::complex<double> (-2.368640642876843e+00 , 1.821334752290886e+03), 
             std::complex<double> (-1.234219488020846e+00 , 1.314509180530722e+03), 
             std::complex<double> (-6.434450266692811e-01 , 9.487355431598161e+02), 
             std::complex<double> (-3.354765466493912e-01 , 6.847657657872270e+02), 
             std::complex<double> (-1.747393013299599e-01 , 4.942752044552273e+02), 
             std::complex<double> (-9.077609439692873e-02 , 3.568234393772188e+02), 
             std::complex<double> (-4.691307387069280e-02 , 2.576611456902768e+02), 
             std::complex<double> (-2.400612203371580e-02 , 1.861476668902309e+02), 
             std::complex<double> (-1.205645833691592e-02 , 1.346094779596603e+02), 
             std::complex<double> (-5.857146461554719e-03 , 9.751721010522917e+01), 
             std::complex<double> (-2.690534497031602e-03 , 7.089331505888089e+01), 
             std::complex<double> (-1.120667793503100e-03 , 5.188689906688285e+01), 
             std::complex<double> (-3.838190108324714e-04 , 3.846186039668368e+01), 
             std::complex<double> (-8.267307848528107e-05 , 2.907134377774091e+01), 
             std::complex<double> (-6.151535848098680e-06 , 2.204915069897227e+01), 
             std::complex<double> (-7.822069571952862e-08 , 1.570887925734013e+01), 
             std::complex<double> (-1.385167044514575e-10 , 9.424780073214315e+00), 
             std::complex<double> (-3.993187134695477e-14 , 3.141592653851127e+00), 
            };
            std::vector< std::complex<double> >zvec (zvec1, zvec1 + sizeof(zvec1)/ sizeof(std::complex<double>));
             std::complex<double> wvec1[] = {
             std::complex<double> (4.531832966708455e+04 , 1.351024653375058e+05), 
             std::complex<double> (1.039632823776562e+05 , 7.610876939306529e+04), 
             std::complex<double> (1.066464408204632e+05 , 1.058135927227278e+04), 
             std::complex<double> (7.948171639602531e+04 , -2.732986374943379e+04), 
             std::complex<double> (5.011144763138166e+04 , -3.905856053208002e+04), 
             std::complex<double> (2.887782913283473e+04 , -3.710378091193706e+04), 
             std::complex<double> (1.586265196640482e+04 , -3.051044613177179e+04), 
             std::complex<double> (8.495553792061870e+03 , -2.350431354118958e+04), 
             std::complex<double> (4.489901910371671e+03 , -1.753816383079628e+04), 
             std::complex<double> (2.356493095971274e+03 , -1.287695899479721e+04), 
             std::complex<double> (1.232322863094649e+03 , -9.376627719296965e+03), 
             std::complex<double> (6.432281766836850e+02 , -6.798625384398039e+03), 
             std::complex<double> (3.354101599974503e+02 , -4.918481341627487e+03), 
             std::complex<double> (1.748056468763063e+02 , -3.554182436618822e+03), 
             std::complex<double> (9.107546029061420e+01 , -2.566775054026677e+03), 
             std::complex<double> (4.744182710202697e+01 , -1.853106608307326e+03), 
             std::complex<double> (2.470913325759803e+01 , -1.337649969309736e+03), 
             std::complex<double> (1.286747439148689e+01 , -9.654896097109040e+02), 
             std::complex<double> (6.699602087418877e+00 , -6.968400460744324e+02), 
             std::complex<double> (3.487466729788453e+00 , -5.029307492645578e+02), 
             std::complex<double> (1.815148146863408e+00 , -3.629758106399501e+02), 
             std::complex<double> (9.447167797920660e-01 , -2.619652136781965e+02), 
             std::complex<double> (4.917018333691266e-01 , -1.890628441161611e+02), 
             std::complex<double> (2.559860843198664e-01 , -1.364469423409842e+02), 
             std::complex<double> (1.333751327278221e-01 , -9.847195708006571e+01), 
             std::complex<double> (6.957555426445142e-02 , -7.106321080662359e+01), 
             std::complex<double> (3.633432104048020e-02 , -5.127978530122873e+01), 
             std::complex<double> (1.898277389764431e-02 , -3.699893350812071e+01), 
             std::complex<double> (9.915188283493496e-03 , -2.668828273779921e+01), 
             std::complex<double> (5.176587200043455e-03 , -1.924140442852947e+01), 
             std::complex<double> (2.695475013557251e-03 , -1.385907231432522e+01), 
             std::complex<double> (1.389246217167925e-03 , -9.963394267487034e+00), 
             std::complex<double> (6.998510096477988e-04 , -7.135551752454193e+00), 
             std::complex<double> (3.385575403274995e-04 , -5.070916359163551e+00), 
             std::complex<double> (1.508153916674418e-04 , -3.554281365658535e+00), 
             std::complex<double> (5.103349355169122e-05 , -2.510253784882813e+00), 
             std::complex<double> (6.634346172152688e-06 , -2.061258155545494e+00), 
             std::complex<double> (1.340336457744053e-07 , -2.001466899739071e+00), 
             std::complex<double> (3.153717353813765e-10 , -2.000004875410597e+00), 
             std::complex<double> (1.206047008490024e-13 , -2.000000000943698e+00), 
            };
            std::vector< std::complex<double> >wvec (wvec1, wvec1 + sizeof(wvec1)/ sizeof(std::complex<double>));
            std::vector< std::complex<double> >wvec2;
            std::vector< std::complex<double> >wvec3;
            wvec2.push_back(NULL);
            wvec3.push_back(NULL);
            fweight.push_back(wvec2);
            eweight.push_back(wvec3);
            method.push_back(m);
            beta.push_back(b);
            numPole.push_back(np);
            zshift.push_back(zvec);
            zweight.push_back(wvec);
        } 

   };
   poleClass::~poleClass(){};

   // give me a set of parameter, return a set of pole and zshift.
   bool poleClass::getPole ( int inputMethod, int inputPole, double inputBeta,  std::vector< std::complex<double> > &out_zshift, std::vector < std::complex<double> > &out_zweight ){
        for ( int i = 0; i < method.size(); i++ ){
            if(inputMethod == method[i] ){ 
               if(numPole[i] >= inputPole && beta[i] >= inputBeta){
                 if(i < method.size() - 1)
                   if((numPole[i+1] >= inputPole) && (beta[i+1] >= inputBeta)){
                      // assume the numPole and beta are in assending or desending order 
                      if( (numPole[i+1] < numPole[i]) || (beta[i+1] < beta[i]))
                          continue;
                   }
                   out_zshift = zshift[i];
                   out_zweight = zweight[i];
                   return true;
               }
            }
        }
        return false;
   }
   

