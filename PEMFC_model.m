%this model uses equations found in Fuel Cell Fundamentals 3rd Edition for
%PEMFC modeling
%we will use V = E_t - E_o - E_c
%V is cell voltage, E_t is thermo, E_o is ohmimc losses, E_c is cathode overpotential

%syms j j_plot;
j = linspace(0,1,1000);
j_plot = j*1000/5;
%j is operating current density, we assume it is in [A/cm^2];

syms t_m p_a p_c p_SAT D_effa D_effc F T R x_wa x_wd x_o t_A t_C a C z;
%t_m is membrane thickness, p_a is pressure of anode, p_c is pressure of 
%cathode, p_SAT is vapor saturation pressure, D_effa is effective 
%hydrogen/water diffusivity (anode), D_effc is effeective oxygen/water 
%diffusivity (cathode), F is Faraday's constant,T is operating tempurature,
%R is gas constant, x_wa is anode water mole fraction, x_wd is cathode 
%water mole fraction, x_o is oxygen mole fraction at the cathode catalyst 
%layer, t_A is anode thickness, t_C is cathode thickness, 
%a and C are unknown variables to be solved


%please note that the following values are assumptions, based on table
%values found on TABLE 6.5, pg. 222 of Fuel Cell Fundamentals textbook

t_m = 125;
%membrane thickness [um]

t_A = 350;
%anode thickness [um]

t_C = 350;
%cathode thickness[um]

p_a = 3;
%anode pressure [atm]

p_c = 3;
%cathode pressure[atm]

p_SAT = 0.307;
%vapor saturation pressure[atm]

F = 96485;
%Faraday's constant [C/mol]

T = 343;
%FC operating temperature[Kelvin]

R = 8.314;
%gas constant [J/mol*K]

x_wa = 0.1;
x_wd = 0.1;
%water mole fraction for anode and cathode

x_o = 0.19;
%oxygen mole fraction at the cathode catalyst layer

D_effa = 0.149;
%Effective hydrogen/water diffusivity (anode) [cm^2/s]

D_effc = 0.0295;
%Effective oxygen/water diffusivity (cathode) [cm^2/s]


%%%%%%%%%%%%%%%for E_t or theoretical value of fuel cell%%%%%%%%%%%%%%%%%%
E_t = 1.2;
% 1.2 volts, ideal voltage


%%%%%%%%%%%%%%%%for E_o or Ohmc Overpotential%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%v = zeros(1,1000);
for j = 0:0.001:1
    

%lambda_b  
%simplify(14*(p_a/p_SAT)*(x_wa - t_A*0.000001*((a*j*(1/0.0001)*8.314*T)/(2*F*p_a*101325*D_effa*0.0001))))

%lambda_c
%simplify( 10 + 4*(p_c/p_SAT)*(x_wd + t_C*0.000001*(((1+a)*j*(1/0.0001)*R*T)/(2*F*p_c*101325*D_effc*0.0001))))
%equations (6.45) & (6.46) derved in p.222;

%4.4*a + C; from (6.49 p. 223)
%4.4*a + C*exp((0.000598*j*t_m*0.0001)/(0.00000381))); from (6.50 p. 223)

eqns = [14*(p_a/p_SAT)*(x_wa - t_A*0.000001*((a*j*(1/0.0001)*8.314*T)/(2*F*p_a*101325*D_effa*0.0001))) - 4.4*a - C == 0 , 
10 + 4*(p_c/p_SAT)*(x_wd + t_C*0.000001*(((1+a)*j*(1/0.0001)*R*T)/(2*F*p_c*101325*D_effc*0.0001))) - 4.4*a - C*exp((0.000598*j*t_m*0.0001)/(0.00000381)) == 0];
%NOTE: 0.0000001 is the conversion from um to m, (1/0.0001) is conversion from
%cm^2 to m^2, 101325 is for atm to Pa, 0.0001 is converstion from cm^2 to
%m^2)


S = solve(eqns, [a C]);
%type S.a and S.C to return values

%now to solve for conductivity profile of the membrane, sigma(z). Equation
%is based off of(6.40) on p.221, it is based on properties of Nafion. f
%represents sigma(z), further, fInt is ASR or area specific resistance

f = (0.005193*(4.4*S.a + S.C*exp((0.000598*j*z)/(0.00000381)))- 0.00326)*exp(1268*((1/303)-(1/T)));
fInt = int(1/f, z, [0 t_m*0.0001]);
%round(fInt*10^3)/ vpa(10^3);

E_o = j*fInt;
%%%%%%%%%%%%%%%%we have ohmic overpotential [Volts]%%%%%%%%%%%%%%%%%%%

%now we must find E_c or the Cathode overpotential 
%we use equation (6.27) p. 217
%we assume the charge transfer coefficient, alpha, is 0.5
syms j_0;
j_0 = 0.0001;
%this is reference current [A/cm^2], it is a constant and shouldn't have to
%be changed

eqn1 = (((R*T)/(4*0.5*F))*log(j/(j_0*p_c*101300*(x_o-t_C*0.000001*(j*10000*R*T)/(4*f*p_c*101325*D_effc*0.0001)))));
%NOTE: 101300 is conversion from atm to Pa, 0.000001 is to convert um to m,
%10000 is to convert A/cm^2 to A/m^2, 101325 is meant to convert atm to Pa, 
%0.0001 is meant to convert cm^2 to m^2
vpa(eqn1)

V = E_t - E_o - eqn1;
V(j) = j;
end
plot(j_plot, V);