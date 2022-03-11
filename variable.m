 global A B C  a1 a2 a3 R1 R2 beta_0 n0 eta I_input I_input0 N_TE0 N_TM0  confin_TE confin_TM ...
     dN_TEbar_dn dN_TMbar_dn alpha_1 c h alpha_2 q d lambda0 n_p  confin_optional Nz  u Nm tolerence S1 S2
 u = symunit;
 A=1.0e8;%*[u.s^-1];
B=2.5e-17;%*[u.m^3]*[u.s^-1];
C=9.4e-41;%*[u.m^6]*[u.s^-1];

a1=2.5e-20;
a2=1.5e19;
a3=2.7e-32;
R1=5.0e-5;
R2=5.0e-5;
beta_0=1.0e-4;
n0=1.1e24;% taghir dadam 1.1e24;
eta=0.6; %60 darsad   carrier injection efficiency
N_TE0=4.03;
N_TM0=3.95;
confin_TE=0.18;
confin_TM=0.12;
confin_optional=0.465;
dN_TEbar_dn=-1.44e-26;
dN_TMbar_dn=-1.20e-26;
alpha_1=3000;%*[u.m^1];
alpha_2=2.7e-21;%unitConvert(2.7e-21*[u.m^2],u.m^-1);
q=1.602e-19;
Nz=40; % number of spatial divisions
Nm=40;
tolerence = 0.1;
lambda0 = 1.3e-6;  % converting to meters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%from universal And soa constants  M.connelly (Wideband Semiconductor Optical Amplifier  Steady-State Numerical Model)
d = 0.4e-6; % m active layer thickness

% dz=(7e10-6);

n_p=2.66e24; % from refrence 14 N_th
 c = 3e8; % speed light
 h = 6.63e-34; % Planck constant (J-s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%polaristion ellepticity
S1=[1 0.68 -0.71 0.03];
S2=[1 0.43 0.42 -0.80];
N_TE0=4.03;
N_TM0=3.95;
confin_TE=0.18;
confin_TM=0.12;

dN_TEbar_dn=-1.44e-26;
dN_TMbar_dn=-1.20e-26;
