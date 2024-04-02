function [] = postprocess()

global u storage Rho_fc
global Lmin Lmax Amin Amax IDAY XLONG XLAT YAW M_pay fbest Cdlift_yes
global Rho_batt P_pay k_pur Rho_en N_lobes Filename Destination ETA_SC noe


M       = u(:,1);                                                          % ratio of distance of max.section from the nose to length
r0      = u(:,2);                                                          % radius of curvature at nose of hull
r1      = u(:,3);                                                          % radius of curvature at tail of hull
cp      = u(:,4);                                                          % prismatic coefficient
l2d     = u(:,5);                                                          % fineness ratio
L       = Lmin + u(:,6)*Lmax;                                              % Length of the airship (m)
XS      = u(:,7)*(L/2);                                                    % Starting point of the solar array (m)
XF      = (L/2) + u(:,8)*(L/2);                                            % Ending point of the array (m)
ALT     = 15 + 5.*u(:,9);                                                  % Operating altitude (km)
Angle   = Amin + Amax*u(:,10);                                             % Angle of array (deg)
D       = L./l2d;                                                          % Width of individual lobe (m)
f       = u(:,11)*D;                                                       % outer lobe distance (m)
g       = u(:,12).*(D/2);                                                  % relative distance between lobes in vertical direction (m)

Desi    = [M r0 r1 cp l2d u(:,11) u(:,12)];

%% Output

Array_length = XF - XS;                                                    % Length of solar array

%% Input Parameters

H         = ALT;                                                           % Altitude (km)
E_pack    = 0.95;                                                          % Packing eff
E_comp    = 0.95;                                                          % Component eff
Eta_conv  = 0.9;                                                           % convert eff
Eta_gear  = 0.98;                                                          % gear Eff
Eta_prop  = 0.72;                                                          % Propulsion Eff
Rho_prop  = 440;                                                           % Power density of propulsion system(W/kg)
Rho_pv    = 0.3;                                                           % Specific mass of the solar array(kg/m^2)
gravity   = 9.8065;                                                        % Acceration due to gravity
PHI       = XLAT;                                                          % Latitude Angle

%% Gertler Shape Coefficients

A = [   1 0 0 0 0 0 
        1 1 1 1 1 1 
     M M^2 M^3 M^4 M^5 M^6
     1 2*M 3*M^2 4*M^3 5*M^4 6*M^5
     1 2 3 4 5 6 
     1/2 1/3 1/4 1/5 1/6 1/7];
 
B = [2*r0 0 1/4 0 -2*r1 1/4*cp]';
a = A\B;
a = round(a*10000)/10000;

%% Meshing Parametrs

m = noe;                                                                   % Along longitudinal
n = noe;                                                                   % Along Circumferential
p = noe;

%% Geometry

[u1,v1] = ndgrid(linspace(pi,-pi, m),linspace(0,L,n));
r_p1    = D*sqrt(a(1)*(v1/L)+a(2)*(v1/L).^2+a(3)*(v1/L).^3+a(4)*(v1/L).^4+a(5)*(v1/L).^5+a(6)*(v1/L).^6);
r_p1    = real(r_p1);
r_p     =  zeros(n,1);

for i = 1:m
    r_p(i) = r_p1(1,i);
end

n1  = round(XS*(m/L));
n2  = round(XF*(m/L));
x1  = v1';
x2  = v1';
x3  = v1';
y1  = zeros(m,m);
y2  = zeros(m,m);
y3  = zeros(m,m);
z1  = zeros(m,m);
z2  = zeros(m,m);
z3  = zeros(m,m);

%% Volume of the envelope

[Volume] = Volume_gertler(L,u,N_lobes);
Width    = D+(2*f);

%% Creating grid dimensions for central lobe 1

t2  = zeros(n2,1);
t3  = zeros(n2,1);
t22 = zeros(n2,1);

if n1 ==0
    n1 = 1;
end

for itr = n1:n2
    t1  = atand(f/g);
    if r_p(itr)+r_p(itr)<=sqrt(f^2+g^2) || r_p(itr)==0 || r_p(itr)==0
        t3(itr)=180;
    else
        t22(itr)=((r_p(itr)^2+f^2+g^2-r_p(itr)^2)/(2*r_p(itr)*sqrt(f^2+g^2)));
        if t22(itr)<-1
            t22(itr)=0;
        end
        t2(itr)=acosd(t22(itr));
        t3(itr)=360-(2*t1+2*t2(itr));
        if t3(itr)>180
            t3(itr)=180;
        end
    end
    theta1 = linspace(90-t3(itr)/2,90+t3(itr)/2,n2-n1+1);
    count=1;
    for jtr=n1:n2
        y1(itr,jtr)=r_p(itr)*cosd(theta1(count));
        z1(itr,jtr)=r_p(itr)*sind(theta1(count));
        count=count+1;
    end
end
AREA_E1=0;
dw1 = zeros(n2,n2);
dw2 = zeros(n2,n2);
dAreaE1 = zeros(n2,n2);
for itr = n1:n2-1
    for jtr = n1:n2-1
        dw1(itr,jtr) = sqrt((x1(itr+1,jtr)-x1(itr,jtr))^2 + (y1(itr+1,jtr)-y1(itr,jtr))^2 + (z1(itr+1,jtr)-z1(itr,jtr))^2);
        dw2(itr,jtr) = sqrt((x1(itr,jtr+1)-x1(itr,jtr))^2 + (y1(itr,jtr+1)-y1(itr,jtr))^2 + (z1(itr,jtr+1)-z1(itr,jtr))^2);
        dAreaE1(itr,jtr)  = dw1(itr,jtr).*dw2(itr,jtr);
        AREA_E1 = AREA_E1+dAreaE1(itr,jtr);
    end
end
A1 = AREA_E1;

%% Creating grid dimensions for outer lobe 1

t5 = zeros(n2,1);
t6 = zeros(n2,1);
t55 = zeros(n2,1);
t66 = zeros(n2,1);
for itr=n1:n2
    t4=atand(g/f);
    if r_p(itr)+r_p(itr)<=sqrt(f^2+g^2) || r_p(itr)==0
        t5(itr)=90;
    else
        t6(itr)=((r_p(itr)^2+f^2+g^2-r_p(itr)^2)/(2*r_p(itr)*sqrt(f^2+g^2)));
        if t6(itr)>1
            t6(itr)=0;
            t66(itr)=acosd(f/r_p(itr));
        else
            t55(itr)=acosd(t6(itr));
            t66(itr)=t4+t55(itr);
        end
        t5(itr)=90-(t66(itr));
        if t5(itr)>90
            t5(itr)=90;
        end
    end
    theta2=linspace(90-Angle,90+t5(itr),n2-n1+1);
    count=1;
    for jtr=n1:n2
        y2(itr,jtr)=r_p(itr)*cosd(theta2(count));
        z2(itr,jtr)=r_p(itr)*sind(theta2(count));
        count=count+1;
    end
end
AREA_E2=0;
dw3 = zeros(n2,n2);
dw4 = zeros(n2,n2);
dAreaE2 = zeros(n2,n2);
for itr = n1:n2-1
    for jtr = n1:n2-1
        dw3(itr,jtr) = sqrt((x2(itr+1,jtr)-x2(itr,jtr))^2 + (y2(itr+1,jtr)-y2(itr,jtr))^2 + (z2(itr+1,jtr)-z2(itr,jtr))^2);
        dw4(itr,jtr) = sqrt((x2(itr,jtr+1)-x2(itr,jtr))^2 + (y2(itr,jtr+1)-y2(itr,jtr))^2 + (z2(itr,jtr+1)-z2(itr,jtr))^2);
        dAreaE2(itr,jtr)  = dw3(itr,jtr).*dw4(itr,jtr);
        AREA_E2 = AREA_E2+dAreaE2(itr,jtr);
    end
end
A2 = AREA_E2;

%% Creating grid dimensions for outer lobe 2

t8 = zeros(n2,1);
t9 = zeros(n2,1);
t88 = zeros(n2,1);
t99 = zeros(n2,1);
for itr=n1:n2
    t7=atand(g/f);
    if r_p(itr)+r_p(itr)<=sqrt(f^2+g^2) || r_p(itr)==0 
        t8(itr)=90;
    else
        t9(itr)=((r_p(itr)^2+f^2+g^2-r_p(itr)^2)/(2*r_p(itr)*sqrt(f^2+g^2)));
        if t9(itr)>1
            t9(itr)=0;
            t99(itr)=acosd(f/r_p(itr));
        else
            t88(itr)=acosd(t9(itr));
            t99(itr)=t7+t88(itr);
        end
        t8(itr)=90-(t99(itr));
        if t8(itr)>90
            t8(itr)=90;
        end
    end
    theta3=linspace(90-t8(itr),90+Angle,n2-n1+1);
    count=1;
    for jtr=n1:n2
        y3(itr,jtr)=r_p(itr)*cosd(theta3(count));
        z3(itr,jtr)=r_p(itr)*sind(theta3(count));
        count=count+1;
    end
end
AREA_E3=0;
dw5 = zeros(n2,n2);
dw6 = zeros(n2,n2);
dAreaE3 = zeros(n2,n2);
for itr = n1:n2-1
    for jtr = n1:n2-1
        dw5(itr,jtr) = sqrt((x3(itr+1,jtr)-x3(itr,jtr))^2 + (y3(itr+1,jtr)-y3(itr,jtr))^2 + (z3(itr+1,jtr)-z3(itr,jtr))^2);
        dw6(itr,jtr) = sqrt((x3(itr,jtr+1)-x3(itr,jtr))^2 + (y3(itr,jtr+1)-y3(itr,jtr))^2 + (z3(itr,jtr+1)-z3(itr,jtr))^2);
        dAreaE3(itr,jtr)  = dw5(itr,jtr).*dw6(itr,jtr);
        AREA_E3 = AREA_E3+dAreaE3(itr,jtr);
    end
end
A3 =AREA_E3;

Area_total = A1+A2+A3;

%% Wetted Area Estimaion

[SA_LOBE,~] = SA_total(u);

%% Wind Speed model

[V_AIR] = wind_speed(IDAY,ALT,XLONG,XLAT);

%% Calculated Atmospheric Parameters

[~, ~, P_OA, rho_OA] = atmosisa(H*1000);                                   % Atmospheric properties at OA

%% Drag Estimation from Surrogate Model

CDhull = CDSurrogate(Desi);
     
%% Density of helium at Operating altitude
% if the density is in slug/ft^3 then multiply by gravity of 32.174 ft/s^2

rho_g   = 0.1692;                                                          % density of helium at sea-level
sigma   = rho_OA/1.225;                                                    % Density ratio
rho_He  = (k_pur*rho_g+(1-k_pur)*1.225)*sigma;                             % density of helium with purity

%% Aerodynamic Lift

Lift        = (rho_OA-rho_He)*gravity*Volume;                              % Buoyancy of Airship
[s_plan,AS] = Proj_Area(a,L,D,f);
AR          = (Width.^2)./s_plan;                                          % Aspect ratio of the envelope
NL          = 2.4;                                                         % Parameter based on volume and surface area
k_lift      = (-0.0145*(1/AR)^4 + 0.182*(1/AR)^3 - 0.514*(1/AR)^2 +...     % drag due-to-lift factor
              0.838*(1/AR) - 0.053)/NL; 
alpha_deg   = 2;                                                           % Angle of attack in deg
C_L_alpha   = ((2*pi*AR)/(2+sqrt(4+AR^2)))*(pi/180);                       % Lift-curve slope in 1/deg
C_LV        = C_L_alpha.*alpha_deg*NL;                                     % Lift coefficient
% C_LV        = CLSurrogate(Desi);                                           % Lift coefficient
L_aero      = 0.5.*rho_OA.*(V_AIR).^2.*(Volume)^(2/3).*C_LV;               % Aerodynamic lift
Cd_lift     = k_lift*(C_LV)^2*Cdlift_yes;                                  % Drag coefficient due to lift

%% Total drag

N        = 1.1;                                                            % Ratio of CD total/CDhull
CD_total = (N*CDhull)+Cd_lift;                                             % total coeff of drag for airship
Dtotal   = 0.5*rho_OA*V_AIR^2*(CD_total)*Volume^(2/3);                     % total drag of airship
P_thrust = (Dtotal*V_AIR)/(Eta_prop*Eta_gear);                             % Thrust Power required
P_total  = P_thrust + P_pay;                                               % Power Required for payload+drag

%% Input for power estimation

q       = 100;
ISC     = 1367;                                                            % solar constant (W/m^2)
Gamma   = 2*pi*(IDAY-1)/365;                                               % Mean Anomaly
delta   = (180/pi)*(0.006918-0.399912*cos(Gamma)+0.070257*sin(Gamma)-...
           0.006758*cos(2*Gamma)+0.000907*sin(2*Gamma)-0.002697*...
           cos(3*Gamma)+0.001480*sin(3*Gamma));                            % Declination Angle                                                  % Declination Angle
Zi      = Gamma + 0.0334*sin(Gamma) + 0.000349*sin(2*Gamma);               % True Anomaly
Eo      = (1.017 + 0.0174*cos(Zi))^2;                                      % Correction factor for solar radiation at top of atmosphere
I_SUN   = ISC*Eo;                                                          % Extra-terrestrial normal solar intensity
tau_m   = 0.2688;                                                          % AOD at measurement site
Omega_m = 1.4160;                                                          % water vapour column at measurement site
tau     = tau_m*exp(-0.691*H);                                             % AOD at different altitude
Omega   = Omega_m*exp(-0.44*H);                                            % water vapour column at different altitude
wSRH    = -acosd(-tand(PHI)*tand(delta));                                  % Sunrise hour angle
wSSH    = -wSRH;                                                           % Sunset hour angle
tSRH    = wSRH/15+12;                                                      % Sunrise time in hour
tSSH    = wSSH/15+12;                                                      % Sunset time in hour
Nd      = (tSSH-tSRH);                                                     % day time duration in hour
Nt      = 24-Nd;                                                           % night time duration in hour
w       = linspace(wSRH,wSSH,q);                                           % array of hour angle
alpha   = asin(sind(delta)*sind(PHI)+cosd(delta)*cosd(PHI)*cosd(w));       % Solar Elevation ANGLE (Sen, 2008)
m_R     = (sin(alpha) + 0.15.*(3.885 + alpha).^-1.253).^-1;                % Relative air mass
m_A     = m_R*(P_OA/101320);                                               % Absolute air mass
I_DN    = I_SUN.*exp(-(0.103.*m_A.^0.571)-(0.081.*(Omega.*m_R).^0.213)-...
         (tau.^0.91.*m_R.^0.87));                                          % Direct Normal radiation
I_dh    = (0.143 + 0.113.*sin(alpha) - 0.0485.*Omega + tau).*...
         (I_SUN-I_DN).*sin(alpha);                                         % Scattered/diffused radiation
In      = I_DN+I_dh;                                                       % Total Incidence

%% Power Estimation

si = zeros(1,q);                                                           % Azimuth angle                            

for i=1:1:q
    if (w(i)<=-90)
        si(i) =  asind(-cosd(delta).*sind(w(i))./cos(alpha(i)));
    elseif (-90<=w(i)) && (w(i)<0)
        si(i) = 180 - asind(-cosd(delta).*sind(w(i))./cos(alpha(i)));
    elseif (0<=w(i)) && (w(i)<=90)
        si(i) = 180 + asind(cosd(delta).*sind(w(i))./cos(alpha(i))); 
    else (90<w(i)); %#ok<SEPEX,VUNUS>
        si(i) =  360 - asind(cosd(delta).*sind(w(i))./cos(alpha(i)));
    end
end

d_si = zeros(1,q);
for i=1:50
    d_si(i) = si(i+1)-si(i);
end
for i=51:q
    d_si(i) = si(i)-si(i-1);
end
% 
for i=1:1:q
    if (d_si(i)<0 && w(i)<=0)
        si(i) = asind(-cosd(delta).*sind(w(i))./cos(alpha(i)));
    elseif (d_si(i)<=0 && 0<=w(i))
        si(i) = 360 + asind(-cosd(delta).*sind(w(i))./cos(alpha(i)));
    end
end

IsunX = (sind(si).*cos(alpha).*cosd(YAW) - cosd(si).*cos(alpha).*sind(YAW));
IsunY = (cosd(si).*cos(alpha).*cosd(YAW) + sind(si).*cos(alpha).*sind(YAW));
IsunZ = (sin(alpha));

%% For Central Lobe
%----------------------------Surface Normal-------------------------------%

[Nx1, Ny1, Nz1] = surfnorm(x1,y1,z1);
Nx_p1    = zeros(p,p);
Ny_p1    = zeros(p,p);
Nz_p1    = zeros(p,p);
dPOWER1  = zeros(p,p);
dT = (w(2)-w(1))/15;

% -------------------------Power Calculation------------------------------%
POWER1   = 0;
for i = n1:n2-1
    for j = n1:n2-1
        Nx_p1(i,j) = (Nx1(i+1,j+1)+Nx1(i,j))/2;
        Ny_p1(i,j) = (Ny1(i+1,j+1)+Ny1(i,j))/2;
        Nz_p1(i,j) = (Nz1(i+1,j+1)+Nz1(i,j))/2;
        for k=1:q
            ANGLE1 = max((Nx_p1(i,j)*IsunX(k)+Ny_p1(i,j)*IsunY(k)+Nz_p1(i,j)*IsunZ(k)),0);
            dPOWER1(i,j) = E_pack*E_comp*ETA_SC*dAreaE1(i,j).*ANGLE1*dT.*In(k);
            POWER1 = POWER1+dPOWER1(i,j);
        end
    end
end
P_sup1 = POWER1;

%% For Outer Lobe 1
%----------------------------Surface Normal-------------------------------%

[Nx2, Ny2, Nz2] = surfnorm(x2,y2,z2);
Nx_p2    = zeros(p,p);
Ny_p2    = zeros(p,p);
Nz_p2    = zeros(p,p);
dPOWER2  = zeros(p,p);
dT = (w(2)-w(1))/15;

% --------------------------Power Calculation-----------------------------%
POWER2   = 0;
for i = n1:n2-1
    for j = n1:n2-1
        Nx_p2(i,j) = (Nx2(i+1,j+1)+Nx2(i,j))/2;
        Ny_p2(i,j) = (Ny2(i+1,j+1)+Ny2(i,j))/2;
        Nz_p2(i,j) = (Nz2(i+1,j+1)+Nz2(i,j))/2;
        for k=1:q
            ANGLE2 = max((Nx_p2(i,j)*IsunX(k)+Ny_p2(i,j)*IsunY(k)+Nz_p2(i,j)*IsunZ(k)),0);
            dPOWER2(i,j) = E_pack*E_comp*ETA_SC*dAreaE2(i,j).*ANGLE2*dT.*In(k);
            POWER2 = POWER2+dPOWER2(i,j);
        end
    end
end
P_sup2 = POWER2;

%% For Outer Lobe 2
%----------------------------Surface Normal-------------------------------%

[Nx3, Ny3, Nz3] = surfnorm(x3,y3,z3);
Nx_p3    = zeros(p,p);
Ny_p3    = zeros(p,p);
Nz_p3    = zeros(p,p);
dPOWER3  = zeros(p,p);
dT = (w(2)-w(1))/15;

% ---------------------------Power Calculation----------------------------%
POWER3   = 0;
for i = n1:n2-1
    for j = n1:n2-1
        Nx_p3(i,j) = (Nx3(i+1,j+1)+Nx3(i,j))/2;
        Ny_p3(i,j) = (Ny3(i+1,j+1)+Ny3(i,j))/2;
        Nz_p3(i,j) = (Nz3(i+1,j+1)+Nz3(i,j))/2;
        for k=1:q
            ANGLE3 = max((Nx_p3(i,j)*IsunX(k)+Ny_p3(i,j)*IsunY(k)+Nz_p3(i,j)*IsunZ(k)),0);
            dPOWER3(i,j) = E_pack*E_comp*ETA_SC*dAreaE3(i,j).*ANGLE3*dT.*In(k);
            POWER3 = POWER3+dPOWER3(i,j);
        end
    end
end
P_sup3 = POWER3;

%% Total Power supplied

P_sup=P_sup1+P_sup2+P_sup3;

%% Efficiency terms related to battery system and Fuel cell

eta_bat = 0.9;                                                             % Battery Efficiency
eta_IL  = 0.98;                                                            % Battery Input Line Efficiency
eta_OL  = 0.98;                                                            % Battery Output Line Efficiency
eta_T   = 0.98;                                                            % Electricity Transmission Efficiency
eta_RTE = 0.9;                                                             % Round-Trip Efficiency
eta_BS  = eta_bat*eta_IL*eta_OL*eta_T*eta_RTE;                             % Total efficiency of battery system 
eta_fc  = 0.55*eta_T;                                                      % Efficiency of fuel cell

%% Mass Estimation

M_prop   = P_thrust/Rho_prop;                                              % Mass of Proplsion System
MW_he    = 4.003;                                                          % Molecular Mass of Helium
MW_air   = 28.97;                                                          % Molecular Mass of Air
M_he     = rho_OA*(MW_he/MW_air)*Volume;                                   % Mass of Lifting gas
M_en     = 1.2*1.26*Rho_en*SA_LOBE;                                        % Mass of hull
Area_sep = 2*0.75*AS;                                                      % Area of the septum in m^2
M_Sep    = 4*Rho_en*Area_sep*1.06;                                         % Weight of the septum in kg
C_ht    = -0.0051*(10^6./(Volume*35.3147)) + 0.0717;
C_vt    = -0.0049*(10^6./(Volume*35.3147)) + 0.0641;
L_arm   = 0.38*((L*3.28084) - 0.45*(L*3.28084));
S_ht_ft = (C_ht*((Volume*35.3147)^(2/3)*(L*3.28084)))/(L_arm);
S_ht_m  = S_ht_ft*0.092903;
S_vt_ft = (C_vt*((Volume*35.3147)^(2/3)*(L*3.28084)))/(L_arm);
S_vt_m  = S_vt_ft*0.092903;
S_ht_cs = 0.02*S_ht_m;
S_vt_cs = 0.02*S_vt_m;
M_fin   = Rho_en*((S_ht_m+S_vt_m )-(S_ht_cs+S_vt_cs))*1.26*2.36;

qmax    = 0.5*rho_OA*V_AIR^2;

if qmax < 4.882
    F_DY = 0.3906;
elseif qmax > 4.882 && qmax < 48.82
    F_DY = 0.3906*(2-qmax);
else
    F_DY = 3.91;
end

M_CS     = F_DY*(S_ht_cs+S_vt_cs);
M_act    = F_DY*1.15*(S_ht_cs+S_vt_cs);

M_tail   = M_fin+M_CS+M_act;
% A_fin    = 2*0.0121*Volume;                                                % Area of fin
% M_fin    = 1.2*Rho_en*A_fin;                                               % Mass of fin
M_array  = 1.3*Rho_pv*Area_total;                                          % Mass of Solar array
switch storage
    case 1
        In_fac = 1.15;                                                     % Installation factor
        M_stor = In_fac*(P_total*Nt/(Rho_batt*eta_BS));                    % Mass of storage system (Battery)
    case 2
        I_fac  = 1.25;                                                     % Installation factor
        M_stor = I_fac*(P_total*Nt/(Rho_fc*eta_fc));                       % Mass of storage system (Fuel cell)
end
M_energy = M_array+M_stor;                                                 % Mass of energy subsystem
M_cont   = 0.25*(M_en+M_tail+M_energy+M_prop);                             % Mass of other components(kg)
M_struc  = M_he+M_en+M_tail+M_cont;                                        % Mass of structure system(kg)
M_total  = 1.006*(M_pay + M_struc + M_energy + M_prop + M_Sep);            % Total Mass(kg)
P_req    = P_total*Nd+(P_total*Nt/Eta_conv);                               % Energy Required
Weight   = M_total*gravity;

%% Empty weight

W_empty = M_total - M_pay;

%% Constraints

Lift_total = Lift + L_aero;                                                % Total lift of airship (N)
P_extra    = abs(real(P_req - P_sup));                                     % Power extra (N)
Lift_extra = abs(real((M_total*gravity)-Lift_total));                      % Lift extra  (N)

%% Best value of objective function 

MinObj = fbest;
    
%% Output
disp(' ');
disp(('Optimized Design Results:'));
disp(' ');
fprintf('The length of the envelope is estimated to be %3.2f m.\n',L);
fprintf('The width of the envelope is estimated to be %2.2f m.\n',Width);
fprintf('The fineness ratio is estimated to be %1.2f \n',L/Width);
fprintf('The starting point of solar array is estimated to be %2.2f m.\n',XS);
fprintf('The ending point of solar array is estimated to be %2.2f m.\n',XF);
fprintf('The angle of array is estimated to be %2.2f deg.\n',Angle);
fprintf('The operating altitude is estimated to be %2.2f km.\n',ALT);
fprintf('The wind speed at operating altitude is estimated to be %2.4f m/s.\n',V_AIR);
fprintf('The length of the array is estimated to be %2.2f m.\n', Array_length);
fprintf('The area of the array is estimated to be %4.0f m^2.\n', Area_total);
fprintf('The wetted area of the envelope is estimated to be %4.0f m^2.\n', SA_LOBE);
fprintf('The hull drag coefficient is estimated to be %1.4f.\n', CDhull);
fprintf('The total drag coefficient is estimated to be %1.4f.\n', CD_total);
fprintf('The drag is estimated to be %4.2f N.\n', Dtotal);
fprintf('The buoyant lift is estimated to be %5.0f N.\n', Lift);
fprintf('The aerodynamic lift is estimated to be %5.0f N.\n', L_aero);
fprintf('The weight is estimated to be %5.0f N.\n',Weight);
fprintf('The array mass is estimated to be %4.0f kg.\n', M_array);
fprintf('The other component mass is estimated to be %4.0f kg.\n', M_cont);
fprintf('The envelope mass is estimated to be %4.0f kg.\n', M_en);
fprintf('The energy mass is estimated to be %4.0f kg.\n', M_energy);
fprintf('The fin mass is estimated to be %4.0f kg.\n', M_fin);
fprintf('The helium mass is estimated to be %4.0f kg.\n', M_he);
fprintf('The energy storage mass is estimated to be %4.0f kg.\n', M_stor);
fprintf('The structure mass is estimated to be %4.0f kg.\n', M_struc);
fprintf('The total mass is estimated to be %4.0f kg.\n', M_total);
fprintf('The empty weight is estimated to be %4.0f kg.\n', W_empty);
fprintf('The objective value is estimated to be %4.0f.\n', MinObj);
fprintf('The total power is estimated to be %3.0f KW.\n', P_total/1000);
fprintf('The power required is estimated to be %3.0f KW.\n', P_req/1000);
fprintf('The excess power during day is estimated to be %3.0f KW.\n', (P_req-P_total)/1000);
fprintf('The power supplied is estimated to be %3.0f KW.\n', P_sup/1000);
fprintf('The extra power is estimated to be %1.4f kW\n', P_extra/1000);
fprintf('The lift extra is estimated to be %1.4f N.\n', Lift_extra);
fprintf('The volume is estimated to be %6.0f m^3.\n', Volume);

%% Envelope and Array Generation

YY = v1;
XX = r_p1.*sin(u1);
ZZ = r_p1.*cos(u1);

figure,
mesh(XX,YY,ZZ,'FaceColor','g','EdgeColor','none')
hold on
mesh(XX+f,YY,ZZ-g,'FaceColor','g','EdgeColor','none')
mesh(XX-f,YY,ZZ-g,'FaceColor','g','EdgeColor','none')

surf(y1(n1:n2,n1:n2),x1(n1:n2,n1:n2),z1(n1:n2,n1:n2),'FaceColor','[0 0 0.2]','EdgeColor','[0 0 0.2]')
hold on
surf(y2(n1:n2,n1:n2)+f,x2(n1:n2,n1:n2),z2(n1:n2,n1:n2)-g,'FaceColor','[0 0 0.2]','EdgeColor','[0 0 0.2]')
surf(y3(n1:n2,n1:n2)-f,x3(n1:n2,n1:n2),z3(n1:n2,n1:n2)-g,'FaceColor','[0 0 0.2]','EdgeColor','[0 0 0.2]')

%% Fin Design

x1fin = f;
y1fin = 7*L/10;
z1fin = 0;
x2fin = f+((2*L/10)*cosd(60));
y2fin = 24*L/28;
z2fin = (2*L/10)*sind(60);
x3fin = f+((2*L/10)*cosd(60));
y3fin = ((24*L/28)+(L/10));
z3fin = (2*L/10)*sind(60);
x4fin = f;
y4fin = ((24*L/28)+(L/10));
z4fin = 0;
x5fin = x1fin;
y5fin = y1fin;
z5fin = z1fin;

xfin_1 = [x1fin x2fin x3fin x4fin x5fin];
yfin_1 = [y1fin y2fin y3fin y4fin y5fin];
zfin_1 = [z1fin z2fin z3fin z4fin z5fin];

xfin_2 = [-x1fin -x2fin -x3fin -x4fin -x5fin];
yfin_2 = [y1fin y2fin y3fin y4fin y5fin];
zfin_2 = zfin_1;

fill3(xfin_1,yfin_1,zfin_1,[0.5 0.5 0.5]);
fill3(xfin_2,yfin_2,zfin_2,[0.5 0.5 0.5]);

x1fin = f;
y1fin = 7*L/10;
z1fin = 0;
x2fin = f+((2*L/10)*sind(150));
y2fin = 24*L/28;
z2fin = (2*L/10)*cosd(150);
x3fin = f+((2*L/10)*sind(150));
y3fin = ((24*L/28)+(L/10));
z3fin = (2*L/10)*cosd(150);
x4fin = f;
y4fin = ((24*L/28)+(L/10));
z4fin = 0;
x5fin = x1fin;
y5fin = y1fin;
z5fin = z1fin;

xfin_3 = [x1fin x2fin x3fin x4fin x5fin];
yfin_3 = [y1fin y2fin y3fin y4fin y5fin];
zfin_3 = [z1fin z2fin z3fin z4fin z5fin];

xfin_4 = [-x1fin -x2fin -x3fin -x4fin -x5fin];
yfin_4 = [y1fin y2fin y3fin y4fin y5fin];
zfin_4 = zfin_3;

fill3(xfin_3,yfin_3,zfin_3,[0.5 0.5 0.5])
fill3(xfin_4,yfin_4,zfin_4,[0.5 0.5 0.5])

camlight(-26,34);
lighting phong
axis equal
ylim([0 round((L+20)/50)*50])
grid on
box on
rotate3d on;

%% To preserve the size of the image when we save it

width  = 3.350394;   % Width in inches
height = 2.75591;    % Height in inches
alw    = 0.75;       % AxesLineWidth
fsz    = 10;         % Fontsize

pos = get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100 height*100]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)-width)/2;
bottom = (papersize(2)-height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition',myfiguresize);

%% Save the file as PNG

% print('trilobe','-dpng','-r600');

eval(['print -dpng Trilobed_' Filename '.png']);

%% Save the image as pdf file

% save2pdf

%% Save the output to an excel file

result = {'Length of the envelope (m)'      ,L;...
          'Width of the envelope (m)'       ,Width;...
          'Fineness ratio'                  ,L/Width;...
          'Starting point of the array (m)' ,XS;...
          'Ending point of the array (m)'   ,XF;...
          'Angle of the array (deg)'        ,Angle;...
          'Operating altitude (km)'         ,ALT;...
          'Wind speed (m/s)'                ,V_AIR;...
          'Length of the array (m)'         ,Array_length;...
          'Area of the array (m^2)'         ,Area_total;...
          'Volume of the envelope (m^3)'    ,Volume;...
          'Wetted area of envelope (m^2)'   ,SA_LOBE;...
          'Hull drag coefficient'           ,CDhull;...
          'Total drag coefficient'          ,CD_total;...
          'Drag (N)'                        ,Dtotal;...
          'Buoyant lift (N)'                ,Lift_total;...
          'Aerodynamic lift (N)'            ,L_aero;...
          'Weight of the airship (N)'       ,Weight;...
          'Mass of the array (kg)'          ,M_array;...
          'Other masses (kg)'               ,M_cont;...
          'Envelope mass (kg)'              ,M_en;...
          'Energy storage mass (kg)'        ,M_stor;...
          'Fin mass (kg)'                   ,M_fin;...
          'Mass of lifting gas'             ,M_he;...
          'Structural mass (kg)'            ,M_struc;...
          'Total mass (kg)'                 ,M_total;...
          'Empty weight (kg)'               ,W_empty;...
          'Power required (KW)'             ,P_req/1000;...
          'Power supplied (KW)'             ,P_sup/1000;...
          'Extra power (KW)'                ,P_extra/1000;...
          'Lift extra (N)'                  ,Lift_extra};
      
%% Creation of folder to store the data

mkdir(Destination)                                                         % Creates folder
baseFileName = 'Optimized_Results.xls';                                    % File name
fullFileName = fullfile(Destination, baseFileName);
col_header={'Name of the parameter','Value'};                              % Row cell array (for column labels)
copyfile('C:\Users\HP\Desktop\Sanjana Files\template.xls',fullFileName)                                      % Copying the given template
xlswrite(fullFileName,result,'Sheet1','A2');                               % Write data
xlswrite(fullFileName,col_header,'Sheet1','A1');                           % Write column header

A = 'C:\Users\HP\Desktop\Project';
f_raw1 = dir(fullfile(A,'*.png'));    
f_raw2 = dir(fullfile(A,'*.txt'));
movefile(fullfile(f_raw1.folder,f_raw1.name),Destination)
movefile(fullfile(f_raw2.folder,f_raw2.name),Destination)
