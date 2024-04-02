function [SA_LOBE, SA_lobe] = SA_total(u)

global Lmin Lmax

M       = u(1);                                                            % ratio of distance of max.section from the nose to length
r0      = u(2);                                                            % radius of curvature at nose of hull
r1      = u(3);                                                            % radius of curvature at tail of hull
cp      = u(4);                                                            % prismatic coefficient
l2d     = u(5);                                                            % fineness ratio
L       = Lmin + u(6).*Lmax;                                               % Length of the airship (m)
D       = L./l2d;                                                          % Width of the single lobe (m)
f       = u(11).*D;                                                        % outer lobe disance (m)
g       = u(12).*(D/2);                                                    % relative distance between lobes in vertical direction (m)

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

m = 201;                                                                   % Along longitudinal
n = 201;                                                                   % Along Circumferential

%% Geometry

[~,v1] = ndgrid(linspace(pi,-pi, m),linspace(0,L,n));
r_p1   = D*sqrt(a(1)*(v1/L)+a(2)*(v1/L).^2+a(3)*(v1/L).^3+a(4)*(v1/L).^4+a(5)*(v1/L).^5+a(6)*(v1/L).^6);
r_p1   = real(r_p1);
r_p    = zeros(n,1);

for i = 1:m
    r_p(i) = r_p1(1,i);
end

y11 = v1';

m1=1;
m2=m;
m3=m1;
m4=m2;
Angle1=90;

%%

x11 = zeros(m2,m2);
z11 = zeros(m2,m2);
tt3 = zeros(1,m2);
tt22 = zeros(1,m2);
tt2 = zeros(1,m2);
for itr=m1:m2
    tt1=atand(f/g);
    if r_p(itr)+r_p(itr)<=sqrt(f^2+g^2) || r_p(itr)==0 || r_p(itr)==0
        tt3(itr)=180;
    else
        tt22(itr)=((r_p(itr)^2+f^2+g^2-r_p(itr)^2)/(2*r_p(itr)*sqrt(f^2+g^2)));
        if tt22(itr)<-1
            tt22(itr)=0;
        end
        tt2(itr)=acosd(tt22(itr));
        tt3(itr)=360-(2*tt1+2*tt2(itr));
        if tt3(itr)>180
            tt3(itr)=180;
        end
    end
    theta11 = linspace(90-tt3(itr)/2,90+tt3(itr)/2,m2-m1+1);
    count=1;

    for jtr=m1:m2
        x11(itr,jtr)=r_p(itr)*cosd(theta11(count));
        z11(itr,jtr)=r_p(itr)*sind(theta11(count));
        count=count+1;
    end
end
AREA_E11=0;
dw11 = zeros(m2,m2);
dw22 = zeros(m2,m2);
dAreaE11 = zeros(m2,m2);
for itr = m1:m2-1
    for jtr = m1:m2-1
        dw11(itr,jtr) = sqrt((x11(itr+1,jtr)-x11(itr,jtr))^2 + (y11(itr+1,jtr)-y11(itr,jtr))^2 + (z11(itr+1,jtr)-z11(itr,jtr))^2);
        dw22(itr,jtr) = sqrt((x11(itr,jtr+1)-x11(itr,jtr))^2 + (y11(itr,jtr+1)-y11(itr,jtr))^2 + (z11(itr,jtr+1)-z11(itr,jtr))^2);
        dAreaE11(itr,jtr)  = dw11(itr,jtr).*dw22(itr,jtr);
        AREA_E11 = AREA_E11+dAreaE11(itr,jtr);
    end
end
A11 = AREA_E11;

%%

x22 = zeros(m4,m4);
z22 = zeros(m4,m4);
tt5 = zeros(1,m4);
tt6 = zeros(1,m4);
tt66 = zeros(1,m4);
tt55 = zeros(1,m4);
for itr=m3:m4
    tt4=atand(g/f);
    if r_p(itr)+r_p(itr)<=sqrt(f^2+g^2) || r_p(itr)==0 || r_p(itr)==0
        tt5(itr)=90;
    else
        tt6(itr)=((r_p(itr)^2+f^2+g^2-r_p(itr)^2)/(2*r_p(itr)*sqrt(f^2+g^2)));
        if tt6(itr)>1
            tt6(itr)=0;
            tt66(itr)=acosd(f/r_p(itr));
        else
            tt55(itr)=acosd(tt6(itr));
            tt66(itr)=tt4+tt55(itr);
        end
        tt5(itr)=90-(tt66(itr));
        if tt5(itr)>90
            tt5(itr)=90;
        end
    end
    theta22=linspace(90-Angle1,90+tt5(itr),m4-m3+1);
    count=1;

    for jtr=m3:m4
        x22(itr,jtr)=r_p(itr)*cosd(theta22(count));
        z22(itr,jtr)=r_p(itr)*sind(theta22(count));
        count=count+1;
    end
end
AREA_E22=0;
dw33 = zeros(m4,m4);
dw44 = zeros(m4,m4);
dAreaE22 = zeros(m4,m4);
for itr = m3:m4-1
    for jtr = m3:m4-1
        dw33(itr,jtr) = sqrt((x22(itr+1,jtr)-x22(itr,jtr))^2 + (y11(itr+1,jtr)-y11(itr,jtr))^2 + (z22(itr+1,jtr)-z22(itr,jtr))^2);
        dw44(itr,jtr) = sqrt((x22(itr,jtr+1)-x22(itr,jtr))^2 + (y11(itr,jtr+1)-y11(itr,jtr))^2 + (z22(itr,jtr+1)-z22(itr,jtr))^2);
        dAreaE22(itr,jtr)  = dw33(itr,jtr).*dw44(itr,jtr);
        AREA_E22 = AREA_E22+dAreaE22(itr,jtr);
    end
end
A22 = AREA_E22;

%% Wetted Area of the envelope

SA_LOBE = 2*(A11+(2*A22));                                                 % Wetted area of the envelope (m^2)
SA_lobe = SA_LOBE*10.7639;                                                 % Wetted area of the envelope (ft^2)
