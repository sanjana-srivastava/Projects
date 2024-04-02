function [Volume] = Volume_gertler(L,u,N_lobes)

M       = u(1);                                                            % ratio of distance of max.section from the nose to length
r0      = u(2);                                                            % radius of curvature at nose of hull
r1      = u(3);                                                            % radius of curvature at tail of hull
cp      = u(4);                                                            % prismatic coefficient
l2d     = u(5);                                                            % fineness ratio
D       = L./l2d;                                                          % Width of the single lobe (m)
f       = u(11).*D;                                                        % outer lobe disance (m)
g       = u(12).*(D/2);                                                    % relative distance between lobes in vertical direction (m)

% Gertler Shape Coefficients

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

m = 501;
n = 501;
p = 501;

%% Geometry

[~,v1] = ndgrid(linspace(pi,-pi, m),linspace(0,L,n));
r_p1   = D*sqrt(a(1)*(v1/L)+a(2)*(v1/L).^2+a(3)*(v1/L).^3+a(4)*(v1/L).^4+a(5)*(v1/L).^5+a(6)*(v1/L).^6);
r_p1   = real(r_p1);
r_p    = zeros(n,1);

for i = 1:m
    r_p(i) = r_p1(1,i);
end                                                               
len = sqrt(f^2+g^2);
Area_section = zeros(p,1);
tc = zeros(p,1);
to = zeros(p,1);
j=0;

%% Area of intersection between lobes

for itr=1:p-1
    if r_p(itr)+r_p(itr)<=len || r_p(itr)==0
        Area_section(itr)=0;
        i=itr;
    else
        tc(itr) = 2*(acosd((r_p(itr)^2+len^2-r_p(itr)^2)/(2*r_p(itr)*len)));
        to(itr) = 2*(acosd((r_p(itr)^2+len^2-r_p(itr)^2)/(2*r_p(itr)*len)));
        Area_section(itr) = 2*((0.5*r_p(itr)^2*((tc(itr)*(pi/180))-(sind(tc(itr)))))+(0.5*r_p(itr)^2*((to(itr)*(pi/180))-sind(to(itr)))));
        j=itr-i;
    end
end

%% Intersection volume between lobes

if j==0
    V_section_net=0;
    disp('There is no intersection');
else
    A_sum=sum(Area_section);
    avg_area=A_sum/j;
    l_section=(L)*(j/p);
    V_section_net=avg_area*l_section;
end

%% Volume of single lobe

Vol     = (cp*pi*D^2*L)/4;                                                 % Volume of hull

%% Total volume of the tri-lobed envelope

Volume  = (N_lobes*Vol)-V_section_net;

