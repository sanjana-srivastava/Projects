function [Area_Top, Area_Side] = Proj_Area(a,L,D,f)

 Y = @(X) real(D*sqrt(a(1)*(X/L)+a(2)*(X/L).^2+a(3)*(X/L).^3+a(4)*(X/L).^4+a(5)*(X/L).^5+a(6)*(X/L).^6));
 
 Y1 = @(X) real(-Y(X));
 
 Y2 = @(X) real(Y(X)+f);
  
 Y3 = @(X) real(-Y2(X)+2*f);
 
 Y4 = @(X) real(Y(X)-f);
 
%  Y5 = @(X) real(-Y4(X)-2*f);

X0    = [5 L];
xint  = zeros(2,1);
xint2 = zeros(2,1);
 for k1 = 1:2
    x0 = X0(k1);
    xint(k1)  = fzero(@(X) Y(X)-Y3(X),   x0);
    xint2(k1) = fzero(@(X) Y1(X)-Y4(X), x0);
 end


reqd1   = xint(1);
reqd2   = xint(2);
II      = 4*integral(Y,0,reqd1);
III     = 4*integral(Y,reqd2,L);
IV      = (reqd2 - reqd1)*2*f;

Area_Side =  2*integral(Y,0,L);                                            %Side Projected Area

I = Area_Side;

Area_Top = I+II+III+IV;                                                    %Top Projected Area

