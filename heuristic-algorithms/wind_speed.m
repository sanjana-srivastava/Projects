function [V_AIR,Diff] = wind_speed(IDAY,ALT,XLONG,XLAT)

XLST    = 0:1:24;
% V_west  = zeros(length(XLST),1);
% V_north = zeros(length(XLST),1);
V_Air   = zeros(length(XLST),1);

for i = 1:length(XLST)
v_air        = atmoshwm(XLAT, XLONG, ALT*1000,'day',IDAY,'seconds',XLST(i)*3600,'apindex',80,'model','total');
% V_west(i,:)  = v_air(1);
% V_north(i,:) = v_air(2);
V_Air(i,:)   = sqrt(v_air(1)^2+v_air(2)^2);
end

min_V_AIR = min(V_Air(:,1));
V_AIR     = max(V_Air(:,1));
Diff      = V_AIR - min_V_AIR;