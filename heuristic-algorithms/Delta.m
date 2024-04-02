function [delta egamma] = Delta(Gamma)
egamma       = -0.00167 + 7.36629*cos(Gamma + 1.49831) - 9.92745*cos(2*Gamma-1.26155) - 0.32123*cos(3*Gamma-1.15710); % Equation of time, min
delta   = (180/pi)*(0.006918-0.399912*cos(Gamma)+0.070257*sin(Gamma)-0.006758*cos(2*Gamma)+0.000907*sin(2*Gamma)-0.002697*cos(3*Gamma)+0.001480*sin(3*Gamma));
end