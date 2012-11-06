function [a,b,c]=boundries(t,y)
if abs(y(2))-(pi/2)>4;
    a=1;
else
    a=0;
end
b=1;
c=0;