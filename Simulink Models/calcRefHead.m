function [ref] = calcRefHead(psi, theta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if(acosd(cosd(psi)*cosd(theta) + sind(psi)*sind(theta))<= acosd(-cosd(psi)*cosd(theta) + sind(psi)*sind(theta)))
    alpha1 = acosd(cosd(psi)*cosd(theta) + sind(psi)*sind(theta))
    alpha2 = acosd(-cosd(psi)*cosd(theta) + sind(psi)*sind(theta))
    ref = theta
else
    alpha1 = acosd(cosd(psi)*cosd(theta) + sind(psi)*sind(theta))
    alpha2 = acosd(-cosd(psi)*cosd(theta) + sind(psi)*sind(theta))
    ref = theta - sign(theta)*180
end


