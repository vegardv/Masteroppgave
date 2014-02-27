function [measuredNodePos, measuredNodeHeading, inFrame] = fcn(nodePos, UAVPos, UAVAttitude)
%#codegen


% Calculate area where camera active
	% Init variables
    x_c = 320;
    y_c = 240;
    
    phi = UAVAttitude(1);
    theta = UAVAttitude(2);
    psi = UAVAttitude(3);
    height = UAVPos(3);
    
    d1 = 0;
    lengthPerPixel = 0.001328497723204;
    
    sinPhi = sin(phi);
	cosPhi = cos(phi);

	sinTheta = sin(theta);
	cosTheta = cos(theta);

	sinPsi = sin(psi);
	cosPsi = cos(psi);
    
    corners = zeros(4,2);
    x = 0;
    y = 0;
    for i = 1:4
        switch(i)
            case 1
                x = 0;
                y = 0;
            case 2
                x = 640;
                y = 0;
            case 3
                x = 640;
                y = 480;
            case 4
                x = 0;
                y = 480;
        end;
    
    deltaX = x - x_c;
 	deltaY = y - y_c;  
   
    alpha = -atan2(deltaY, -deltaX);
 	d = sqrt(deltaX^2 + deltaY^2) * lengthPerPixel;
 	beta = -atan(d);
 	sinAlpha = sin(alpha);
 	cosAlpha = cos(alpha);
 	sinBeta = sin(beta);
 	cosBeta = cos(beta);
    
    z1 = sinAlpha*sinBeta*cosPsi*cosTheta - cosAlpha*sinBeta*(-sinPsi*cosPhi + cosPsi*sinTheta*sinPhi) + cosBeta*(sinPsi*sinPhi + cosPsi*sinTheta*cosPhi);
	z2 = sinAlpha*sinBeta*sinPsi*cosTheta - cosAlpha*sinBeta*(cosPsi*cosPhi + sinPsi*sinTheta*sinPhi) + cosBeta * (-cosPsi*sinPhi + sinPsi*sinTheta*sinPhi);
 	z3 = -sinAlpha*sinBeta*sinTheta - cosAlpha*sinBeta*cosTheta*sinPhi + cosBeta*cosTheta*cosPhi;
    
    d2 = (height - d1*cosTheta*cosPhi)/z3;
    
 	corners(i,1) = d2*z1 + d1*(sinPsi*sinPhi + cosPsi*sinTheta*cosPhi);
 	corners(i,2) = d2*z2 + d1*(-cosPsi*sinPhi + sinPsi*sinTheta*cosPhi);    
    end;
 

measuredNodePos = 1; 
measuredNodeHeading = 1;
inFrame = 1;