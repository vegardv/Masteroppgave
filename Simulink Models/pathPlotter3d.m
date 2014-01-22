function pathPlotter3d(x, y,z,  psi, tsamp, dec)
%PATHPLOTTER draws the path of the ship MS Fartoystyring used in TTK4190 
%Guidance and Control Assignment 3, Tasks 2.3 to 2.7, inclusive.
%
%PATHPLOTTER outputs a single xy-plot of the MS Fartoystyring's trajectory 
%and, depending upon the input, either waypoints and path or the trajectory
%of the target. The MS Fartoystyring is yellow and the target green. MS
%Fartoystyring's trajectory is blue, and the waypoints/target trajectory
%red.
%
%x      is the ship's north position (in NED). x is a vector in R^T, T is 
%          number of samples
%y      is the ship's east position (in NED). y is a vector in R^T 
%psi    is the ship's yaw angle (in NED). psi is a vector in R^T 
%tsamp  is the sampling time
%dec>=1 is how much of the data should be used in the plot. dec should 
%          be a natural number. E.g. dec=2 reduces number of data points by
%          a factor of 2. Too low value of dec makes the plot illegible
%tstart is simulation start time
%tstop  is simulation stop time
%track  is a boolean value indicating whether or not the ship is following 
%          waypoints (Tasks 2.3-2.6) (track=0) or tracking a target (Task
%          2.7) (track=1).
%
%Input must be sampled with a fixed time step.
%
%Author   : Christian Holden
%Date     : 2008-02-22
%Revisions: 2009-01-15 C. Holden. Updated for new semester. 
%           2010-02-22 C. Holden. Updated for new semester.
%           2010-03-02 C. Holden. Minor bug fix
%           2011-03-11 A. Lekkas  Updated for new semester
%           2012-03-03 A. Lekkas  Updated for new semester
%           2012-03-15 A. Lekkas  Updated for new semester
%You are free to modify the code as necessary.
%
%Bugs should be reported to the TA.
close all

hold on;

plot3(y, x, z)
% L=0.2;
% 
% tnow=0;
% for now=1:dec:length(x)
% 
%     b = L/8;
%     %MS Fartoystyring
%     tmpR=[cos(psi(now)) -sin(psi(now)); sin(psi(now)) cos(psi(now))];
%     boat = tmpR*[L/2, .9*L/2, .5*L/2, b/2, b/2, - b/2, - b/2, -L/2, -L/2,- b/2, - b/2, b/2, b/2, .5*L/2, .9*L/2, L/2; 
%               %0 10 20 20 -20 -20 -10 0]
%               0, b/4, b/2, b/2, L/2, L/2, b/2,  b/2, -b/2, -b/2, -L/2, -L/2, -b/2, -b/2, -b/4, 0];
%     plot(y(now)+boat(2,:),x(now)+boat(1,:),'b');
%     patch(y(now)+boat(2,:),x(now)+boat(1,:),'y');
%     
%     north = tmpR*[L/2, .9*L/2, .5*L/2, b/2, b/2, .5*L/2, .9*L/2, L/2; 
%               %0 10 20 20 -20 -20 -10 0]
%                  0, b/4, b/2, b/2, -b/2, -b/2, -b/4, 0];
%     patch(y(now)+north(2,:),x(now)+north(1,:),'r');
%     tnow=tnow+tsamp*dec;
% end
% hold off
% xlabel('East [m]')
% ylabel('North [m]')
% axis equal
% grid on
hold off