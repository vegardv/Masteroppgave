close all;
d1 = 0.002886751345948; % I bruk
d2 = 0.003340893189596; % e_cam
d3 = 0.008660254037844;
d4 = 0.004423626322780; % 

x = 320/2;
y = 240/2;

cam1 = [d1*(-x), d1*y;
        d1*(-x), -d1*y;
        d1*(x), -d1*y;
        d1*(x), d1*y;
        d1*(-x), d1*y];
cam2 = [d2*(-x), d2*y;
        d2*(-x), -d2*y;
        d2*(x), -d2*y;
        d2*(x), d2*y;
        d2*(-x), d2*y];
    cam3 = [d3*(-x), d3*y;
        d3*(-x), -d3*y;
        d3*(x), -d3*y;
        d3*(x), d3*y;
        d3*(-x), d3*y];
    
        cam4 = [d4*(-x), d4*y;
        d4*(-x), -d4*y;
        d4*(x), -d4*y;
        d4*(x), d4*y;
        d4*(-x), d4*y];
    
figure(1)
plot(cam1(:,1), cam1(:,2), 'b')
axis equal
hold on
plot(cam2(:,1), cam2(:,2), 'r')
grid on
plot(cam3(:,1), cam3(:,2), 'g')
grid on
plot(cam4(:,1), cam4(:,2), 'm')
grid on
legend('e-CAM51\_USB (DFOV 60^{\circ})', 'e-CAM51\_44x (DFOV 67.5^{\circ})', 'Logitech C910 (DFOV 83^{\circ})', 'Genius WideCam F100 (DFOV 120^{\circ})')
hold off;
xlabel('Length [m]')
ylabel('Width [m]')
%saveas(1, 'C:\Users\vegardvo\Dropbox\Masteroppgave\Report\fig\camFrame.eps', 'eps2c')

