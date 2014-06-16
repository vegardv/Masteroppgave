function pathPlotterNode(x, y,  cornersX, cornersY, tsamp, dec, inframe, nodePos)
hold on;

tnow=0;
cornersX = [cornersX, cornersX(:,1)];
cornersY = [cornersY, cornersY(:,1)];
for now=1:dec:length(x)
    if inframe(now)
       h1 = plot(cornersY(now, :),cornersX(now, :),'g');
    else
        h2 =plot(cornersY(now, :),cornersX(now, :),'r');
    end

    tnow=tnow+tsamp*dec;
end
h3 = plot(nodePos(:,2), nodePos(:,1), 'xk', 'MarkerSize', 15, 'LineWidth', 2);
%h3 = plot(nodePos(:,2), nodePos(:,1), 'k')%, 'MarkerSize', 15, 'LineWidth', 2);

h = plot(y, x)
h4 = scatter(0,0, 'o', 'm', 'LineWidth', 3);
legend([h4,h,h2,h1,h3], 'Start position of UAV', 'Position of UAV', 'Camera frame does not contain node','Camera frame containing node', 'Node position')
hold off
xlabel('East [m]')
ylabel('North [m]')
axis equal
grid on

hold off