
l1 = -2.13593394;
l2 = -2.13593394;

mu1 = 0.65641833;
mu2 = 0.65641833;

x = (l2 - l1) * cos(mu1 + mu2);
y = mu2 - mu1;

R = 6371000;
d = R * sqrt(x^2 + y^2);

theta = atan2(sin(l2 - l1) * cos(mu2), cos(mu1) * sin(mu2) - sin(mu1) * cos(mu2) * cos(l2 - l1))

N = d * cos(theta)
E = d * sin(theta)

