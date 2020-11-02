%part a
%three precision points chebysev spacing%
for j = 1:3
    theta2(j) = 0.5*(15 + 165) - 0.5 * (165 - 15) * cosd((180*(2*j - 1))/6);
end

%get theta4%
theta4 = 65 + (0.43*theta2);

%solve link length ratios
syms k1 k2 k3
eq1 = k1 * cosd(theta4(1)) - k2 * cosd(theta2(1)) + k3 == cosd(theta2(1) - theta4(1));
eq2 = k1 * cosd(theta4(2)) - k2 * cosd(theta2(2)) + k3 == cosd(theta2(2) - theta4(2));
eq3 = k1 * cosd(theta4(3)) - k2 * cosd(theta2(3)) + k3 == cosd(theta2(3) - theta4(3));
S=solve(eq1,eq2,eq3,k1,k2,k3);
k1 = double(S.k1);
k2 = double(S.k2);
k3 = double(S.k3);
k1A = k1;
k2A = k2;
k3A = k3;

% find missing link lengths
d = 410;
a = abs(d/k1);
c = abs(d/k2);
b = abs(sqrt(a * a + c * c + d * d - 2 * a * c * k3));
fprintf("a %d b %d c %d", a,b,c)

% transmission angle a
inputAngle = 15:5:165;
transAngleA = acosd((b*b + c*c - (a*a+d*d) + 2 * a * d * cosd(inputAngle))/ (2 * b * c));
plot(inputAngle, transAngleA, 'r', 'LineWidth',2)
title('Transmission angle vs Input angle')
xlabel('Input angle (degrees)')
ylabel('Transmission angle (degrees)')
%set(gcf,'Position', get(0,'Screensize'));
%set(gcf,'Visible','on')
set(gca, 'FontSize', 20)

%part b
%five precision points chebysev spacing%
for j = 1:5
    theta2(j) = 0.5*(15 + 165) - 0.5 * (165 - 15) * cosd((180*(2*j - 1))/10);
end

%get theta4%
theta4 = 65 + (0.43*theta2);

%get k1 k2 k3
A = zeros(3,3)
B = zeros(3,1)
ab = 0
for i = 1:5
    A(1,1) = A(1,1) + cosd(theta4(i)).^2
    A(1,2) = A(1,2) -cosd(theta4(i)) * cosd(theta2(i))
    A(1,3) = A(1,3) + cosd(theta4(i))
    A(2,1) = A(2,1) + cosd(theta2(i)) * cosd(theta4(i))
    A(2,2) = A(2,2) -cosd(theta2(i)).^2
    A(2,3) = A(2,3) + cosd(theta2(i))
    A(3,1) = A(3,1) + cosd(theta4(i))
    A(3,2) = A(3,2) - cosd(theta2(i))
    A(3,3) = i
    B(1,1) = B(1,1) + cos(theta4(i)) * cos(theta2(i) - theta4(i))
    B(2,1) = B(2,1) + cos(theta2(i)) * cos(theta2(i) - theta4(i))
    B(3,1) = B(3,1) + cos(theta2(i) - theta4(i))
end
sol = linsolve(A,B);
k1 = sol(1,1)
k2 = sol(2,1)
k3 = sol(3,1)

%get link lengths
d = 410;
a = abs(d/k1);
c = abs(d/k2);
b = abs(sqrt(a * a + c * c + d * d - 2 * a * c * k3));
fprintf("a %d b %d c %d", a,b,c)

% transmission angle b
inputAngle = 15:5:165;
transAngleB = acosd((b*b + c*c - (a*a+d*d) + 2 * a * d * cosd(inputAngle))/ (2 * b * c));
plot(inputAngle, transAngleB, 'r', 'LineWidth',2)
title('Transmission angle vs Input angle')
xlabel('Input angle (degrees)')
ylabel('Transmission angle (degrees)')
%set(gcf,'Position', get(0,'Screensize'));
%set(gcf,'Visible','on')
set(gca, 'FontSize', 20)

%comparison graph
plotyy(inputAngle, transAngleA, inputAngle, transAngleB)
title('Comparison of structural Error in relation to input angle')
xlabel('Input Angle (degrees)')
ylabel('Transmission Angle (degrees)')
set(gcf,'Position', get(0,'Screensize'));
%set(gcf,'Visible','on')
set(gca, 'FontSize', 20)

%part c
theta2 = 15:5:165
theta4 = 65 + (0.43*theta2);
eA = zeros(1,31)
eB = zeros(1,31)
for j = 1:length(theta2)
    eA(j) = k1A * cosd(theta4(j)) - k2A * cosd(theta2(j)) + k3A - cosd(theta2(j) - theta4(j))
    eB(j) = k1 * cosd(theta4(j)) - k2 * cosd(theta2(j)) + k3 - cosd(theta2(j) - theta4(j))
end
ax = plotyy(theta2,eA,theta2,eB)
title('Comparison of structural Error in relation to input angle')
xlabel('Input Angle (degrees)')
ylabel(ax(1),'Structural Error in A')
ylabel(ax(2),'Structural Error in B')
set(gcf,'Position', get(0,'Screensize'));
set(gcf,'Visible','on')
set(gca, 'FontSize', 20)
