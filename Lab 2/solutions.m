%part a
%define angles given
theta2 = 40:5:60;
theta4 = [70, 76, 83, 91, 100];

%get k1 k2 k3
A = []
B = []
for i = 1:length(theta2)
    temp1 = [(cosd(theta4(i))) (-1 * (cosd(theta2(i)))) (1)];
    temp2 = cosd(theta2(i)-theta4(i));
    A = [A; temp1];
    B = [B; temp2];
end
%link_ratios = lsqr(A,B); alternative solution
sol = linsolve(A,B);
k1 = sol(1,1)
k2 = sol(2,1)
k3 = sol(3,1)

%get link lengths
d = 180;
a = abs(d/k1);
c = abs(d/k2);
b = abs(sqrt(a * a + c * c + d * d - 2 * a * c * k3));
fprintf("a %d b %d c %d", a,b,c)

%part b
% transmission angle
inputAngle = 40:1:60;
transAngle = acosd((b*b + c*c - (a*a+d*d) + 2 * a * d * cosd(inputAngle))/ (2 * b * c));
plot(inputAngle, transAngle, 'r', 'LineWidth',2)
figure(1)
title('Transmission angle vs Input angle')
xlabel('Input angle (degrees)')
ylabel('Transmission angle (degrees)')
set(gcf,'Position', get(0,'Screensize'));
set(gcf,'Visible','on')
set(gca, 'FontSize', 20)

%part c
%relationship of theta2 to theta4
x = theta2
y = theta4
for i= 1:5
    xy(i) = x(i) * y(i)
end
ySquared = y.^2
xSquared = x.^2
aR = (sum(y)*sum(xSquared) - sum(x)*sum(xy)) / (5 * sum(xSquared) - sum(x).^2)
bR = (5 * sum(xy) - sum(x) * sum(y)) / (5 * sum(xSquared) - sum(x).^2)

theta2C = 40:1:60
theta4C = aR + bR * theta2C
sError = zeros(1,21)
%get structural errors
for j = 1:length(theta2C)
    sError(j) = k1 * cosd(theta4C(j)) - k2 * cosd(theta2C(j)) + k3 - cosd(theta2C(j) - theta4C(j))
end
figure(2)
plot(theta2C, sError, 'r', 'LineWidth',2)
title('Input angle vs Structural Errors')
xlabel('Input angle (degrees)')
ylabel('Structural Errors')
set(gcf,'Position', get(0,'Screensize'));
set(gcf,'Visible','on')
set(gca, 'FontSize', 20)
