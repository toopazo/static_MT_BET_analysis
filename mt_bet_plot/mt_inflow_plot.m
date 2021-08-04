
omega = 2000 * pi / 30;         % 2000 RPM
R = 0.2;
A = pi * R^2;
Vtip = omega * R;

Thover = 5 * 9.8;               % vehicle mass kg
rho = 1.225;                    % air density kg/m3
CT = Thover / ( rho*A*Vtip^2 );
mu = 0;
alpha = pi/2;

Vmax = 50;                      % vertical speed m/s
lambda_c_final = Vmax / Vtip;
lambda_arr = [];
lambda_c_arr = 0:0.01:lambda_c_final;
for lambda_c = lambda_c_arr
    lambda_arr = [lambda_arr; mt_inflow(CT, mu, alpha, Vtip, lambda_c)];
end
lambda_h = sqrt(CT/2)
%lambda_h = lambda_arr(1)
%lambda_h_err = lambda_arr(1) - sqrt(CT/2)
%return

fig = figure(1);
hold on;
plot(lambda_c_arr ./ lambda_h, lambda_arr ./ lambda_h);
plot(lambda_c_arr ./ lambda_h, lambda_c_arr ./ lambda_h, '-.');
title('MT predicted axial climb inflow')
xlabel('\lambda_c / \lambda_h');
ylabel('\lambda / \lambda_h');
textstr = ['CT ' num2str(CT) ' \lambda_h ' num2str(lambda_h)];
text(0, 3, textstr)
grid on;
axis equal
