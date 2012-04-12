%% Capacity / Quasi-TEM Character Resistor Calculation Using MoM

clear

w1 = 2e-3;
w2 = 10e-3;
d = 0e-3;

eps_r = 1;
h = 0.79e-3;

mesh_step = min(w1, w2) / 10;


INT_STEP = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c0 = 3e8;
z0 = 120*pi;
mu_0 =  z0 / c0;
eps_0 = 1/(z0 * c0);

width_org = max(w1*0.5, w2*0.5 - d) + d + max(w1*0.5 - d, w2*0.5);
port.left = -width_org*0.5;
port.right = width_org*0.5;
port.mesh_left = 1;
port.mesh_right = round((port.right - port.left) / mesh_step);

strip1.left = width_org*0.5 - max(w1*0.5, w2*0.5 - d) - w1*0.5;
strip1.right = width_org*0.5 - max(w1*0.5, w2*0.5 - d) + w1*0.5;
strip1.mesh_left = 1 + round((strip1.left - port.left) / mesh_step);
strip1.mesh_right = round((strip1.right - port.left) / mesh_step);

strip2.left = -0.5*width_org + max(w1*0.5 - d, w2*0.5) - w2*0.5;
strip2.right = -0.5*width_org + max(w1*0.5 - d, w2*0.5) + w2*0.5;
strip2.mesh_left = 1 + round((strip2.left - port.left) / mesh_step);
strip2.mesh_right = round((strip2.right - port.left) / mesh_step);

mesh_1 = (strip1.mesh_right - strip1.mesh_left + 1);
mesh_2 = (strip2.mesh_right - strip2.mesh_left + 1);

mesh_num = mesh_1 + mesh_2;

%%%%%%%%%%

disp('Filling Matrix');

A = zeros(mesh_num, mesh_num);
eqn_counter = 1;

%%%%%%%%%%
%% Potential on the Conductor

%% Target : Strip1's Conductor
for i = strip1.mesh_left : strip1.mesh_right
	for j = strip1.mesh_left : strip1.mesh_right
		%% upper side
		if i == j
			%% upper side
			A(eqn_counter, (j - strip1.mesh_left + 1)) = mesh_step *(log(0.5*mesh_step) - 1);
		else
			a = mesh_step * (abs(j - i) - 0.5);
			b = mesh_step * (abs(j - i) + 0.5);
			X = linspace(a, b, INT_STEP);
			Y = log(X);
			A(eqn_counter, (j - strip1.mesh_left + 1)) = trapz(X, Y);
		end
	end

	for j = strip2.mesh_left : strip2.mesh_right
		a = mesh_step * (abs(j - i) - 0.5);
		b = mesh_step * (abs(j - i) + 0.5);
		X = linspace(a, b, INT_STEP);

		%% under side
		Y = log(sqrt(X.*X + h*h));
		A(eqn_counter, mesh_1 + j - strip2.mesh_left + 1) = trapz(X, Y);
	end
	eqn_counter = eqn_counter + 1;
end

%% Target : Strip2's Conductor
for i = strip2.mesh_left : strip2.mesh_right
	for j = strip2.mesh_left : strip2.mesh_right
		%% under side
		if i == j
			%% under side
			A(eqn_counter, mesh_1 + j - strip2.mesh_left + 1) = mesh_step *(log(0.5*mesh_step) - 1);
		else
			a = mesh_step * (abs(j - i) - 0.5);
			b = mesh_step * (abs(j - i) + 0.5);
			X = linspace(a, b, INT_STEP);
			Y = log(X);
			A(eqn_counter, mesh_1 + j - strip2.mesh_left + 1) = trapz(X, Y);
		end
	end

	for j = strip1.mesh_left : strip1.mesh_right
		a = mesh_step * (abs(j - i) - 0.5);
		b = mesh_step * (abs(j - i) + 0.5);
		X = linspace(a, b, INT_STEP);

		%% upper side
		Y = log(sqrt(X.*X + h*h));
		A(eqn_counter, j - strip1.mesh_left + 1) = trapz(X, Y);

	end
	eqn_counter = eqn_counter + 1;
end

A = -A ./ (2*pi*eps_0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% equivalent Charge

disp('Inversing Matrix');

	%%  k1 * Potential1 + k2 *Potential2 = 0
	A1 = (A^-1);
	B = A1;

	i = [1 : strip1.mesh_right - strip1.mesh_left + 1];
	k1 = sum(sum(B(i, 1:mesh_1)));

	i = mesh_1 + [1 : strip2.mesh_right - strip2.mesh_left + 1];
	k2 = sum(sum(B(i, (mesh_1 + 1) : mesh_num)));

	if k1 == 0
		p2 = 0;
		p1 = 1;
	elseif k2 == 0
		p1 = 0;
		p2 = 1;
	else
		p1 = 1;
		p2 = -k1/k2;
	end

	potential = p2*ones(mesh_num, 1);
	potential(1 : mesh_1, 1) = p1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

charge_free = A1 * potential * mesh_step;

Q(1) = sum( charge_free(1 : mesh_1));
Q(2) = sum( charge_free((mesh_1 + 1) : mesh_num));

C = (Q(1) - Q(2))*0.5 / (p1 - p2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output

mesh_num
Q
percent_Q = sum(Q)/Q(1)
C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

disp('Plotting');

figure;
title('free charge');

subplot(211);
bar([zeros(1, strip1.mesh_left) , charge_free(1:mesh_1)' , zeros(1, port.mesh_right - strip1.mesh_right)]);

subplot(212);
bar([zeros(1, strip2.mesh_left) , charge_free((mesh_1 + 1):mesh_num)' , zeros(1, port.mesh_right - strip2.mesh_right)]);
