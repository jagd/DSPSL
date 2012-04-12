%% Capacity / Quasi-TEM Character Resistor Calculation Using MoM

clear

w1 = 2.4e-3;
w2 = 10e-3;
d = 0e-3;

eps_r = 2.2;
% eps_r = 1;
h = 0.79e-3;

port_ext = 2e-3;
mesh_step = min(w1, w2) / 10;


INT_STEP = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mu_0 = 120*pi / 3e8;
eps_0 = 1/(120*pi*3e8);

width_org = max(w1*0.5, w2*0.5 - d) + d + max(w1*0.5 - d, w2*0.5);
port.left = -width_org*0.5 - port_ext;
port.right = width_org*0.5 + port_ext;
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

mesh_num_half = (port.mesh_right - port.mesh_left + 1);
mesh_num = mesh_num_half + mesh_num_half;

%%%%%%%%%%

disp('Filling Matrix');

A = zeros(mesh_num, mesh_num);
eqn_counter = 1;

%%%%%%%%%%
%% Potential on the Conductor

%% Target : Strip1's Conductor
for i = strip1.mesh_left : strip1.mesh_right
	for j = port.mesh_left : port.mesh_right
		a = mesh_step * (abs(j - i) - 0.5);
		b = mesh_step * (abs(j - i) + 0.5);
		X = linspace(a, b, INT_STEP);

		%% upper side
		if i == j
			%% upper side
			A(eqn_counter, j) = mesh_step *(log(0.5*mesh_step) - 1);
		else
			Y = log(X);
			A(eqn_counter, j) = trapz(X, Y);
		end

		%% under side
		Y = log(sqrt(X.*X + h*h));
		A(eqn_counter, mesh_num_half + j) = trapz(X, Y);
	end
	eqn_counter = eqn_counter + 1;
end

%% Target : Strip2's Conductor
for i = strip2.mesh_left : strip2.mesh_right
	for j = port.mesh_left : port.mesh_right
		a = mesh_step * (abs(j - i) - 0.5);
		b = mesh_step * (abs(j - i) + 0.5);
		X = linspace(a, b, INT_STEP);

		%% upper side
		Y = log(sqrt(X.*X + h*h));
		A(eqn_counter, j) = trapz(X, Y);

		%% under side
		if i == j
			%% under side
			A(eqn_counter, mesh_num_half + j) = mesh_step *(log(0.5*mesh_step) - 1);
		else
			Y = log(X);
			A(eqn_counter, mesh_num_half + j) = trapz(X, Y);
		end
	end
	eqn_counter = eqn_counter + 1;
end

A = -A ./ (2*pi*eps_0);

%%%%%%%%%%
%% Electric Field on the Non-Conductor

%% Strip1's Non-Conductor
for i = [port.mesh_left : (strip1.mesh_left - 1)...
	, (strip1.mesh_right + 1) : port.mesh_right]

	%% Only Strip2, because Strip1 do not give
	%% any <y> Component on the Electric Field

	for j = port.mesh_left : port.mesh_right
		a = mesh_step * (abs(j - i) - 0.5);
		b = mesh_step * (abs(j - i) + 0.5);
		X = linspace(a, b, INT_STEP);
		Y = h ./ (X.*X + h*h);
		A(eqn_counter, mesh_num_half + j) = trapz(X, Y);
	end
	A(eqn_counter, i) = +pi*(1+eps_r)/(1 - eps_r); %% fixed?
	eqn_counter = eqn_counter + 1;
end

%% Strip2's Non-Conductor
for i = [port.mesh_left : (strip2.mesh_left - 1)...
	, (strip2.mesh_right + 1) : port.mesh_right]

	%% Only Strip1, because Strip2 do not give
	%% any <y> Component on the Electric Field

	for j = port.mesh_left : port.mesh_right
		a = mesh_step * (abs(j - i) - 0.5);
		b = mesh_step * (abs(j - i) + 0.5);
		X = linspace(a, b, INT_STEP);
		Y = h ./ (X.*X + h*h);
		A(eqn_counter, j) = -trapz(X, Y);
	end
	A(eqn_counter, mesh_num_half + i) = -pi*(1+eps_r)/(1 - eps_r); %% fixed?
	eqn_counter = eqn_counter + 1;
end


%%%%%%%%%%
%% Free Charge

%% K is the factor := (free Charge) / (all Charge)

K = zeros(mesh_num_half, mesh_num_half);
k_aux = (1 - eps_r)/2/pi;

for i = port.mesh_left : port.mesh_right
	%% only conside the other interface
	for j = port.mesh_left : port.mesh_right
		a = mesh_step * (abs(j - i) - 0.5);
		b = mesh_step * (abs(j - i) + 0.5);
		X = linspace(a, b, INT_STEP);
		Y = h ./ (X.*X + h*h);
		K(i, mesh_num_half + j) = k_aux*trapz(X, Y);
		K(mesh_num_half + i, j) = K(i, mesh_num_half + j);
		if i == j
			K(i, i) = (eps_r + 1)*0.5;
			K(mesh_num_half + i, mesh_num_half + i) = K(i, i);
		end
	end
end

clear k_aux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% equivalent Charge

disp('Inversing Matrix');

	%%  k1 * Potential1 + k2 *Potential2 = 0
	A1 = (A^-1);
	B = K*A1;

	tmp = (strip1.mesh_right - strip1.mesh_left + 1);
	i = [strip1.mesh_left : strip1.mesh_right];
	k1 = sum(sum(B(i, 1:tmp)));

	i = mesh_num_half + [strip2.mesh_left : strip2.mesh_right];
	k2 = sum(sum(B(i, (tmp + 1) : ...
		(tmp + strip2.mesh_right - strip2.mesh_left + 1))));

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

	tmp = strip1.mesh_right - strip1.mesh_left + 1;
	potential = zeros(mesh_num, 1);
	potential(1 : tmp, 1) = p1;
	potential((tmp + 1) : (tmp + strip2.mesh_right - strip2.mesh_left + 1) , 1) = p2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

charge_all = A1 * potential * mesh_step;
charge_free = K * charge_all;

Q(1) = sum( charge_free(strip1.mesh_left : strip1.mesh_right) );
Q(2) = sum( charge_free(mesh_num_half + [strip2.mesh_left : strip2.mesh_right]) );

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
bar(charge_free(1:mesh_num_half));

subplot(212);
bar(charge_free((mesh_num_half + 1):mesh_num));

%%%%%%%%%%

figure;
title('all charge');

subplot(211);
bar(charge_all(1:mesh_num_half));

subplot(212);
bar(charge_all((mesh_num_half + 1):mesh_num));
