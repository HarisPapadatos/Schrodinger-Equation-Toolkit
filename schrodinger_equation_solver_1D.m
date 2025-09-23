clc; clear; close all;

% ==========================================
% ||          required scripts            ||
% ==========================================

% Calculus;
% Units;
% Modes;

Units.initialize();

% ==========================================
% ||               input                  ||
% ==========================================

% units
% ---------------------------------------
energy_units = 'ev';
distance_units = 'nm';
time_units = 's';
mass_units = 'me';

% system
% ---------------------------------------
mass = 1;
periodic = false;
bounds = [-0.25-0.125 0.375+0.125];


% simulation quaity
% -----------------------------------------
resolution = 3000;


% potential energy
% ---------------------------------------

descrete_potential = true; 
% for descrete potential
potential_dividers = [-0.25 -0.125 0 0.125 0.25 0.375];
potential_values = [50 0 50 0 50 0 50];

% for continuous potential
potential_function = @(x) 0.5*mass*2*pi*5e10*Units.hz*x^2;


% display
% -----------------------------------------
display_states = [3];
display_energies = true;
display_wavefunctions = true;
display_wavefunction_multiplier = 5;

mode = Mode.eigenstate;
show_energy_levels = true;

%% Calculation

% calculating units;
bounds=bounds                                   *Units.value(distance_units); 
potential_dividers=potential_dividers           *Units.value(distance_units);
potential_values=potential_values               *Units.value(energy_units);
potential_function = @(x) potential_function(x) *Units.value(energy_units);
mass = mass                                     *Units.value(mass_units);



x = linspace(bounds(1), bounds(2), resolution+2);
dx = x(2)-x(1);
x_interior = x(2:end-1);
box_size = bounds(2)-bounds(1);
deviders = [bounds(1), potential_dividers, bounds(2)];

if length(deviders) ~= size(potential_values)+1
    error("There is not exactly one potential value for every region of space");
end

if bounds(2) <= bounds(1)
    error("The second bound must be to the right");
end

KE_opr = -Units.hbar^2/(2*mass)*Calculus.laplacian1D(resolution, box_size, periodic); 
PE_opr = zeros(resolution, resolution);

current_region = 0;

if descrete_potential
    for i = 1:resolution
        
        if x(i) >= deviders(current_region+1)
            current_region = current_region + 1;
        end
        
        potential(i) = potential_values(current_region);
    
    end
else
    potential = arrayfun(potential_function,x);
end

PE_opr = diag(potential(1:end));

hamiltonian = KE_opr + PE_opr;

[V, D] = eig(hamiltonian);



%% Plotting

% returning back the original values
bounds=bounds                                   /Units.value(distance_units); 
potential = potential                           /Units.value(energy_units);
mass = mass                                     /Units.value(mass_units);

% sizing
eigen_energies = diag(D)';

display_wavefunction_multiplier = display_wavefunction_multiplier*0.04;

max_energy_displayed = max(eigen_energies(display_states));
min_energy_displayed = min(eigen_energies(display_states));
energy_height = max_energy_displayed - min_energy_displayed;

if descrete_potential
    max_potential_displayed = max(potential_values);
    min_potential_displayed = min(potential_values);
    potential_left = potential_values(1);
    potential_right = potential_values(end);
else
    max_potential_displayed = max(potential);
    min_potential_displayed = min(potential);
    potential_left = potential(1);
    potential_right = potential(end);
end

potential_height = max_potential_displayed-min_potential_displayed;

min_point = min([min_energy_displayed min_potential_displayed]);
max_point = max([max_energy_displayed max_potential_displayed]);

amount_displayed = length(display_states);

if display_energies
    height = max_point + potential_height;
else
    height = 1 * max(V(:,display_states));
end


if max_energy_displayed-min_point>(max_potential_displayed-min_point)/1.5
     max_point_displayed =(3*max_energy_displayed+height)/4;
else
    max_point_displayed = (20*max_energy_displayed+max_potential_displayed)/21;
end

screenheight = max_point_displayed-min_point;


V=V./ sqrt(sum(abs(V).^2, 1)/resolution);


figure; hold on;

if mode == Mode.eigenstate || mode == Mode.probability 
    for i = display_states
        if ~periodic
            psi = [0;V(:,i);0];
        else
            psi = V(:,i);
        end

        func = psi*screenheight/sqrt(amount_displayed)*display_wavefunction_multiplier;

        if mode == Mode.probability
            func = psi.^2*screenheight/sqrt(amount_displayed)*display_wavefunction_multiplier;
        elseif i==1
            func = abs(func);
        end

        if display_energies
            energy_line = eigen_energies(i);
            text(bounds(2)+box_size*0.02, energy_line, sprintf('E_{%i}= %.3e %s',i ,eigen_energies(i), Units.display(energy_units)), 'Color', 'r', 'FontSize', 10);
            line([bounds(1) bounds(2)],[1 1]*energy_line,'Color','r','LineStyle','--');
            func = func + energy_line;
        end

        if display_wavefunctions
            if periodic
                plot(x_interior,func,'b');
            else
                plot(x,func,'b');
            end
        end

    end
end

if ~periodic
    line([bounds(1), bounds(1)],[potential_left height],'Color','k');
    line([bounds(2), bounds(2)],[potential_right height],'Color','k');
else
    line([bounds(1), bounds(1)],[potential_left potential_right],'Color','k');
    line([bounds(2), bounds(2)],[potential_right potential_left],'Color','k');
end



plot(x(2:end-1),potential(1:end),'k');


%xlim([bounds(1)-box_size*0.1 bounds(2)+box_size*0.1]);
%ylim([min_point max_point_displayed]);

%xlim([-0.3 0.375+0.05]);
%ylim([-5 55]);

grid on;
title("Wavefunction plots over the potential");
xlabel(sprintf('Position: x (%s)', Units.display(distance_units)));
ylabel(sprintf('Energy: E (%s)', Units.display(energy_units)));
