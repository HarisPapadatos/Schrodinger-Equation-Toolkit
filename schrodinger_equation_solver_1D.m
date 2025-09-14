clc; close all;

% functions & constants

function L = laplacian1D(n, width)
    h = width/(n+1);
    e = ones(n,1);
    L = spdiags([e -2*e e], -1:1, n, n) / h^2;
end

SI_units = {1.0545e-34 9.109e-31 'J'};
eV_units = {6.582e-16 5.109e5 'eV'};
dimensionless = {1 0.5 ''};

% input

%system
[hbar, m0, energy_units] = dimensionless{:};
resolution = 400;
bounds = [-50 50]; % width = 1e-10 m 

x = linspace(bounds(1), bounds(2), resolution+2);

descrete_potential = true; 
% for descrete potential
potential_deviders = [10 20 30];
potential_values = [0 0.2 0 0.2];
% for continuous potential
potential = x.^2/1000;


% display
display_states = [4:8];
display_energies = true;
display_wavefunctions = true;
display_wavefunction_multiplier = 0.05;

mode = Mode.eigenstate;

%% Setting up
dx = x(2)-x(1);
x_interior = x(2:end-1);
box_size = bounds(2)-bounds(1);
deviders = [bounds(1), potential_deviders, bounds(2)];

if length(deviders) ~= size(potential_values)+1
    error("There is not exactly one potential value for every region of space");
end

if bounds(2) <= bounds(1)
    error("The second bound must be to the right");
end

KE_opr = -hbar^2/(2*m0)*laplacian1D(resolution, box_size); 
PE_opr = zeros(resolution, resolution);

current_region = 0;

if descrete_potential
    for i = 1:resolution
        
        if x(i) >= deviders(current_region+1)
            current_region = current_region + 1;
        end
        
        PE_opr(i, i) = potential_values(current_region);
    
    end
else
    PE_opr = diag(potential(2:(end-1)));
end

hamiltonian = KE_opr + PE_opr;

[V, D] = eig(hamiltonian);



%% Plotting
clc; close all;

eigen_energies = diag(D)';

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
        psi = [0;V(:,i);0];

        func = psi*screenheight/sqrt(amount_displayed)*display_wavefunction_multiplier;

        if mode == Mode.probability
            func = psi.^2*screenheight/sqrt(amount_displayed)*display_wavefunction_multiplier;
        elseif i==1
            func = abs(func);
        end

        if display_energies
            energy_line = eigen_energies(i);
            text(bounds(2)+box_size*0.02, energy_line, sprintf('E_{%i}= %.3e %s',i ,eigen_energies(i), energy_units), 'Color', 'r', 'FontSize', 10);
            line([bounds(1) bounds(2)],[1 1]*energy_line,'Color','r','LineStyle','--');
            func = func + energy_line;
        end

        if display_wavefunctions
                plot(x,func,'b');
        end

    end
end

line([bounds(1), bounds(1)],[potential_left height],'Color','k');
line([bounds(2), bounds(2)],[potential_right height],'Color','k');

if descrete_potential
    for i = 1:(length(deviders)-1)
        line([deviders(i) deviders(i+1)],[1 1]*potential_values(i),'Color','k');
    
        if i ~=length(deviders)-1
            line([deviders(i+1) deviders(i+1)],[potential_values(i) potential_values(i+1)],'Color','k');
        end
    end
else
    plot(x,potential,'k');
end


xlim([bounds(1)-box_size*0.1 bounds(2)+box_size*0.1]);
ylim([min_point max_point_displayed]);

