clc; clear; close all;

% ==========================================
% ||          required scripts            ||
% ==========================================

% Units;
Units.initialize();

% ==========================================
% ||         auxiliary functions          ||
% ==========================================

function [pt, en, dir] = transmission_plot(min_energy, max_energy, resolution, direction)
    pt = 0;
    en = linspace(min_energy,max_energy,resolution);
    dir = repmat(direction,1,length(en));
end

function [pt, en, dir] = wavefunction_plot(energies, directions)
    pt = 1;
    en = energies;
    dir = directions;
end


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
bounds = [-1 1-0.375+0.5];

% simulation quaity
% -----------------------------------------
resolution = 1200;


% potential energy
% ---------------------------------------
discrete_potential = true; 
% for descrete potential
potential_dividers = [-0.375 -0.25 -0.125 0 0.125 0.25 0.375 0.5];
potential_values = [0 50 0 50 0 50 0 50 0];

% for continuous potential
potential_function = @(x) 0.5*mass*2*pi*5e10*Units.hz*x^2;


% display
% -----------------------------------------
wavefunction_multiplier = 2;

%[plot_type, energies, directions] =transmission_plot(10, 105, 600, "right");
[plot_type, energies, directions] =wavefunction_plot([11.46], ["right"]);


%%
% ==========================================
% ||            Calculation               ||
% ==========================================

% unit conversion
bounds=bounds                                   *Units.value(distance_units); 
potential_dividers=potential_dividers           *Units.value(distance_units);
potential_values=potential_values               *Units.value(energy_units);
potential_function = @(x) potential_function(x) *Units.value(energy_units);
mass = mass                                     *Units.value(mass_units);
energies = energies                             *Units.value(energy_units);

x = linspace(bounds(1),bounds(2),resolution);


potential = zeros(1,resolution);
dividers = [bounds(1),potential_dividers,bounds(2)];
dx = x(2)-x(1);

current_region = 0;

if discrete_potential
    for i = 1:resolution
        if i == resolution
            break;
        elseif x(i) >= dividers(current_region + 1 )
            current_region = current_region + 1;
        end
        potential(i) = potential_values(current_region);
    
    end
else
    potential = potential(x);
end

% modifying energies such that they represent excess energy above the
% potential
for e = 1:length(energies)
    if strcmp(directions(e), "left") && energies(e)<potential(end)
        fprintf("A particle of energy %.3f %s at index %d doesn't propagate towards the left!\n", abs(energies(e)),energy_units, e);
        energies(e)=[];
        directions(e)=[];

        e = e-1;
    elseif strcmp(directions(e), "right") && energies(e)<potential(1)
        fprintf("A particle of energy %.3f %s at index %d doesn't propagate towards the right!\n", abs(energies(e)),energy_units, e);
        energies(e)=[];
        directions(e)=[];
        e=e-1;
    elseif ~strcmp(directions(e),"left") && ~strcmp(directions(e),"right")
        fprintf("Direction '%s' at index $d is not valid. You can only choose ""left"" or ""right""!",directions(e), e);
        energies(e)=[];
        directions(e)=[];
        e=e-1;
    end
end

% calculating wave numbers
%----------------------------------------------

wavenumbers = zeros(length(energies),resolution);
wave_terms = zeros(length(energies),resolution-1);
P_matrix = zeros(length(energies),2,2);

for e = 1:length(energies)
    wavenumbers(e,:) = sqrt(2*mass*(energies(e)-potential))/Units.hbar;

    wave_terms(e,:) = wavenumbers(e,1:end-1)./wavenumbers(e,2:end);
end

% calculating reflection and transmission
%----------------------------------------------


for e = 1:length(energies)
    P_temp =eye(2);
    for i=1:resolution-1
        % calculates p matrix for every potential region
        P_propagation = [exp(1i*wavenumbers(e,i)*dx) 0 ; 0 exp(-1i*wavenumbers(e,i)*dx)];
        P_step = 0.5*[1+wave_terms(e,i) 1-wave_terms(e,i); 1-wave_terms(e,i) 1+wave_terms(e,i)];
    
        P_temp = P_step * P_propagation * P_temp;
    end
    P_matrix(e,:,:) = P_temp;
end

reflection_coefficients = zeros(1,length(energies));
transmission_coefficients = zeros(1,length(energies));


reflection_coefficients(:) = -P_matrix(:,2,1)./P_matrix(:,2,2);
transmission_coefficients(:) = P_matrix(:,1,1)+reflection_coefficients(:).*P_matrix(:,1,2);

transmittance = abs(transmission_coefficients).^2;
reflectance = abs(reflection_coefficients).^2;


% calculating wave function
if plot_type == 1
    wavefunction_coefficients = zeros (length(energies), 2, resolution);


    for e = 1:length(energies)

        if strcmp(directions(e),"right")
            wavefunction_coefficients(:,1,1) = 1;
            wavefunction_coefficients(:,2,1) = reflection_coefficients(:);
        else
            wavefunction_coefficients(:,1,1) = transmission_coefficients(:);
            wavefunction_coefficients(:,2,1) = 0;
        end
        
        for i=1:resolution
            wavefunction(e,i)= wavefunction_coefficients(e,1,i)*exp(1i*wavenumbers(e,i)*(dx))+ ...
            wavefunction_coefficients(e,2,i)*exp(-1i*wavenumbers(e,i)*dx);

            if i == resolution
                break;
            end
    
            % calculates p matrix for every potential region
            P_propagation = [exp(1i*wavenumbers(e,i)*dx) 0 ; 0 exp(-1i*wavenumbers(e,i)*dx)];
            P_step = 0.5*[1+wave_terms(e,i) 1-wave_terms(e,i); 1-wave_terms(e,i) 1+wave_terms(e,i)];
            

            wavefunction_coefficients(e,:,i+1) = ( (P_step* P_propagation ) * wavefunction_coefficients(e,:,i).' ).';

        end

    end

end



%%
% ==========================================
% ||              Plotting                ||
% ==========================================

% converting quantities back to original units

bounds=bounds                                   /Units.value(distance_units); 
potential = potential                           /Units.value(energy_units);
mass = mass                                     /Units.value(mass_units);
energies = energies                             /Units.value(energy_units);


if plot_type==0
    f=figure; hold on;
    plot(energies,transmittance);
    plot(energies,reflectance);
    line([max(potential) max(potential)], [0 1], Color ='r',LineStyle='--' );

    grid on;
    xlabel(sprintf("Energy: E (%s) ", Units.display(energy_units)))
    xlim([energies(1) energies(end)]);
    title("Reflection and Transmission of Particle incoming at Potential Barrier");
    f.Name = "Reflection/Transmission";

    legend("Transmittance", "Reflectance","Potential peak", Location="east");

elseif plot_type ==1

    figure; hold on;
    plot(x(1:end-1), potential(1:end-1),'k');

    for e = 1:length(energies)
        plot(x, real(wavefunction(e,:))*wavefunction_multiplier+energies(e),'b');
        plot(x, imag(wavefunction(e,:))*wavefunction_multiplier+energies(e), Color='#e677f1');
        plot(x, abs(wavefunction(e,:)).^2*wavefunction_multiplier+energies(e), ':k');
        line([bounds(1) bounds(2)], [energies(e) energies(e)], Color='r', LineStyle='--');
        text(bounds(2)+0.02, energies(e), sprintf('E_{%i}= %.3e %s',e ,energies(e), Units.display(energy_units)), 'Color', 'r', 'FontSize', 10);
    end
    title("Wavefunction scattering states due to potential barrier");
    xlabel(sprintf('Position: x (%s)', Units.display(distance_units)));
    ylabel(sprintf('Energy: E (%s)', Units.display(energy_units)));

    h1 = plot(nan, nan, 'b');                       % real part
    h2 = plot(nan, nan, 'Color', '#e677f1');        % imag part
    h3 = plot(nan, nan, ':k');                     % probability
    h4 = plot(nan, nan, '--r');                     % energy line

    legend([h1 h2 h3 h4], {'Real part','Imaginary part','Probability','Energy'})

end