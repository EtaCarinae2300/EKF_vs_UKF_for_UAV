function control_surface = AircraftEOMPulsed(time, pulse_magnitude, control_surface, control_surface_index, pulse_type)
    
if strcmp(pulse_type, 'pulse')
    if time < 1
        del_pulse = pulse_magnitude;
    else
        del_pulse = 0;
    end
    
elseif strcmp(pulse_type, 'doublet')
    if time < 0.5
        del_pulse = pulse_magnitude;
    elseif time < 1  
        del_pulse = -pulse_magnitude;
    else
        del_pulse = 0;
    end
end

    control_surface(control_surface_index) = control_surface(control_surface_index) + del_pulse;