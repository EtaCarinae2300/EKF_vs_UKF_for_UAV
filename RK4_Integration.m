function [x_next] = RK4_Integration(f, t, x, dt)

    k1 = f(t, x);
    k2 = f(t + dt/2, x + k1*dt/2);
    k3 = f(t + dt/2, x + k2*dt/2);
    k4 = f(t + dt, x + k3*dt);
    
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    
    % --- THE FOOLPROOF TRAP ---
    if any(isnan(x_next)) || any(isinf(x_next))
        disp('--- NaN/Inf DETECTED IN RK4 OUTPUT ---');
        disp('The shadow drone state going IN:');
        disp(['Pos (NED): ', num2str(x(1:3)')]);
        disp(['Euler:   ', num2str(x(4:6)')]);
        disp(['Vel:     ', num2str(x(7:9)')]);
        disp(['Rates:   ', num2str(x(10:12)')]);
        
        disp('The corrupted state coming OUT:');
        disp(x_next');
        error('Physics singularity hit during RK4 integration.');
    end

end