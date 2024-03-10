function [t, y] = pend(t_intv, ICs, t_step, m1, m2, L1, L2, gv, method)
    % Function to simulate the motion of a double pendulum
    
    % Define the ODEs for the double pendulum system
    dydt = @(t,y) YDOT(t, y, m1, m2, L1, L2, gv);
    
    % Choose the integration method
    if strcmpi(method, 'rk')
        [t, y] = rk4(dydt, t_intv, ICs, t_step);
    elseif strcmpi(method, 'trp')
        [t, y] = trapezoid(dydt, t_intv, ICs, t_step);
    elseif strcmpi(method, 'eul')
        [t, y] = euler(dydt, t_intv, ICs, t_step);
    else
        error('Invalid integration method. Choose either ''rk'', ''trp'', or ''eul''.');
    end
end

% Runge-Kutta 4th Order Method
function [t, y] = rk4(dydt, t_intv, y0, h)
    t = linspace(t_intv(1), t_intv(2), h);
    n = length(t);
    y = zeros(n, length(y0));
    y(1, :) = y0;
    for i = 1:n-1
        k1 = dydt(t(i), y(i, :));
        k2 = dydt(t(i) + h/2, y(i, :) + h/2 * k1);
        k3 = dydt(t(i) + h/2, y(i, :) + h/2 * k2);
        k4 = dydt(t(i) + h, y(i, :) + h * k3);
        disp('Size of k1 is:'); disp(size(k1));
        disp('Size of k2 is:'); disp(size(k2));
        disp('Size of k3 is:'); disp(size(k3));
        disp('Size of y(i) is:'); disp(size(y(i, :)));  % Debugging statement
        y(i+1, :) = y(i, :) + h/6 * (k1 + 2*k2 + 2*k3 + k4);

    end
end

% Trapezoidal Method
function [t, y] = trapezoid(dydt, t_intv, y0, h)
    t = linspace(t_intv(1), t_intv(2), h);
    n = length(t);
    y = zeros(n, length(y0));
    y(1, :) = y0;
    for i = 1:n-1
        f = dydt(t(i), y(i, :));
        y(i+1, :) = y(i, :) + h/2 * (f + dydt(t(i+1), y(i, :) + h * f));
    end
end

% Euler Method
function [t, y] = euler(dydt, t_intv, y0, h)
    t = linspace(t_intv(1), t_intv(2), h);
    n = length(t);
    y = zeros(n, length(y0));
    y(1, :) = y0;
    for i = 1:n-1
        y(i+1, :) = y(i, :) + h * dydt(t(i), y(i, :));
    end
end
