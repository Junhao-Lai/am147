function dydt = YDOT(t, y, m1, m2, L1, L2, g)
    % Extract the state variables
    theta1 = y(1);
    omega1 = y(2);
    theta2 = y(3);
    omega2 = y(4);

    % Compute derivatives
    dydt = zeros(4, 1);
    
    % Equation 1: d(theta1)/dt = omega1
    dydt(1) = omega1;

    % Equation 2: d(theta2)/dt = omega2
    dydt(3) = omega2;

    % Intermediate terms
    d = theta1 - theta2;
    c = cos(d);
    s = sin(d);

    % Equations 3 and 4: Compute angular accelerations
    num1 = -g * (2*m1 + m2) * sin(theta1) - m2 * g * sin(theta1 - 2*theta2);
    num2 = -2 * sin(theta1 - theta2) * m2 * (omega2^2 * L2 + omega1^2 * L1 * cos(theta1 - theta2));
    den1 = L1 * (2 * m1 + m2 - m2 * cos(2 * (theta1 - theta2)));
    dydt(2) = (num1 + num2) / den1;

    num3 = 2 * sin(theta1 - theta2) * (omega1^2 * L1 * (m1 + m2) + g * (m1 + m2) * cos(theta1) + omega2^2 * L2 * m2 * cos(theta1 - theta2));
    den2 = L2 * (2 * m1 + m2 - m2 * cos(2 * (theta1 - theta2)));
    dydt(4) = num3 / den2;
end
