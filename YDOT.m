function dydt = YDOT(t, y, m1, m2, L1, L2, gv)
    % Function to compute the derivatives of the double pendulum system
    
    % Extract state variables
    theta1 = y(1);
    theta2 = y(3);
    theta1_dot = y(2);
    theta2_dot = y(4);
    
    % Equations of motion for the double pendulum system
    % Using Lagrangian Dynamics
    A = (m1 + m2) * L1^2;
    B = m2 * L1 * L2 * cos(theta1 - theta2);
    C = m2 * L1 * L2;
    D = m2 * L2^2;
    E = -m2 * L1 * L2 * theta2_dot^2 * sin(theta1 - theta2) - (m1 + m2) * gv * L1 * sin(theta1);
    F = m2 * L1 * L2 * theta1_dot^2 * sin(theta1 - theta2) - m2 * gv * L2 * sin(theta2);
    
    % Compute accelerations
    theta1_ddot = (E*D - B*F) / (A*D - B*C);
    theta2_ddot = (A*F - C*E) / (A*D - B*C);
    
    % Output derivatives
    dydt = [theta1_dot; theta1_ddot; theta2_dot; theta2_ddot];
end
