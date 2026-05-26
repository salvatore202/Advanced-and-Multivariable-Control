% Integrate xdot = (A-B*K(t)) x  from 0 to T
function xdot = closed_loop_ode(t, x, A, B, K_of_t)
    Kt = K_of_t(t)';  % 1x2
    xdot = (A - B*Kt)*x;
end