function dxdt = dinamica_nonlineare(t, x, K, u_bar)
    x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
    
    u = -K * x + u_bar;
    
    dx1 = x2;
    dx2 = -0.8*sin(x1) - 0.5*x2 + 0.3*x4 + u(1);
    dx3 = x4;
    dx4 = -0.6*sin(x3) - 0.4*x4 + 0.2*x2 + u(2) + 0.1*cos(x1);
    
    dxdt = [dx1; dx2; dx3; dx4];
end