function dvecPdt = riccati_ode(t, vecP, A, B, Q, R)
    P = reshape(vecP, 2, 2);
    % Riccati RHS: -Pdot = A'*P + P*A - P*B*R^{-1}*B'*P + Q
    Pdot = -(A'*P + P*A - P*(B*(R\ (B'))) * P + Q);
    dvecPdt = Pdot(:);
end
