A=diag([-3,-1,-2]);

cvx_begin sdp
variable P(3,3) symmetric
A'*P + P*A <= -eye(3)
P >= eye(3)
cvx_end