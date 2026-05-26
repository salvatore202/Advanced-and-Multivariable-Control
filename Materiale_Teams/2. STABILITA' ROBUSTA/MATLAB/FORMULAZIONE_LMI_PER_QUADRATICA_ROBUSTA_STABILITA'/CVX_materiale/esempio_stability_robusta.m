clear all
close all
clc

alpha=[0 3.82];
A=[0 1 ; -2 -1];
A1=A+alpha(1)*[0 0 ; -1 0];
A2=A+alpha(2)*[0 0 ; -1 0];



cvx_begin sdp
variable P(2,2) symmetric
A1'*P + P*A1 <= -eye(2)
A2'*P + P*A2 <= -eye(2)
P >= eye(2)
cvx_end