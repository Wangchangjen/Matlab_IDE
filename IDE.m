% =========================================================================
% Iterative discrete estimation (IDE)
%   -- inputs:
%       - par: struct of simulation parameters
%       - s: Ux1 complex-valued symbol vector
%       - H: UxB complex-valued channel matrix
%       - noise_var: noise power spectral density (scalar)
%   -- outputs: 
%       - x: Bx1 complex-valued precoded vector
%       - beta: precoding factor (scalar)
%   -- paprmeters:
%       - beta: precoding factor (scalar)
%       - r: penalty parameters (scalar)
%       - alpha: damping factor (scalar)
% -------------------------------------------------------------------------
% (c) 2018 Chang-Jen Wang and Chao-Kai Wen
% e-mail: dkman0988@gmail.com and chaokai.wen@mail.nsysu.edu.tw
% =========================================================================
function [x, beta] = IDE(par,s,H,noise_var)
    % convert to real-valued channel
    U_r = [real(s);imag(s)];
    H_r = [real(H),-1*imag(H);imag(H),real(H)];
    x2 = 0*sqrt(1/(2*par.B))*ones(2*par.B,1);
    x2_old = x2;
    iteration =100;
    r=1;
    r_old = r;
   alpha = 0.95;
    beta=1;
    H_r_b = beta*H_r;
    %% IDE loop
    for t=1:iteration
        b = x2 +  ...
            diag(diag(inv(H_r_b.'*H_r_b+r/2*eye(2*par.B))*H_r_b.'*H_r_b).^-1)*inv(H_r_b.'*H_r_b+r/2*eye(2*par.B))*H_r'*(U_r-H_r_b*(x2));
    % DAC-level
    if par.L == 2
        % 1-bit
        x2= sqrt(1/(2*par.B))*(sign(b));
        x = x2;
    else  
        % multi-bit
        b = b(1:par.B) +sqrt(-1)*b(par.B+1:2*par.B);
        % use QAM-symbol mapping to replace the quantization output
        % cx2= qammod(qamdemod(sqrt(2/3*(P-1))*b/norm(b)*sqrt(N),P),P)/sqrt(2/3*(P-1))/sqrt(N)*norm(b);
        cx2 = par.quantizer(b);
        x2 = [real(cx2);imag(cx2)];
        x = x2;
    end
        r = (trace(H_r_b.'*H_r_b))/(norm(U_r-H_r_b*x2).^2);
    if (mod(t+1,10) == 0) && t < 0.9*iteration;
        x3_t= x2(1:par.B)+sqrt(-1)*x2(par.B+1:2*par.B);
        beta = (real(s'*H*x3_t))/(norm(H*x3_t,2).^2+length(s)*noise_var);
        H_r_b = beta*H_r;
    end
    % damp
        r = (1-alpha)*r + alpha*r_old;
        x2 = (1-alpha)*x2 + alpha*x2_old;
        r_old = r;
        x2_old = x2;
    end
        x = x(1:par.B)+sqrt(-1)*x(par.B+1:2*par.B);

end



