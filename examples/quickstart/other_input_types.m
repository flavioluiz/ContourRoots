%% Symbolic and Control System Toolbox inputs (optional toolboxes)
% Run setup_contourroots once per session before this script.

%% 1. Symbolic expressions (Symbolic Math Toolbox)
% Symbolic input is checked automatically: 'AssumeAnalytic' is not needed
% when the expression is built from polynomials, exp, sin, cos, sinh, cosh.
if license('test','Symbolic_Toolbox')
    syms s
    Delta = 1 + s + s^2 + (2*s + 3)*exp(-s);
    [r,info] = croots(Delta, [-8 2 -20 20]);
    fprintf('%d roots, status %s\n', numel(r), info.status);

    p = cpoles(cosh(s)/sinh(s), [-1 1 -10 10])   % i*pi*k, k = -3..3
end

%% 2. tf models with a transport delay (Control System Toolbox)
if license('test','Control_Toolbox')
    s = tf('s');
    G = exp(-2*s)/(s + 1);
    p = cpoles(G, [-3 1 -3 3])          % -1
end
