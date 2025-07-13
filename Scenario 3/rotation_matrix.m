function R = rotation_matrix(Omega, i, omega)
    Rz_Omega = [cos(Omega) -sin(Omega) 0;
                sin(Omega)  cos(Omega) 0;
                     0           0     1];
    Rx_i = [1     0           0;
            0 cos(i) -sin(i);
            0 sin(i)  cos(i)];
    Rz_omega = [cos(omega) -sin(omega) 0;
                sin(omega)  cos(omega) 0;
                     0           0     1];
    R = Rz_Omega * Rx_i * Rz_omega;
end

