nspins

Comp_I=[1 0; 0 0];
Comp_II=[0 0; 0 1];
Iden_I=eye(4);

L_pulse_90_H1x=kron(Comp_I,kron(expm(-complex(0,1)*(HxNi)*pi/2),conj(expm(-complex(0,1)*(HxNi)*pi/2))))+kron(Comp_II,kron(Iden_I,Iden_I));
L_pulse_90_H2x=kron(Comp_II,kron(expm(-complex(0,1)*(HxNi)*pi/2),conj(expm(-complex(0,1)*(HxNi)*pi/2))))+kron(Comp_I,kron(Iden_I,Iden_I));
L_pulse_90_H12x=L_pulse_90_H1x*L_pulse_90_H2x;

L_pulse_90_H1y=kron(Comp_I,kron(expm(-complex(0,1)*(HyNi)*pi/2),conj(expm(-complex(0,1)*(HyNi)*pi/2))))+kron(Comp_II,kron(Iden_I,Iden_I));
L_pulse_90_H2y=kron(Comp_II,kron(expm(-complex(0,1)*(HyNi)*pi/2),conj(expm(-complex(0,1)*(HyNi)*pi/2))))+kron(Comp_I,kron(Iden_I,Iden_I));
L_pulse_90_H12y=L_pulse_90_H1y*L_pulse_90_H2y;

L_pulse_180_H1x=kron(Comp_I,kron(expm(-complex(0,1)*(HxNi)*pi),conj(expm(-complex(0,1)*(HxNi)*pi))))+kron(Comp_II,kron(Iden_I,Iden_I));
L_pulse_180_H2x=kron(Comp_II,kron(expm(-complex(0,1)*(HxNi)*pi),conj(expm(-complex(0,1)*(HxNi)*pi))))+kron(Comp_I,kron(Iden_I,Iden_I));
L_pulse_180_H12x=L_pulse_180_H1x*L_pulse_180_H2x;

L_pulse_180_H1y=kron(Comp_I,kron(expm(-complex(0,1)*(HyNi)*pi),conj(expm(-complex(0,1)*(HyNi)*pi))))+kron(Comp_II,kron(Iden_I,Iden_I));
L_pulse_180_H2y=kron(Comp_II,kron(expm(-complex(0,1)*(HyNi)*pi),conj(expm(-complex(0,1)*(HyNi)*pi))))+kron(Comp_I,kron(Iden_I,Iden_I));
L_pulse_180_H12y=L_pulse_180_H1y*L_pulse_180_H2y;

L_pulse_90_N1x=kron(Comp_I,kron(expm(-complex(0,1)*(HiNx)*pi/2),conj(expm(-complex(0,1)*(HiNx)*pi/2))))+kron(Comp_II,kron(Iden_I,Iden_I));
L_pulse_90_N2x=kron(Comp_II,kron(expm(-complex(0,1)*(HiNx)*pi/2),conj(expm(-complex(0,1)*(HiNx)*pi/2))))+kron(Comp_I,kron(Iden_I,Iden_I));
L_pulse_90_N12x=L_pulse_90_N1x*L_pulse_90_N2x;

L_pulse_90_N1y=kron(Comp_I,kron(expm(-complex(0,1)*(HiNy)*pi/2),conj(expm(-complex(0,1)*(HiNy)*pi/2))))+kron(Comp_II,kron(Iden_I,Iden_I));
L_pulse_90_N2y=kron(Comp_II,kron(expm(-complex(0,1)*(HiNy)*pi/2),conj(expm(-complex(0,1)*(HiNy)*pi/2))))+kron(Comp_I,kron(Iden_I,Iden_I));
L_pulse_90_N12y=L_pulse_90_N1y*L_pulse_90_N2y;

L_pulse_180_N1x=kron(Comp_I,kron(expm(-complex(0,1)*(HiNx)*pi),conj(expm(-complex(0,1)*(HiNx)*pi))))+kron(Comp_II,kron(Iden_I,Iden_I));
L_pulse_180_N2x=kron(Comp_II,kron(expm(-complex(0,1)*(HiNx)*pi),conj(expm(-complex(0,1)*(HiNx)*pi))))+kron(Comp_I,kron(Iden_I,Iden_I));
L_pulse_180_N12x=L_pulse_180_N1x*L_pulse_180_N2x;

L_pulse_180_N1y=kron(Comp_I,kron(expm(-complex(0,1)*(HiNy)*pi),conj(expm(-complex(0,1)*(HiNy)*pi))))+kron(Comp_II,kron(Iden_I,Iden_I));
L_pulse_180_N2y=kron(Comp_II,kron(expm(-complex(0,1)*(HiNy)*pi),conj(expm(-complex(0,1)*(HiNy)*pi))))+kron(Comp_I,kron(Iden_I,Iden_I));
L_pulse_180_N12y=L_pulse_180_N1y*L_pulse_180_N2y;


