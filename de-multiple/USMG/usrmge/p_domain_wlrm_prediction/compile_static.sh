#icc -o p_domain_wlrm_pre_angle_constraints.e   p_domain_wlrm_pre_angle_constraints.cpp -O3 -I/home/swq/fftw/include -L/home/swq/fftw/lib -lfftw3f -lm -g  LS_TAU_P.o -lgfortran

#icc -o local_tau_p.e  -I/home/swq/fftw/include  local_tau_p.cpp  LS_TAU_P.o -lgfortran -L/home/swq/fftw/lib -lfftw3f -lm


gfortran  -c   LS_TAU_P.f -lgfortran -I/home/swq/fftw/include  -L/home/swq/fftw/lib -lfftw3f

