function esp = ExpGrow_1Cmp_fun(coeff, x)


k = coeff(1);
a = coeff(2);


esp = a*( 1 - (exp(-k*x)) );