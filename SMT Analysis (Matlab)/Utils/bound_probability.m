%Probability that a free molecule would be accounted as bound

function bound_error=bound_probability(Bound_threshold, Db,Frame_interval);

bound_error=(1-exp((-Bound_threshold(1)^2)/(4*Db*Frame_interval)))^Bound_threshold(2);

end

