
function probability_unimodal = test_bimodality(sample);

nboot = 500;

[xpdf, n, b] = compute_xpdf(sample);
[dip, p_value, xlow, xup] = HartigansDipSignifTest(xpdf, nboot); 

% figure()
% bar(n, b);
% title(sprintf('Probability of unimodal %.2f', p_value));
probability_unimodal = p_value;
end