function vec = sampleVectorUniformDistribution(wnorm, nsamples)

q = sampleQuaternionUniformDistribution(nsamples);
vec = wnorm.*rotframe([1; 0; 0], q);

end
