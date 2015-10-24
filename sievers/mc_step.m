function[value]=mc_step(sigs)
if min(size(sigs))==1,
  value=randn(size(sigs)).*sigs;
else
  assert(issquare(sigs));
  value=create_fake_data(sigs)';
end

