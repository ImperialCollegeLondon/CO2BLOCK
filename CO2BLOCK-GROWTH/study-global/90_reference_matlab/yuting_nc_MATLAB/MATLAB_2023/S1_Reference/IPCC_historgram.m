load('IPCC_targets_2050.txt')
targets = IPCC_targets_2050(:,1)

Histo = targets;
h = histogram(Histo,30)
xlabel('IPCC compiled storage rate for 2050 [Gt/year]')
ylabel('Frequency')
