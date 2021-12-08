

program simk
	replace x = rnormal()
	ksmirnovk x, by(y)
end

simulate pvaluesKS=e(pKS) pvaluesCVM=e(pCVM), reps(2000): simk
