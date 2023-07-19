*M005 Mouse Proteomics: Peptide rollup to protein abundances
*Jack Wiedrick, Adam Burns - 2023.07.07

* Import data
import delimited "lcproteomics_m005_formatted.csv", stringcols(1) clear

* Setup output
tempname lcproteomics_m005_rollup
postfile `lcproteomics_m005_rollup' ///
	protein ///
	sampleid ///
	log2ym ///
	log2ymu ///
	using lcproteomics_m005_rollup.dta, replace

* Transform peptide intensities
generate log2yp=log(cond(y>0,y+1,.))/log(2)
	
*** Fit models and perform rollup
quietly levelsof protein, local(proteins)
foreach protein of local proteins {
	preserve
	keep if protein=="`protein'"
	*Calculate summary statistics 
	quietly summarize log2yp, detail
	*Extract median/50th percentile
	local gm=r(p50)
	if (pepperprot[1]==1) {
		generate mu=log2yp
	}
	else if (pepperprot[1]==2) {
		areg log2yp, absorb(peptide)
		predict dpep, d
		by sampleid (peptide), sort: egen mu=mean(log2yp-dpep)
	}
	else {
		*Mixed effects model with log2yp as the dependent variable, peptideidn as a random effect with one group comprising all observations, and subject id as a fixed effect.
		*Fit model exclusively using EM (emonly)
		mixed log2yp || _all: R.peptide || sampleid:, emonly
		*Calculate best linear unbiased predictions (BLUPs) of the random effects at the level of all observations
		predict dpep, reffects relevel(_all)
		*Calculate fitted values (fixed portion plus contribution of random effects)
		predict mu, fitted
		quietly replace mu=mu-dpep
		by sampleid (mu), sort: replace mu=mu[_n-1] if missing(mu)
	}
	quietly summarize mu, detail
	quietly replace mu=mu-r(p50)+`gm'
	collapse (mean) log2ym=log2yp log2mu=mu, by(sampleid protein)
	capture noisily append using lcproteomics_m005_rollup
	save lcproteomics_m005_rollup, replace
	restore
}

use lcproteomics_m005_rollup, clear
order sampleid protein log2ym log2mu
sort sampleid protein
compress
save "lcproteomics_m005_rollup.dta", replace
