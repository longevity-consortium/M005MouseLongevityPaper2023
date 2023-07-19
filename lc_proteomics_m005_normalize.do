* M005 Mouse Proteomics: Normalize data using the protein roll-up 
* Jack Wiedrick, Adam Burns - 2023-07-07

* Load and save (in STATA format) sample metadata
import delimited "lcproteomics_M005_metadata.csv", clear

save "lcproteomics_M005_metadata.dta", replace

*Load data
use "lcproteomics_m005_rollup.dta", clear

*Merge with metadata
merge m:1 sampleid using "lcproteomics_M005_metadata", keepusing(tissue sex background intervention cage dissectbatch plate run repeat) assert(match) nogenerate


*Generate tissuegroupsex variable
egen tisgroup=group(tissue intervention sex)

*Normalize
areg log2mu, absorb(tisgroup)
predict log2mur, residuals
egen tissueprotein=group(tissue protein)
areg log2mur, absorb(tissueprotein)
predict log2murr, residuals
sort tissue run sampleid

generate runorder=run
mkspline runorder = runorder, cubic nknots(7) displayknots
mixed log2murr runorder1 runorder2 runorder3 runorder4 runorder5 runorder6 ibn.dissectbatch, noconstant || _all: R.cage || plate: || run:, reml emonly emiterate(300) stddeviations
predict correction, fitted
generate log2munorm=log2mu-correction


*Save results
save "lcproteomics_m005_normalized.dta", replace

