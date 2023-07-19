* M005 Mouse Proteomics: Hurdle and linear models to estimate effects of longevity promoting interventions on mouse proteomes.
* Jack Wiedrick, Alicia Feryn, Adam Burns - 2023-07-07

* Load data
use "lcproteomics_m005_normalized.dta", clear

* Create numeric encoding of proteins 
encode protein, generate(protnum)

* Save file with protein names and protein number for downstream merging
keep protein protnum
duplicates drop
save "lcproteomics_M005_protnameskey.dta", replace

* Load data
use "lcproteomics_m005_normalized.dta", clear

* Create numeric encoding of proteins 
encode protein, generate(protnum)


********************************************************************************************
************************** Drug and Diet Treatments ****************************************
********************************************************************************************

* Select tissue (change and rerun for each tissue)
keep if tissue == "liver"

* Do drug and diet treatments first (only HET3 background)
keep if background == "UM-HET3"

* Only keep treatment groups from HET3 background (not including young)
keep if intervention == "control" | intervention == "acarbose" | intervention == "17aE2" | intervention == "canagliflozin" | intervention == "rapamycin" | intervention == "calorierestriction"

* Label treatments
gen treatment = ""
replace treatment="C" if intervention == "control"
replace treatment="Aca" if intervention == "acarbose"
replace treatment="Est" if intervention == "17aE2"
replace treatment="Rapa" if intervention == "rapamycin"
replace treatment="Cana" if intervention == "canagliflozin"
replace treatment="CR" if intervention == "calorierestriction"

* Drop unneeded columns
keep sampleid protein treatment sex log2munorm

* Convert string variables to numeric
destring log2munorm, replace force

* Encode treatment
rename treatment TREATMENT
label define treatment 0 "C" 1 "Aca" 2 "Est" 3 "Rapa" 4 "Cana" 5 "CR"
encode TREATMENT, gen(treatment)
drop TREATMENT

* Encode sex
rename sex SEX
label define sex 0 "F" 1 "M"
encode SEX, gen(sex)
drop SEX

* Tally number of samples/observations per protein
bysort protein: egen n=count(log2munorm)

* Count # of observations per protein per sex and arm group 
* Total samples: 186; Control: 52; Acarbose: 27; Estradiol: 27; Rapamycin: 27; Cana: 27; CR: 26
bysort protein sex treatment: egen npergroup=count(log2munorm)

* Exclude any proteins with groups with less than 2 observations 
bysort protein: egen min_npergroup = min(npergroup)
drop if min_npergroup <2

* Convert all missing to 0
replace log2munorm = 0 if log2munorm == .

* Save formated data
save "lcproteomics_m005_curated.dta", replace


************************** Hurdle model ****************************************

* Reload full data
use "lcproteomics_m005_curated.dta", clear

* Setup file for hurdle model output
tempname lcproteomics_m005_hurdle
postfile `lcproteomics_m005_hurdle' ///
	protnum ///
	converge ///
	b_C_F ///
	se_C_F ///
	b_C_M ///
	se_C_M ///
	b_MvF_C ///
	se_MvF_C ///
	b_CR_F ///
	se_CR_F ///
	b_CR_M ///
	se_CR_M ///
	b_CRvC ///
	se_CRvC ///
	b_CRvC_F ///
	se_CRvC_F ///
	b_CRvC_M ///
	se_CRvC_M ///
	b_CRFvC_nosex ///
	se_CRFvC_nosex ///	
	b_CRMvC_nosex ///
	se_CRMvC_nosex ///
	b_CRFvC_sex ///
	se_CRFvC_sex ///
	b_CRMvC_sex ///
	se_CRMvC_sex ///
	b_CRvC_MvF ///
	se_CRvC_MvF ///
	b_Aca_F ///
	se_Aca_F ///
	b_Aca_M ///
	se_Aca_M ///	
	b_AcavC ///
	se_AcavC ///
	b_AcavC_F ///
	se_AcavC_F ///
	b_AcavC_M ///
	se_AcavC_M ///
	b_AcaFvC_nosex ///
	se_AcaFvC_nosex ///		
	b_AcaMvC_nosex ///
	se_AcaMvC_nosex ///
	b_AcaFvC_sex ///
	se_AcaFvC_sex ///	
	b_AcaMvC_sex ///
	se_AcaMvC_sex ///
	b_AcavC_MvF ///
	se_AcavC_MvF ///
	b_Cana_F ///
	se_Cana_F ///
	b_Cana_M ///
	se_Cana_M ///	
	b_CanavC ///
	se_CanavC ///
	b_CanavC_F ///
	se_CanavC_F ///
	b_CanavC_M ///
	se_CanavC_M ///
	b_CanaFvC_nosex ///
	se_CanaFvC_nosex ///		
	b_CanaMvC_nosex ///
	se_CanaMvC_nosex ///
	b_CanaFvC_sex ///
	se_CanaFvC_sex ///	
	b_CanaMvC_sex ///
	se_CanaMvC_sex ///
	b_CanavC_MvF ///
	se_CanavC_MvF ///
	b_Est_F ///
	se_Est_F ///
	b_Est_M ///
	se_Est_M ///	
	b_EstvC ///
	se_EstvC ///
	b_EstvC_F ///
	se_EstvC_F ///
	b_EstvC_M ///
	se_EstvC_M ///
	b_EstFvC_nosex ///
	se_EstFvC_nosex ///		
	b_EstMvC_nosex ///
	se_EstMvC_nosex ///
	b_EstFvC_sex ///
	se_EstFvC_sex ///	
	b_EstMvC_sex ///
	se_EstMvC_sex ///
	b_EstvC_MvF ///
	se_EstvC_MvF ///
	b_Rapa_F ///
	se_Rapa_F ///
	b_Rapa_M ///
	se_Rapa_M ///	
	b_RapavC ///
	se_RapavC ///
	b_RapavC_F ///
	se_RapavC_F ///
	b_RapavC_M ///
	se_RapavC_M ///
	b_RapaFvC_nosex ///
	se_RapaFvC_nosex ///		
	b_RapaMvC_nosex ///
	se_RapaMvC_nosex ///
	b_RapaFvC_sex ///
	se_RapaFvC_sex ///	
	b_RapaMvC_sex ///
	se_RapaMvC_sex ///
	b_RapavC_MvF ///
	se_RapavC_MvF ///
	using lcproteomics_m005_hurdle.dta, replace

levelsof protnum if n<186, local(prot)
foreach i of local prot {
	di `i'
	
	*Subset data to protein i
	keep if protnum==`i' 
	
	capture noisily churdle exponential log2munorm i.treatment##i.sex, select(i.treatment i.sex) ll(0) vce(oim) iterate(20)
	
	capture noisily churdle exponential log2munorm i.treatment##i.sex if protnum==`i', select(i.treatment i.sex) ll(0) vce(oim) iterate(20)
	if (!_rc & e(converged)) {
		local converge = e(converged)
		
		margins, predict(e(0,.)) at(treatment=(0/5) sex=(0/1))  post
		
		nlcom ///
			(C_F: _b[1._at]) ///
			(C_M: _b[2._at]) ///
			(MvF_C: _b[2._at]-_b[1._at]) ///
			(CR_F: _b[3._at]) ///
			(CR_M: _b[4._at]) ///
			(CRvC: ((_b[3._at]-_b[1._at]) + (_b[4._at]-_b[2._at]))/2) ///
			(CRvC_F: (_b[3._at]-_b[1._at])) ///
			(CRvC_M: (_b[4._at]-_b[2._at])) ///
			(CRFvC_nosex: _b[3._at]-((_b[1._at]+_b[2._at])/2)) ///
			(CRMvC_nosex: _b[4._at]-((_b[1._at]+_b[2._at])/2)) ///
			(CRFvC_sex: ((11*_b[3._at]-_b[5._at]-_b[7._at]-_b[9._at]-_b[11._at]-7*_b[1._at])+(_b[4._at]+_b[6._at]+_b[8._at]+_b[10._at]+_b[12._at]-5*_b[2._at]))/12) ///
			(CRMvC_sex: ((11*_b[4._at]-_b[6._at]-_b[8._at]-_b[10._at]-_b[12._at]-7*_b[2._at])+(_b[3._at]+_b[5._at]+_b[7._at]+_b[9._at]+_b[11._at]-5*_b[1._at]))/12) ///
			(CRvC_MvF: _b[4._at]-_b[2._at]-(_b[3._at]-_b[1._at])) ///
			(Aca_F: _b[5._at]) ///
			(Aca_M: _b[6._at]) ///
			(AcavC: ((_b[5._at]-_b[1._at]) + (_b[6._at]-_b[2._at]))/2) ///
			(AcavC_F: (_b[5._at]-_b[1._at])) ///
			(AcavC_M: (_b[6._at]-_b[2._at])) ///
			(AcaFvC_nosex: _b[5._at]-((_b[1._at]+_b[2._at])/2)) ///
			(AcaMvC_nosex: _b[6._at]-((_b[1._at]+_b[2._at])/2)) ///
			(AcaFvC_sex: ((11*_b[5._at]-_b[3._at]-_b[7._at]-_b[9._at]-_b[11._at]-7*_b[1._at])+(_b[4._at]+_b[6._at]+_b[8._at]+_b[10._at]+_b[12._at]-5*_b[2._at]))/12) ///
			(AcaMvC_sex: ((11*_b[6._at]-_b[4._at]-_b[8._at]-_b[10._at]-_b[12._at]-7*_b[2._at])+(_b[3._at]+_b[5._at]+_b[7._at]+_b[9._at]+_b[11._at]-5*_b[1._at]))/12) ///
			(AcavC_MvF: _b[6._at]-_b[2._at]-(_b[5._at]-_b[1._at])) ///
			(Cana_F: _b[7._at]) ///
			(Cana_M: _b[8._at]) ///
			(CanavC: ((_b[7._at]-_b[1._at]) + (_b[8._at]-_b[2._at]))/2) ///
			(CanavC_F: (_b[7._at]-_b[1._at])) ///
			(CanavC_M: (_b[8._at]-_b[2._at])) ///
			(CanaFvC_nosex: _b[7._at]-((_b[1._at]+_b[2._at])/2)) ///
			(CanaMvC_nosex: _b[8._at]-((_b[1._at]+_b[2._at])/2)) ///
			(CanaFvC_sex: ((11*_b[7._at]-_b[3._at]-_b[5._at]-_b[9._at]-_b[11._at]-7*_b[1._at])+(_b[4._at]+_b[6._at]+_b[8._at]+_b[10._at]+_b[12._at]-5*_b[2._at]))/12) ///
			(CanaMvC_sex: ((11*_b[8._at]-_b[4._at]-_b[6._at]-_b[10._at]-_b[12._at]-7*_b[2._at])+(_b[3._at]+_b[5._at]+_b[7._at]+_b[9._at]+_b[11._at]-5*_b[1._at]))/12) ///
			(CanavC_MvF: _b[8._at]-_b[2._at]-(_b[7._at]-_b[1._at])) ///
			(Est_F: _b[9._at]) ///
			(Est_M: _b[10._at]) ///
			(EstvC: ((_b[9._at]-_b[1._at]) + (_b[10._at]-_b[2._at]))/2) ///
			(EstvC_F: (_b[9._at]-_b[1._at])) ///
			(EstvC_M: (_b[10._at]-_b[2._at])) ///
			(EstFvC_nosex: _b[9._at]-((_b[1._at]+_b[2._at])/2)) ///
			(EstMvC_nosex: _b[10._at]-((_b[1._at]+_b[2._at])/2)) ///
			(EstFvC_sex: ((11*_b[9._at]-_b[3._at]-_b[5._at]-_b[7._at]-_b[11._at]-7*_b[1._at])+(_b[4._at]+_b[6._at]+_b[8._at]+_b[10._at]+_b[12._at]-5*_b[2._at]))/12) ///
			(EstMvC_sex: ((11*_b[10._at]-_b[4._at]-_b[6._at]-_b[8._at]-_b[12._at]-7*_b[2._at])+(_b[3._at]+_b[5._at]+_b[7._at]+_b[9._at]+_b[11._at]-5*_b[1._at]))/12) ///
			(EstvC_MvF: _b[10._at]-_b[2._at]-(_b[9._at]-_b[1._at])) ///
			(Rapa_F: _b[11._at]) ///
			(Rapa_M: _b[12._at]) ///
			(RapavC: ((_b[11._at]-_b[1._at]) + (_b[12._at]-_b[2._at]))/2) ///
			(RapavC_F: (_b[11._at]-_b[1._at])) ///
			(RapavC_M: (_b[12._at]-_b[2._at])) ///
			(RapaFvC_nosex: _b[11._at]-((_b[1._at]+_b[2._at])/2)) ///
			(RapaMvC_nosex: _b[12._at]-((_b[1._at]+_b[2._at])/2)) ///
			(RapaFvC_sex: ((11*_b[11._at]-_b[3._at]-_b[5._at]-_b[7._at]-_b[9._at]-7*_b[1._at])+(_b[4._at]+_b[6._at]+_b[8._at]+_b[10._at]+_b[12._at]-5*_b[2._at]))/12) ///
			(RapaMvC_sex: ((11*_b[12._at]-_b[4._at]-_b[6._at]-_b[8._at]-_b[10._at]-7*_b[2._at])+(_b[3._at]+_b[5._at]+_b[7._at]+_b[9._at]+_b[11._at]-5*_b[1._at]))/12) ///
			(RapavC_MvF: _b[12._at]-_b[2._at]-(_b[11._at]-_b[1._at])), post

		post `lcproteomics_m005_hurdle' ///
			(`i') ///
			(`converge') ///
			(`=_b[C_F]') ///
			(`=_se[C_F]') ///
			(`=_b[C_M]') ///
			(`=_se[C_M]') ///
			(`=_b[MvF_C]') ///
			(`=_se[MvF_C]') ///
			(`=_b[CR_F]') ///
			(`=_se[CR_F]') ///
			(`=_b[CR_M]') ///
			(`=_se[CR_M]') ///
			(`=_b[CRvC]') ///
			(`=_se[CRvC]') ///
			(`=_b[CRvC_F]') ///
			(`=_se[CRvC_F]') ///
			(`=_b[CRvC_M]') ///
			(`=_se[CRvC_M]') ///
			(`=_b[CRFvC_nosex]') ///
			(`=_se[CRFvC_nosex]') ///
			(`=_b[CRMvC_nosex]') ///
			(`=_se[CRMvC_nosex]') ///
			(`=_b[CRFvC_sex]') ///
			(`=_se[CRFvC_sex]') ///
			(`=_b[CRMvC_sex]') ///
			(`=_se[CRMvC_sex]') ///
			(`=_b[CRvC_MvF]') ///
			(`=_se[CRvC_MvF]') ///
			(`=_b[Aca_F]') ///
			(`=_se[Aca_F]') ///
			(`=_b[Aca_M]') ///
			(`=_se[Aca_M]') ///
			(`=_b[AcavC]') ///
			(`=_se[AcavC]') ///
			(`=_b[AcavC_F]') ///
			(`=_se[AcavC_F]') ///
			(`=_b[AcavC_M]') ///
			(`=_se[AcavC_M]') ///
			(`=_b[AcaFvC_nosex]') ///
			(`=_se[AcaFvC_nosex]') ///	
			(`=_b[AcaMvC_nosex]') ///
			(`=_se[AcaMvC_nosex]') ///
			(`=_b[AcaFvC_sex]') ///
			(`=_se[AcaFvC_sex]') ///
			(`=_b[AcaMvC_sex]') ///
			(`=_se[AcaMvC_sex]') ///
			(`=_b[AcavC_MvF]') ///
			(`=_se[AcavC_MvF]') ///
			(`=_b[Cana_F]') ///
			(`=_se[Cana_F]') ///
			(`=_b[Cana_M]') ///
			(`=_se[Cana_M]') ///
			(`=_b[CanavC]') ///
			(`=_se[CanavC]') ///
			(`=_b[CanavC_F]') ///
			(`=_se[CanavC_F]') ///
			(`=_b[CanavC_M]') ///
			(`=_se[CanavC_M]') ///
			(`=_b[CanaFvC_nosex]') ///
			(`=_se[CanaFvC_nosex]') ///	
			(`=_b[CanaMvC_nosex]') ///
			(`=_se[CanaMvC_nosex]') ///
			(`=_b[CanaFvC_sex]') ///
			(`=_se[CanaFvC_sex]') ///
			(`=_b[CanaMvC_sex]') ///
			(`=_se[CanaMvC_sex]') ///
			(`=_b[CanavC_MvF]') ///
			(`=_se[CanavC_MvF]') ///
			(`=_b[Est_F]') ///
			(`=_se[Est_F]') ///
			(`=_b[Est_M]') ///
			(`=_se[Est_M]') ///
			(`=_b[EstvC]') ///
			(`=_se[EstvC]') ///
			(`=_b[EstvC_F]') ///
			(`=_se[EstvC_F]') ///
			(`=_b[EstvC_M]') ///
			(`=_se[EstvC_M]') ///
			(`=_b[EstFvC_nosex]') ///
			(`=_se[EstFvC_nosex]') ///	
			(`=_b[EstMvC_nosex]') ///
			(`=_se[EstMvC_nosex]') ///
			(`=_b[EstFvC_sex]') ///
			(`=_se[EstFvC_sex]') ///
			(`=_b[EstMvC_sex]') ///
			(`=_se[EstMvC_sex]') ///
			(`=_b[EstvC_MvF]') ///
			(`=_se[EstvC_MvF]') ///
			(`=_b[Rapa_F]') ///
			(`=_se[Rapa_F]') ///
			(`=_b[Rapa_M]') ///
			(`=_se[Rapa_M]') ///
			(`=_b[RapavC]') ///
			(`=_se[RapavC]') ///
			(`=_b[RapavC_F]') ///
			(`=_se[RapavC_F]') ///
			(`=_b[RapavC_M]') ///
			(`=_se[RapavC_M]') ///
			(`=_b[RapaFvC_nosex]') ///
			(`=_se[RapaFvC_nosex]') ///	
			(`=_b[RapaMvC_nosex]') ///
			(`=_se[RapaMvC_nosex]') ///
			(`=_b[RapaFvC_sex]') ///
			(`=_se[RapaFvC_sex]') ///
			(`=_b[RapaMvC_sex]') ///
			(`=_se[RapaMvC_sex]') ///
			(`=_b[RapavC_MvF]') ///
			(`=_se[RapavC_MvF]')
	}
		else {
		    post `lcproteomics_m005_hurdle' (`i') (0) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
	}
	
	*Reload full data
	use "lcproteomics_m005_curated.dta", clear
}
postclose `lcproteomics_m005_hurdle'


************************** Linear model ****************************************

* Reload full data
use "lcproteomics_m005_curated.dta", clear

* Setup file for linear model output
tempname lcproteomics_m005_linear
postfile `lcproteomics_m005_linear' ///
	protnum ///
	converge ///
	b_C_F ///
	se_C_F ///
	b_C_M ///
	se_C_M ///
	b_MvF_C ///
	se_MvF_C ///
	b_CR_F ///
	se_CR_F ///
	b_CR_M ///
	se_CR_M ///
	b_CRvC ///
	se_CRvC ///
	b_CRvC_F ///
	se_CRvC_F ///
	b_CRvC_M ///
	se_CRvC_M ///
	b_CRFvC_nosex ///
	se_CRFvC_nosex ///	
	b_CRMvC_nosex ///
	se_CRMvC_nosex ///
	b_CRFvC_sex ///
	se_CRFvC_sex ///
	b_CRMvC_sex ///
	se_CRMvC_sex ///
	b_CRvC_MvF ///
	se_CRvC_MvF ///
	b_Aca_F ///
	se_Aca_F ///
	b_Aca_M ///
	se_Aca_M ///	
	b_AcavC ///
	se_AcavC ///
	b_AcavC_F ///
	se_AcavC_F ///
	b_AcavC_M ///
	se_AcavC_M ///
	b_AcaFvC_nosex ///
	se_AcaFvC_nosex ///		
	b_AcaMvC_nosex ///
	se_AcaMvC_nosex ///
	b_AcaFvC_sex ///
	se_AcaFvC_sex ///	
	b_AcaMvC_sex ///
	se_AcaMvC_sex ///
	b_AcavC_MvF ///
	se_AcavC_MvF ///
	b_Cana_F ///
	se_Cana_F ///
	b_Cana_M ///
	se_Cana_M ///	
	b_CanavC ///
	se_CanavC ///
	b_CanavC_F ///
	se_CanavC_F ///
	b_CanavC_M ///
	se_CanavC_M ///
	b_CanaFvC_nosex ///
	se_CanaFvC_nosex ///		
	b_CanaMvC_nosex ///
	se_CanaMvC_nosex ///
	b_CanaFvC_sex ///
	se_CanaFvC_sex ///	
	b_CanaMvC_sex ///
	se_CanaMvC_sex ///
	b_CanavC_MvF ///
	se_CanavC_MvF ///
	b_Est_F ///
	se_Est_F ///
	b_Est_M ///
	se_Est_M ///	
	b_EstvC ///
	se_EstvC ///
	b_EstvC_F ///
	se_EstvC_F ///
	b_EstvC_M ///
	se_EstvC_M ///
	b_EstFvC_nosex ///
	se_EstFvC_nosex ///		
	b_EstMvC_nosex ///
	se_EstMvC_nosex ///
	b_EstFvC_sex ///
	se_EstFvC_sex ///	
	b_EstMvC_sex ///
	se_EstMvC_sex ///
	b_EstvC_MvF ///
	se_EstvC_MvF ///
	b_Rapa_F ///
	se_Rapa_F ///
	b_Rapa_M ///
	se_Rapa_M ///	
	b_RapavC ///
	se_RapavC ///
	b_RapavC_F ///
	se_RapavC_F ///
	b_RapavC_M ///
	se_RapavC_M ///
	b_RapaFvC_nosex ///
	se_RapaFvC_nosex ///		
	b_RapaMvC_nosex ///
	se_RapaMvC_nosex ///
	b_RapaFvC_sex ///
	se_RapaFvC_sex ///	
	b_RapaMvC_sex ///
	se_RapaMvC_sex ///
	b_RapavC_MvF ///
	se_RapavC_MvF ///
	using lcproteomics_m005_linear.dta, replace

levelsof protnum if n==186, local(prot)
foreach i of local prot {
	di `i'
	
	*Subset data to protein i
	keep if protnum==`i' 
	
	capture noisily regress log2munorm i.treatment##i.sex if protnum==`i', vce(robust)
	if (!_rc & e(converged)) {
		
		local b_CRvC_MvF = `=_b[1.treatment#1.sex]'
		local se_CRvC_MvF = `=_se[1.treatment#1.sex]'
		local b_AcavC_MvF = `=_b[2.treatment#1.sex]'
		local se_AcavC_MvF = `=_se[2.treatment#1.sex]'
		local b_CanavC_MvF = `=_b[3.treatment#1.sex]'
		local se_CanavC_MvF = `=_se[3.treatment#1.sex]'
		local b_EstvC_MvF = `=_b[4.treatment#1.sex]'
		local se_EstvC_MvF = `=_se[4.treatment#1.sex]'
		local b_RapavC_MvF = `=_b[5.treatment#1.sex]'
		local se_RapavC_MvF = `=_se[5.treatment#1.sex]'
		local converge = e(converged)

		margins, at(treatment=(0/5) sex=(0/1))  post
		
		nlcom ///
			(C_F: _b[1._at]) ///
			(C_M: _b[2._at]) ///
			(MvF_C: _b[2._at]-_b[1._at]) ///
			(CR_F: _b[3._at]) ///
			(CR_M: _b[4._at]) ///
			(CRvC: ((_b[3._at]-_b[1._at]) + (_b[4._at]-_b[2._at]))/2) ///
			(CRvC_F: (_b[3._at]-_b[1._at])) ///
			(CRvC_M: (_b[4._at]-_b[2._at])) ///
			(CRFvC_nosex: _b[3._at]-((_b[1._at]+_b[2._at])/2)) ///
			(CRMvC_nosex: _b[4._at]-((_b[1._at]+_b[2._at])/2)) ///
			(CRFvC_sex: ((11*_b[3._at]-_b[5._at]-_b[7._at]-_b[9._at]-_b[11._at]-7*_b[1._at])+(_b[4._at]+_b[6._at]+_b[8._at]+_b[10._at]+_b[12._at]-5*_b[2._at]))/12) ///
			(CRMvC_sex: ((11*_b[4._at]-_b[6._at]-_b[8._at]-_b[10._at]-_b[12._at]-7*_b[2._at])+(_b[3._at]+_b[5._at]+_b[7._at]+_b[9._at]+_b[11._at]-5*_b[1._at]))/12) ///
			(CRvC_MvF: _b[4._at]-_b[2._at]-(_b[3._at]-_b[1._at])) ///
			(Aca_F: _b[5._at]) ///
			(Aca_M: _b[6._at]) ///
			(AcavC: ((_b[5._at]-_b[1._at]) + (_b[6._at]-_b[2._at]))/2) ///
			(AcavC_F: (_b[5._at]-_b[1._at])) ///
			(AcavC_M: (_b[6._at]-_b[2._at])) ///
			(AcaFvC_nosex: _b[5._at]-((_b[1._at]+_b[2._at])/2)) ///
			(AcaMvC_nosex: _b[6._at]-((_b[1._at]+_b[2._at])/2)) ///
			(AcaFvC_sex: ((11*_b[5._at]-_b[3._at]-_b[7._at]-_b[9._at]-_b[11._at]-7*_b[1._at])+(_b[4._at]+_b[6._at]+_b[8._at]+_b[10._at]+_b[12._at]-5*_b[2._at]))/12) ///
			(AcaMvC_sex: ((11*_b[6._at]-_b[4._at]-_b[8._at]-_b[10._at]-_b[12._at]-7*_b[2._at])+(_b[3._at]+_b[5._at]+_b[7._at]+_b[9._at]+_b[11._at]-5*_b[1._at]))/12) ///
			(AcavC_MvF: _b[6._at]-_b[2._at]-(_b[5._at]-_b[1._at])) ///
			(Cana_F: _b[7._at]) ///
			(Cana_M: _b[8._at]) ///
			(CanavC: ((_b[7._at]-_b[1._at]) + (_b[8._at]-_b[2._at]))/2) ///
			(CanavC_F: (_b[7._at]-_b[1._at])) ///
			(CanavC_M: (_b[8._at]-_b[2._at])) ///
			(CanaFvC_nosex: _b[7._at]-((_b[1._at]+_b[2._at])/2)) ///
			(CanaMvC_nosex: _b[8._at]-((_b[1._at]+_b[2._at])/2)) ///
			(CanaFvC_sex: ((11*_b[7._at]-_b[3._at]-_b[5._at]-_b[9._at]-_b[11._at]-7*_b[1._at])+(_b[4._at]+_b[6._at]+_b[8._at]+_b[10._at]+_b[12._at]-5*_b[2._at]))/12) ///
			(CanaMvC_sex: ((11*_b[8._at]-_b[4._at]-_b[6._at]-_b[10._at]-_b[12._at]-7*_b[2._at])+(_b[3._at]+_b[5._at]+_b[7._at]+_b[9._at]+_b[11._at]-5*_b[1._at]))/12) ///
			(CanavC_MvF: _b[8._at]-_b[2._at]-(_b[7._at]-_b[1._at])) ///
			(Est_F: _b[9._at]) ///
			(Est_M: _b[10._at]) ///
			(EstvC: ((_b[9._at]-_b[1._at]) + (_b[10._at]-_b[2._at]))/2) ///
			(EstvC_F: (_b[9._at]-_b[1._at])) ///
			(EstvC_M: (_b[10._at]-_b[2._at])) ///
			(EstFvC_nosex: _b[9._at]-((_b[1._at]+_b[2._at])/2)) ///
			(EstMvC_nosex: _b[10._at]-((_b[1._at]+_b[2._at])/2)) ///
			(EstFvC_sex: ((11*_b[9._at]-_b[3._at]-_b[5._at]-_b[7._at]-_b[11._at]-7*_b[1._at])+(_b[4._at]+_b[6._at]+_b[8._at]+_b[10._at]+_b[12._at]-5*_b[2._at]))/12) ///
			(EstMvC_sex: ((11*_b[10._at]-_b[4._at]-_b[6._at]-_b[8._at]-_b[12._at]-7*_b[2._at])+(_b[3._at]+_b[5._at]+_b[7._at]+_b[9._at]+_b[11._at]-5*_b[1._at]))/12) ///
			(EstvC_MvF: _b[10._at]-_b[2._at]-(_b[9._at]-_b[1._at])) ///
			(Rapa_F: _b[11._at]) ///
			(Rapa_M: _b[12._at]) ///
			(RapavC: ((_b[11._at]-_b[1._at]) + (_b[12._at]-_b[2._at]))/2) ///
			(RapavC_F: (_b[11._at]-_b[1._at])) ///
			(RapavC_M: (_b[12._at]-_b[2._at])) ///
			(RapaFvC_nosex: _b[11._at]-((_b[1._at]+_b[2._at])/2)) ///
			(RapaMvC_nosex: _b[12._at]-((_b[1._at]+_b[2._at])/2)) ///
			(RapaFvC_sex: ((11*_b[11._at]-_b[3._at]-_b[5._at]-_b[7._at]-_b[9._at]-7*_b[1._at])+(_b[4._at]+_b[6._at]+_b[8._at]+_b[10._at]+_b[12._at]-5*_b[2._at]))/12) ///
			(RapaMvC_sex: ((11*_b[12._at]-_b[4._at]-_b[6._at]-_b[8._at]-_b[10._at]-7*_b[2._at])+(_b[3._at]+_b[5._at]+_b[7._at]+_b[9._at]+_b[11._at]-5*_b[1._at]))/12) ///
			(RapavC_MvF: _b[12._at]-_b[2._at]-(_b[11._at]-_b[1._at])), post

		post `lcproteomics_m005_linear' ///
			(`i') ///
			(`converge') ///
			(`=_b[C_F]') ///
			(`=_se[C_F]') ///
			(`=_b[C_M]') ///
			(`=_se[C_M]') ///
			(`=_b[MvF_C]') ///
			(`=_se[MvF_C]') ///
			(`=_b[CR_F]') ///
			(`=_se[CR_F]') ///
			(`=_b[CR_M]') ///
			(`=_se[CR_M]') ///
			(`=_b[CRvC]') ///
			(`=_se[CRvC]') ///
			(`=_b[CRvC_F]') ///
			(`=_se[CRvC_F]') ///
			(`=_b[CRvC_M]') ///
			(`=_se[CRvC_M]') ///
			(`=_b[CRFvC_nosex]') ///
			(`=_se[CRFvC_nosex]') ///
			(`=_b[CRMvC_nosex]') ///
			(`=_se[CRMvC_nosex]') ///
			(`=_b[CRFvC_sex]') ///
			(`=_se[CRFvC_sex]') ///
			(`=_b[CRMvC_sex]') ///
			(`=_se[CRMvC_sex]') ///
			(`=_b[CRvC_MvF]') ///
			(`=_se[CRvC_MvF]') ///
			(`=_b[Aca_F]') ///
			(`=_se[Aca_F]') ///
			(`=_b[Aca_M]') ///
			(`=_se[Aca_M]') ///
			(`=_b[AcavC]') ///
			(`=_se[AcavC]') ///
			(`=_b[AcavC_F]') ///
			(`=_se[AcavC_F]') ///
			(`=_b[AcavC_M]') ///
			(`=_se[AcavC_M]') ///
			(`=_b[AcaFvC_nosex]') ///
			(`=_se[AcaFvC_nosex]') ///	
			(`=_b[AcaMvC_nosex]') ///
			(`=_se[AcaMvC_nosex]') ///
			(`=_b[AcaFvC_sex]') ///
			(`=_se[AcaFvC_sex]') ///
			(`=_b[AcaMvC_sex]') ///
			(`=_se[AcaMvC_sex]') ///
			(`=_b[AcavC_MvF]') ///
			(`=_se[AcavC_MvF]') ///
			(`=_b[Cana_F]') ///
			(`=_se[Cana_F]') ///
			(`=_b[Cana_M]') ///
			(`=_se[Cana_M]') ///
			(`=_b[CanavC]') ///
			(`=_se[CanavC]') ///
			(`=_b[CanavC_F]') ///
			(`=_se[CanavC_F]') ///
			(`=_b[CanavC_M]') ///
			(`=_se[CanavC_M]') ///
			(`=_b[CanaFvC_nosex]') ///
			(`=_se[CanaFvC_nosex]') ///	
			(`=_b[CanaMvC_nosex]') ///
			(`=_se[CanaMvC_nosex]') ///
			(`=_b[CanaFvC_sex]') ///
			(`=_se[CanaFvC_sex]') ///
			(`=_b[CanaMvC_sex]') ///
			(`=_se[CanaMvC_sex]') ///
			(`=_b[CanavC_MvF]') ///
			(`=_se[CanavC_MvF]') ///
			(`=_b[Est_F]') ///
			(`=_se[Est_F]') ///
			(`=_b[Est_M]') ///
			(`=_se[Est_M]') ///
			(`=_b[EstvC]') ///
			(`=_se[EstvC]') ///
			(`=_b[EstvC_F]') ///
			(`=_se[EstvC_F]') ///
			(`=_b[EstvC_M]') ///
			(`=_se[EstvC_M]') ///
			(`=_b[EstFvC_nosex]') ///
			(`=_se[EstFvC_nosex]') ///	
			(`=_b[EstMvC_nosex]') ///
			(`=_se[EstMvC_nosex]') ///
			(`=_b[EstFvC_sex]') ///
			(`=_se[EstFvC_sex]') ///
			(`=_b[EstMvC_sex]') ///
			(`=_se[EstMvC_sex]') ///
			(`=_b[EstvC_MvF]') ///
			(`=_se[EstvC_MvF]') ///
			(`=_b[Rapa_F]') ///
			(`=_se[Rapa_F]') ///
			(`=_b[Rapa_M]') ///
			(`=_se[Rapa_M]') ///
			(`=_b[RapavC]') ///
			(`=_se[RapavC]') ///
			(`=_b[RapavC_F]') ///
			(`=_se[RapavC_F]') ///
			(`=_b[RapavC_M]') ///
			(`=_se[RapavC_M]') ///
			(`=_b[RapaFvC_nosex]') ///
			(`=_se[RapaFvC_nosex]') ///	
			(`=_b[RapaMvC_nosex]') ///
			(`=_se[RapaMvC_nosex]') ///
			(`=_b[RapaFvC_sex]') ///
			(`=_se[RapaFvC_sex]') ///
			(`=_b[RapaMvC_sex]') ///
			(`=_se[RapaMvC_sex]') ///
			(`=_b[RapavC_MvF]') ///
			(`=_se[RapavC_MvF]')
	}
		else {
		    post `lcproteomics_m005_linear' (`i') (0) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
	}
	
	* Reload full data
	use "lcproteomics_m005_curated.dta", clear	
}
postclose `lcproteomics_m005_linear'


************** Combine linear and hurdle results *********

 use "lcproteomics_m005_hurdle.dta", clear

gen model = "hurdle"

append using "lcproteomics_m005_linear.dta"

replace model = "linear" if model == ""

save "lcproteomics_m005_modelres1.dta", replace


******************************* Protein names **********************************

*Merge protein names by protein number 
use "lcproteomics_M005_protnameskey.dta", clear

*Remove duplicate rows
duplicates drop protein, force

*Merge
merge 1:m protnum using "lcproteomics_m005_modelres1.dta", nogenerate keep(match)

*Save
save "lcproteomics_m005_modelres_het3.dta", replace 


*************************************************************************
************************** GHRKO ****************************************
*************************************************************************

* Load data
use "lcproteomics_m005_normalized.dta", clear

* Create numeric encoding of proteins 
encode protein, generate(protnum)

* Select tissue (change and rerun for each tissue)
keep if tissue == "liver"

* Only keep treatment groups from HET3 background (not including young)
keep if intervention == "GHRWT" | intervention == "GHRKO"

* Label treatments
gen treatment = ""
replace treatment="GHRWT" if intervention == "GHRWT"
replace treatment="KO" if intervention == "GHRKO"

* Drop unneeded columns
keep sampleid protein treatment sex log2munorm

* Convert string variables to numeric
destring log2munorm, replace force

* Encode treatment
rename treatment TREATMENT
label define treatment 0 "GHRWT" 1 "KO"
encode TREATMENT, gen(treatment)
drop TREATMENT

* Encode sex
rename sex SEX
label define sex 0 "F" 1 "M"
encode SEX, gen(sex)
drop SEX

* Tally number of samples/observations per protein
bysort protein: egen n=count(log2munorm)

* Count # of observations per protein per sex and arm group 
* Total: 38; GHRWT: 20 samples; GHRKO: 18 
bysort protein sex treatment: egen npergroup=count(log2munorm)

* Exclude any proteins with groups with less than 2 observations 
bysort protein: egen min_npergroup = min(npergroup)
drop if min_npergroup <2

* Convert all missing to 0
replace log2munorm = 0 if log2munorm == .

* Save formated data
save "lcproteomics_m005_curated.dta", replace


************************** Hurdle model ****************************************

* Reload full data
use "lcproteomics_m005_curated.dta", clear

* Setup file for hurdle model output
tempname lcproteomics_m005_hurdle
postfile `lcproteomics_m005_hurdle' ///
	protnum ///
	converge ///
	b_GHRWT_F ///
	se_GHRWT_F ///
	b_GHRWT_M ///
	se_GHRWT_M ///
	b_KO_F ///
	se_KO_F ///
	b_KO_M ///
	se_KO_M ///
	b_KOvGHRWT ///
	se_KOvGHRWT ///
	b_KOvGHRWT_F ///
	se_KOvGHRWT_F ///
	b_KOvGHRWT_M ///
	se_KOvGHRWT_M ///
	b_KOFvGHRWT_nosex ///
	se_KOFvGHRWT_nosex ///
	b_KOMvGHRWT_nosex ///
	se_KOMvGHRWT_nosex ///
	b_KOFvGHRWT_sex ///
	se_KOFvGHRWT_sex ///
	b_KOMvGHRWT_sex ///
	se_KOMvGHRWT_sex ///
	b_MvF_GHRWT ///
	se_MvF_GHRWT ///
	b_KOvGHRWT_MvF ///
	se_KOvGHRWT_MvF ///
	using lcproteomics_m005_hurdle.dta, replace 

levelsof protnum if n<38, local(prot)
foreach i of local prot {
	di `i'
	
	*Subset data to protein i
	keep if protnum==`i' 
	
	capture noisily churdle exponential log2munormalt i.treatment##i.sex if protnum==`i', select(i.treatment i.sex) ll(0) vce(oim) iterate(20)
	*If the command completed successfully (_rc = 0, so !_rc = 1) and the model converged (e(converged) = 0) then...
	if (!_rc & e(converged)) {
		
		local converge = e(converged)
		
		*Estimate margins of all variables on the truncated expected value of y, e(0,.)
		margins, predict(e(0,.)) at(treatment=(0/1) sex=(0/1))  post
		
		nlcom ///
			(GHRWT_F: _b[1._at]) ///
			(GHRWT_M: _b[2._at]) ///
			(KO_F: _b[3._at]) ///
			(KO_M: _b[4._at]) ///
			(KOvGHRWT: ((_b[3._at]-_b[1._at]) + (_b[4._at]-_b[2._at]))/2) ///
			(KOvGHRWT_F: (_b[3._at]-_b[1._at])) ///
			(KOvGHRWT_M: (_b[4._at]-_b[2._at])) ///
			(KOFvGHRWT_nosex: _b[3._at]-((_b[1._at]+_b[2._at])/2)) ///
			(KOMvGHRWT_nosex: _b[4._at]-((_b[1._at]+_b[2._at])/2)) ///
			(KOFvGHRWT_sex: (3*_b[3._at]-3*_b[1._at]+(_b[4._at]-_b[2._at]))/4) ///
			(KOMvGHRWT_sex: (3*_b[4._at]-3*_b[2._at]+(_b[3._at]-_b[1._at]))/4) ///
			(MvF_GHRWT: _b[2._at]-_b[1._at]) ///
			(KOvGHRWT_MvF: _b[4._at]-_b[2._at]-(_b[3._at]-_b[1._at])), post
		
		post `lcproteomics_m005_hurdle' ///
			(`i') ///
			(`converge') ///
			(`=_b[GHRWT_F]') ///
			(`=_se[GHRWT_F]') ///
			(`=_b[GHRWT_M]') ///
			(`=_se[GHRWT_M]') ///
			(`=_b[KO_F]') ///
			(`=_se[KO_F]') ///
			(`=_b[KO_M]') ///
			(`=_se[KO_M]') ///
			(`=_b[KOvGHRWT]') ///
			(`=_se[KOvGHRWT]') ///
			(`=_b[KOvGHRWT_F]') ///
			(`=_se[KOvGHRWT_F]') ///
			(`=_b[KOvGHRWT_M]') ///
			(`=_se[KOvGHRWT_M]') ///
			(`=_b[KOFvGHRWT_nosex]') ///
			(`=_se[KOFvGHRWT_nosex]') ///
			(`=_b[KOMvGHRWT_nosex]') ///
			(`=_se[KOMvGHRWT_nosex]') ///
			(`=_b[KOFvGHRWT_sex]') ///
			(`=_se[KOFvGHRWT_sex]') ///
			(`=_b[KOMvGHRWT_sex]') ///
			(`=_se[KOMvGHRWT_sex]') ///
			(`=_b[MvF_GHRWT]') ///
			(`=_se[MvF_GHRWT]') ///
			(`=_b[KOvGHRWT_MvF]') ///
			(`=_se[KOvGHRWT_MvF]')
}

	else {
		post `lcproteomics_m005_hurdle' (`i') (0) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
	}
	
	*Reload full data
	use "lcproteomics_m005_curated.dta", clear
}
postclose `lcproteomics_m005_hurdle'


************************** Linear model ****************************************

* Reload full data
use "lcproteomics_m005_curated.dta", clear

* Setup file for linear model output
tempname lcproteomics_m005_linear
postfile `lcproteomics_m005_linear' ///
	protnum ///
	converge ///
	b_GHRWT_F ///
	se_GHRWT_F ///
	b_GHRWT_M ///
	se_GHRWT_M ///
	b_KO_F ///
	se_KO_F ///
	b_KO_M ///
	se_KO_M ///
	b_KOvGHRWT ///
	se_KOvGHRWT ///
	b_KOvGHRWT_F ///
	se_KOvGHRWT_F ///
	b_KOvGHRWT_M ///
	se_KOvGHRWT_M ///
	b_KOFvGHRWT_nosex ///
	se_KOFvGHRWT_nosex ///
	b_KOMvGHRWT_nosex ///
	se_KOMvGHRWT_nosex ///
	b_KOFvGHRWT_sex ///
	se_KOFvGHRWT_sex ///
	b_KOMvGHRWT_sex ///
	se_KOMvGHRWT_sex ///
	b_MvF_GHRWT ///
	se_MvF_GHRWT ///
	b_KOvGHRWT_MvF ///
	se_KOvGHRWT_MvF ///
	using lcproteomics_m005_linear.dta , replace

levelsof protnum if n==38, local(prot)
foreach i of local prot {
	di `i'
	
	*Subset data to protein i
	keep if protnum==`i' 
	
	capture noisily regress log2munormalt i.treatment##i.sex if protnum==`i', vce(robust)
	if (!_rc & e(converged)) {
		
		local b_KOvGHRWT_MvF = `=_b[1.treatment#1.sex]'
		local se_KOvGHRWT_MvF = `=_se[1.treatment#1.sex]'
		local converge = e(converged)
		
		margins, at(treatment=(0/1) sex=(0/1))  post
		
		nlcom ///
			(GHRWT_F: _b[1._at]) ///
			(GHRWT_M: _b[2._at]) ///
			(KO_F: _b[3._at]) ///
			(KO_M: _b[4._at]) ///
			(KOvGHRWT: ((_b[3._at]-_b[1._at]) + (_b[4._at]-_b[2._at]))/2) ///
			(KOvGHRWT_F: (_b[3._at]-_b[1._at])) ///
			(KOvGHRWT_M: (_b[4._at]-_b[2._at])) ///
			(KOFvGHRWT_nosex: _b[3._at]-((_b[1._at]+_b[2._at])/2)) ///
			(KOMvGHRWT_nosex: _b[4._at]-((_b[1._at]+_b[2._at])/2)) ///
			(KOFvGHRWT_sex: (3*_b[3._at]-3*_b[1._at]+(_b[4._at]-_b[2._at]))/4) ///
			(KOMvGHRWT_sex: (3*_b[4._at]-3*_b[2._at]+(_b[3._at]-_b[1._at]))/4) ///
			(MvF_GHRWT: _b[2._at]-_b[1._at]) ///
			(KOvGHRWT_MvF: _b[4._at]-_b[2._at]-(_b[3._at]-_b[1._at])), post

		post `lcproteomics_m005_linear' ///
			(`i') ///
			(`converge') ///
			(`=_b[GHRWT_F]') ///
			(`=_se[GHRWT_F]') ///
			(`=_b[GHRWT_M]') ///
			(`=_se[GHRWT_M]') ///
			(`=_b[KO_F]') ///
			(`=_se[KO_F]') ///
			(`=_b[KO_M]') ///
			(`=_se[KO_M]') ///
			(`=_b[KOvGHRWT]') ///
			(`=_se[KOvGHRWT]') ///
			(`=_b[KOvGHRWT_F]') ///
			(`=_se[KOvGHRWT_F]') ///
			(`=_b[KOvGHRWT_M]') ///
			(`=_se[KOvGHRWT_M]') ///
			(`=_b[KOFvGHRWT_nosex]') ///
			(`=_se[KOFvGHRWT_nosex]') ///
			(`=_b[KOMvGHRWT_nosex]') ///
			(`=_se[KOMvGHRWT_nosex]') ///
			(`=_b[KOFvGHRWT_sex]') ///
			(`=_se[KOFvGHRWT_sex]') ///
			(`=_b[KOMvGHRWT_sex]') ///
			(`=_se[KOMvGHRWT_sex]') ///
			(`=_b[MvF_GHRWT]') ///
			(`=_se[MvF_GHRWT]') ///
			(`=_b[KOvGHRWT_MvF]') ///
			(`=_se[KOvGHRWT_MvF]') 
}

	else {
		post `lcproteomics_m005_linear' (`i') (0) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
	}
	
	*Reload full data
	use "lcproteomics_m005_curated.dta", clear
}
postclose `lcproteomics_m005_linear'


************** Combine linear and hurdle results *********

 use "lcproteomics_m005_hurdle.dta", clear

gen model = "hurdle"

append using "lcproteomics_m005_linear.dta"

replace model = "linear" if model == ""

save "lcproteomics_m005_modelres1.dta", replace


******************************* Protein names **********************************

*Merge protein names by protein number 
use "lcproteomics_M005_protnameskey.dta", clear

*Remove duplicate rows
duplicates drop protein, force

*Merge
merge 1:m protnum using "lcproteomics_m005_modelres1.dta", nogenerate keep(match)

*Save
save "lcproteomics_m005_modelres_ghr.dta", replace 



*******************************************************************************
************************** Snell Dwarf ****************************************
*******************************************************************************

* Load data
use "lcproteomics_m005_normalized.dta", clear

* Create numeric encoding of proteins 
encode protein, generate(protnum)

* Select tissue (change and rerun for each tissue)
keep if tissue == "liver"

* Only keep treatment groups from HET3 background (not including young)
keep if intervention == "SnellWT" | intervention == "SnellDW"

* Label treatments
gen treatment = ""
replace treatment="SnellWT" if intervention == "SnellWT"
replace treatment="DW" if intervention == "SnellDW"

* Drop unneeded columns
keep sampleid protein treatment sex log2munorm

* Convert string variables to numeric
destring log2munorm, replace force

* Encode treatment
rename treatment TREATMENT
label define treatment 0 "SnellWT" 1 "DW"
encode TREATMENT, gen(treatment)
drop TREATMENT

* Encode sex
rename sex SEX
label define sex 0 "F" 1 "M"
encode SEX, gen(sex)
drop SEX

* Tally number of samples/observations per protein
bysort protein: egen n=count(log2munorm)

* Count # of observations per protein per sex and arm group 
* Total: 20; SnellWT: 21 samples; SnellDW: 19 
bysort protein sex treatment: egen npergroup=count(log2munorm)

* Exclude any proteins with groups with less than 2 observations 
bysort protein: egen min_npergroup = min(npergroup)
drop if min_npergroup <2

* Convert all missing to 0
replace log2munorm = 0 if log2munorm == .

* Save formated data
save "lcproteomics_m005_curated.dta", replace


************************** Hurdle model ****************************************

* Reload full data
use "lcproteomics_m005_curated.dta", clear

* Setup file for hurdle model output
tempname lcproteomics_m005_hurdle
postfile `lcproteomics_m005_hurdle' ///
	protnum ///
	converge ///
	b_SnellWT_F ///
	se_SnellWT_F ///
	b_SnellWT_M ///
	se_SnellWT_M ///
	b_DW_F ///
	se_DW_F ///
	b_DW_M ///
	se_DW_M ///
	b_DWvSnellWT ///
	se_DWvSnellWT ///
	b_DWvSnellWT_F ///
	se_DWvSnellWT_F ///
	b_DWvSnellWT_M ///
	se_DWvSnellWT_M ///
	b_DWFvSnellWT_nosex ///
	se_DWFvSnellWT_nosex ///
	b_DWMvSnellWT_nosex ///
	se_DWMvSnellWT_nosex ///
	b_DWFvSnellWT_sex ///
	se_DWFvSnellWT_sex ///
	b_DWMvSnellWT_sex ///
	se_DWMvSnellWT_sex ///
	b_MvF_SnellWT ///
	se_MvF_SnellWT ///
	b_DWvSnellWT_MvF ///
	se_DWvSnellWT_MvF ///
	using lcproteomics_m005_hurdle.dta, replace 

levelsof protnum if n<40, local(prot)
foreach i of local prot {
	di `i'
	
	*Subset data to protein i
	keep if protnum==`i' 
	
	capture noisily churdle exponential log2munormalt i.treatment##i.sex if protnum==`i', select(i.treatment i.sex) ll(0) vce(oim) iterate(20)
	*If the command completed successfully (_rc = 0, so !_rc = 1) and the model converged (e(converged) = 0) then...
	if (!_rc & e(converged)) {
		
		local converge = e(converged)
		
		*Estimate margins of all variables on the truncated expected value of y, e(0,.)
		margins, predict(e(0,.)) at(treatment=(0/1) sex=(0/1))  post
		
		nlcom ///
			(SnellWT_F: _b[1._at]) ///
			(SnellWT_M: _b[2._at]) ///
			(DW_F: _b[3._at]) ///
			(DW_M: _b[4._at]) ///
			(DWvSnellWT: ((_b[3._at]-_b[1._at]) + (_b[4._at]-_b[2._at]))/2) ///
			(DWvSnellWT_F: (_b[3._at]-_b[1._at])) ///
			(DWvSnellWT_M: (_b[4._at]-_b[2._at])) ///
			(DWFvSnellWT_nosex: _b[3._at]-((_b[1._at]+_b[2._at])/2)) ///
			(DWMvSnellWT_nosex: _b[4._at]-((_b[1._at]+_b[2._at])/2)) ///
			(DWFvSnellWT_sex: (3*_b[3._at]-3*_b[1._at]+(_b[4._at]-_b[2._at]))/4) ///
			(DWMvSnellWT_sex: (3*_b[4._at]-3*_b[2._at]+(_b[3._at]-_b[1._at]))/4) ///
			(MvF_SnellWT: _b[2._at]-_b[1._at]) ///
			(DWvSnellWT_MvF: _b[4._at]-_b[2._at]-(_b[3._at]-_b[1._at])), post
		
		post `lcproteomics_m005_hurdle' ///
			(`i') ///
			(`converge') ///
			(`=_b[SnellWT_F]') ///
			(`=_se[SnellWT_F]') ///
			(`=_b[SnellWT_M]') ///
			(`=_se[SnellWT_M]') ///
			(`=_b[DW_F]') ///
			(`=_se[DW_F]') ///
			(`=_b[DW_M]') ///
			(`=_se[DW_M]') ///
			(`=_b[DWvSnellWT]') ///
			(`=_se[DWvSnellWT]') ///
			(`=_b[DWvSnellWT_F]') ///
			(`=_se[DWvSnellWT_F]') ///
			(`=_b[DWvSnellWT_M]') ///
			(`=_se[DWvSnellWT_M]') ///
			(`=_b[DWFvSnellWT_nosex]') ///
			(`=_se[DWFvSnellWT_nosex]') ///
			(`=_b[DWMvSnellWT_nosex]') ///
			(`=_se[DWMvSnellWT_nosex]') ///
			(`=_b[DWFvSnellWT_sex]') ///
			(`=_se[DWFvSnellWT_sex]') ///
			(`=_b[DWMvSnellWT_sex]') ///
			(`=_se[DWMvSnellWT_sex]') ///
			(`=_b[MvF_SnellWT]') ///
			(`=_se[MvF_SnellWT]') ///
			(`=_b[DWvSnellWT_MvF]') ///
			(`=_se[DWvSnellWT_MvF]')
}

	else {
		post `lcproteomics_m005_hurdle' (`i') (0) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
	}
	
	*Reload full data
	use "lcproteomics_m005_curated.dta", clear
}
postclose `lcproteomics_m005_hurdle'


************************** Linear model ****************************************

* Reload full data
use "lcproteomics_m005_curated.dta", clear

* Setup file for linear model output
tempname lcproteomics_m005_linear
postfile `lcproteomics_m005_linear' ///
	protnum ///
	converge ///
	b_SnellWT_F ///
	se_SnellWT_F ///
	b_SnellWT_M ///
	se_SnellWT_M ///
	b_DW_F ///
	se_DW_F ///
	b_DW_M ///
	se_DW_M ///
	b_DWvSnellWT ///
	se_DWvSnellWT ///
	b_DWvSnellWT_F ///
	se_DWvSnellWT_F ///
	b_DWvSnellWT_M ///
	se_DWvSnellWT_M ///
	b_DWFvSnellWT_nosex ///
	se_DWFvSnellWT_nosex ///
	b_DWMvSnellWT_nosex ///
	se_DWMvSnellWT_nosex ///
	b_DWFvSnellWT_sex ///
	se_DWFvSnellWT_sex ///
	b_DWMvSnellWT_sex ///
	se_DWMvSnellWT_sex ///
	b_MvF_SnellWT ///
	se_MvF_SnellWT ///
	b_DWvSnellWT_MvF ///
	se_DWvSnellWT_MvF ///
	using lcproteomics_m005_linear.dta , replace

levelsof protnum if n==40, local(prot)
foreach i of local prot {
	di `i'
	
	*Subset data to protein i
	keep if protnum==`i' 
	
	capture noisily regress log2munormalt i.treatment##i.sex if protnum==`i', vce(robust)
	if (!_rc & e(converged)) {
		
		local b_DWvSnellWT_MvF = `=_b[1.treatment#1.sex]'
		local se_DWvSnellWT_MvF = `=_se[1.treatment#1.sex]'
		local converge = e(converged)
		
		margins, at(treatment=(0/1) sex=(0/1))  post
		
		nlcom ///
			(SnellWT_F: _b[1._at]) ///
			(SnellWT_M: _b[2._at]) ///
			(DW_F: _b[3._at]) ///
			(DW_M: _b[4._at]) ///
			(DWvSnellWT: ((_b[3._at]-_b[1._at]) + (_b[4._at]-_b[2._at]))/2) ///
			(DWvSnellWT_F: (_b[3._at]-_b[1._at])) ///
			(DWvSnellWT_M: (_b[4._at]-_b[2._at])) ///
			(DWFvSnellWT_nosex: _b[3._at]-((_b[1._at]+_b[2._at])/2)) ///
			(DWMvSnellWT_nosex: _b[4._at]-((_b[1._at]+_b[2._at])/2)) ///
			(DWFvSnellWT_sex: (3*_b[3._at]-3*_b[1._at]+(_b[4._at]-_b[2._at]))/4) ///
			(DWMvSnellWT_sex: (3*_b[4._at]-3*_b[2._at]+(_b[3._at]-_b[1._at]))/4) ///
			(MvF_SnellWT: _b[2._at]-_b[1._at]) ///
			(DWvSnellWT_MvF: _b[4._at]-_b[2._at]-(_b[3._at]-_b[1._at])), post

		post `lcproteomics_m005_linear' ///
			(`i') ///
			(`converge') ///
			(`=_b[SnellWT_F]') ///
			(`=_se[SnellWT_F]') ///
			(`=_b[SnellWT_M]') ///
			(`=_se[SnellWT_M]') ///
			(`=_b[DW_F]') ///
			(`=_se[DW_F]') ///
			(`=_b[DW_M]') ///
			(`=_se[DW_M]') ///
			(`=_b[DWvSnellWT]') ///
			(`=_se[DWvSnellWT]') ///
			(`=_b[DWvSnellWT_F]') ///
			(`=_se[DWvSnellWT_F]') ///
			(`=_b[DWvSnellWT_M]') ///
			(`=_se[DWvSnellWT_M]') ///
			(`=_b[DWFvSnellWT_nosex]') ///
			(`=_se[DWFvSnellWT_nosex]') ///
			(`=_b[DWMvSnellWT_nosex]') ///
			(`=_se[DWMvSnellWT_nosex]') ///
			(`=_b[DWFvSnellWT_sex]') ///
			(`=_se[DWFvSnellWT_sex]') ///
			(`=_b[DWMvSnellWT_sex]') ///
			(`=_se[DWMvSnellWT_sex]') ///
			(`=_b[MvF_SnellWT]') ///
			(`=_se[MvF_SnellWT]') ///
			(`=_b[DWvSnellWT_MvF]') ///
			(`=_se[DWvSnellWT_MvF]') 
}

	else {
		post `lcproteomics_m005_linear' (`i') (0) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
	}
	
	*Reload full data
	use "lcproteomics_m005_curated.dta", clear
}
postclose `lcproteomics_m005_linear'


************** Combine linear and hurdle results *********

 use "lcproteomics_m005_hurdle.dta", clear

gen model = "hurdle"

append using "lcproteomics_m005_linear.dta"

replace model = "linear" if model == ""

save "lcproteomics_m005_modelres1.dta", replace


******************************* Protein names **********************************

*Merge protein names by protein number 
use "lcproteomics_M005_protnameskey.dta", clear

*Remove duplicate rows
duplicates drop protein, force

*Merge
merge 1:m protnum using "lcproteomics_m005_modelres1.dta", nogenerate keep(match)

*Save
save "lcproteomics_m005_modelres_dw.dta", replace 












