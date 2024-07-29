* Hypercapnia TriNetX Computable Phenotype Analysis
* Updated 4/14/2023 BWL

/* ---------
OVERALL TODO LIST


-----------*/

capture log close

* Load data
clear

//Replace with your working directory**
//cd "C:\Users\reblo\Box\Residency Personal Files\Scholarly Work\Locke Research Projects\TriNetX Code" 
cd "/Users/blocke/Box Sync/Residency Personal Files/Scholarly Work/Locke Research Projects/TriNetX Code"


program define datetime 
end

capture mkdir "Results and Figures"
capture mkdir "Results and Figures/$S_DATE/" //make new folder for figure output if needed
capture mkdir "Results and Figures/$S_DATE/Logs/" //new folder for stata logs
local a1=substr(c(current_time),1,2)
local a2=substr(c(current_time),4,2)
local a3=substr(c(current_time),7,2)
local b = "Hypercapnia Case Defintions.do" // do file name
copy "`b'" "Results and Figures/$S_DATE/Logs/(`a1'_`a2'_`a3')`b'"

set scheme cleanplots
graph set window fontface "Helvetica"
log using temp.log, replace

clear
cd "Data"
use full_db
//use subsample_db_5_perc
cd ..

/* -----------------------

Define program to write Se, Sp, +LR, PPV to table
CSV file
Takes a list of ICD codes 

-----------------------*/ 
capture program drop test_char_from_icd_list
program define test_char_from_icd_list
    syntax, ref_std(string) icdcodes(string)

    * Loop through each ICD code in the list provided in the syntax
    foreach icdcode of local icdcodes {
		quietly diagt `ref_std' `icdcode'
		* Store the results in local macros with two decimal formatting
		local sens = string(r(sens), "%9.2f")
		local sens_lb = string(r(sens_lb), "%9.2f")
		local sens_ub = string(r(sens_ub), "%9.2f")
		local spec = string(r(spec), "%9.2f")
		local spec_lb = string(r(spec_lb), "%9.2f")
		local spec_ub = string(r(spec_ub), "%9.2f")
		local ppv = string(r(ppv), "%9.2f")
		local ppv_lb = string(r(ppv_lb), "%9.2f")
		local ppv_ub = string(r(ppv_ub), "%9.2f")
		local lrp = string(r(lrp), "%9.2f")
		local lrp_lb = string(r(lrp_lb), "%9.2f")
		local lrp_ub = string(r(lrp_ub), "%9.2f")
        
        * Calculate the area under the curve (AUC)
        quietly roctab `ref_std' `icdcode', nograph
        local auc = string(r(area), "%9.3f")
        local auc_lb = string(r(lb), "%9.3f")
        local auc_ub = string(r(ub), "%9.3f")

        di "----------------"
        di "`icdcode' (95%% CI)"
        di "Sensitivity: `sens' (`sens_lb' - `sens_ub')" 
        di "Specificity: `spec' (`spec_lb' - `spec_ub')" 
        di "Positive LR: `lrp' (`lrp_lb' - `lrp_ub')" 
        di "Positive PV: `ppv' (`ppv_lb' - `ppv_ub')" 
        di "AUC: `auc' (`auc_lb' - `auc_ub')"
        di "----------------"
    }
end




/* ------------------
   Pre-processing 
--------------------*/ 



/* Impute missing abg_ph from serum_hco3 if there is no pH*/ 
//pH = 6.1 + log (HCO3-/ (0.03 x PCO2))
gen predictor_hco3 = abg_hco3
replace predictor_hco3 = serum_hco3 if missing(abg_hco3)
replace abg_ph = 6.1 + log10(predictor_hco3/(0.03 * paco2)) if !missing(predictor_hco3) & !missing(paco2) & missing(abg_ph)
drop predictor_hco3

gen predictor_hco3 = vbg_hco3
replace predictor_hco3 = serum_hco3 if missing(vbg_hco3)
replace vbg_ph = 6.1 + log10(predictor_hco3/(0.03 * vbg_co2)) if !missing(predictor_hco3) & !missing(vbg_co2) & missing(vbg_ph)
drop predictor_hco3

gen abg_vbg_confusion_matrix = 0
replace abg_vbg_confusion_matrix = 1 if has_abg == 1
replace abg_vbg_confusion_matrix = 2 if has_vbg == 1
replace abg_vbg_confusion_matrix = 3 if has_abg == 1 & has_vbg == 1
label variable abg_vbg_confusion_matrix "ABG or VBG Obtained?"
label define abg_vbg_confusion_matrix_label 0 "Neither ABG or VBG obtained on Day 1" 1 "Only ABG obtained on day 1" 2 "Only VBG obtained on day 1" 3 "Both ABG and VBG obtained on day 1"
label values abg_vbg_confusion_matrix abg_vbg_confusion_matrix_label 

tab abg_vbg_confusion_matrix, missing


///-----
//[ ] remove above to preprocessing
////////

replace is_inp = 0 if missing(is_inp) //not sure why is_emer has this missing but is_inp doesn't
replace paco2_flag = 0 if !missing(paco2) & paco2 >= 45 & prim_met_alk == 1
//[ ] TODO: adjust ABG-VBG flag as well

// Restrict to ONLY inpatient and emergency
keep if is_inp == 1 | is_emer == 1



//-----------------
/* Main Analysis */ 
//-----------------

//Summarizing Demographics of the Entire Emulation Cohort

tab first_encounter, missing

table1_mc,  ///
		vars( ///
		age_at_encounter contn %4.1f \ ///
		female bin %4.1f \ ///
		black_race bin %4.1f \ ///
		hisp_eth bin %4.1f \ ///
		bmi contn %4.1f \ ///
		chf bin %4.1f \ ///
		ckd bin %4.1f \ ///
		copd bin %4.1f \ ///
		nmd bin %4.1f \ ///
		osa bin %4.1f \ ///
		has_abg bin %4.1f \ ///
		paco2_flag bin %4.1f \ ///
		has_vbg_and_cat cat %4.1f \ ///
		paco2 contn %4.1f \ ///
		paco2_flag bin %4.1f \ ///
		vbg_or_abg_co2_flag bin %4.1f \ ///
		serum_hco3 contn %4.1f \ ///
		abg_ph contn %4.2f \ ///
		acidemia bin %4.1f \ ///
		is_emer bin %4.1f \ ///
		is_inp bin %4.1f \ ///
		cc_time bin %4.1f \ ///
		vent_proc bin %4.1f \ ///
		died bin %4.1f \ ///
		) ///
		percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") missing onecol saving("Results and Figures/$S_DATE/Overall Cohort chars.xlsx", replace)

//--------------
/* Agreement Analysis */ 
//--------------


// -----------------------------
// SECTION 2: A review of currently available cohort definitions: 
// Systematic Review methods: 
// -----------------------------

/* 
Generate Flags for various study definitions in use 
*/

/* ICU */
//Adler PaCO2 over 47.25 mmHg. NIV or IMV treatment	 Exclude: neuromuscular disease, prognosis < 3 months, iatrogenic hypercapnia, persistent confusion
gen def1 = (paco2 >= 47.25 & vent_proc == 1) if !missing(paco2, vent_proc)
replace def1 = 0 if missing(def1)
label variable def1 "Adler"
label define adler_lab 1 "Adler"
label values def1 adler_lab

//Thille pH<7.35 and PaCO2 45 mmHg treated with NIV or IMV - 0 if any of those are missing.
gen def2 = (paco2 >= 45 & abg_ph < 7.35 & vent_proc == 1) if !missing(paco2, abg_ph, vent_proc)
replace def2 = 0 if missing(def2)
label variable def2 "Thille"
label define thille_lab 1 "Thille"
label values def2 thille_lab

//Ouanes-Besbes pH < 7.35 and PaCO2 > 45 mmHg (no OSA workup, but that's because they're going to be worked up... hm)
gen def3 = (paco2 >= 45 & abg_ph < 7.35) if !missing(paco2, abg_ph)
replace def3 = 0 if missing(def3)
label variable def3 "Ouanes-Besbes, Bülbül"
label define ouanes_lab 1 "Ouanes-Besbes, Bülbül"
label values def3 ouanes_lab

/* Hospitalized / ED  (No ICU requirement) */ 
//Segreelles-Calvo PaCO2 over 45 and pH <7.35  & Received NIV
gen def4 = (def3 == 1 & niv_proc == 1) if !missing(paco2, abg_ph, niv_proc)
replace def4 = 0 if missing(def4)
label variable def4 "Calvo"
label define segreelles_lab 1 "Calvo"
label values def4 segreelles_lab

//Bulbul is the same as Ouanes-Besbes


//Meservey et al 2020:  Admit with code for hypercapnic respiratory failure, exclude Advanced cancer, trauma, stroke, seizure, cardiac arrest, advanced neurologic disease, serious non-pulmonary illness. 
gen def5 = hypercap_resp_failure
replace def5 = 0 if missing(def5)
label variable def5 "Meservey"
label define meservey_lab 1 "Meservey" 
label values def5 meservey_lab

//Vonderbank We preferred capillary blood gas analysis but also included patients with arterial blood gas analysis and some patients with only venous blood gas analysis. (Arterial blood gas analysis is the gold standard in the measurement of blood gases. However, the procedure to obtain arterial blood gas data is painful. Arterialized capillary gases sampled at the ear lobe give similar results for pH and pCO2.2 The interpretation of venous blood gas data is more difficult. The pH is slightly lower (0.02–0.04 pH units) and the pCO2 is slightly higher (3–8 mmHg). However, differences can be greater in patients with hypotension and they depend on local metabolism. We only used venous blood gases if the pCO2 was <45 mmHg and pH was >7.35, which allowed hypercapnia and acidosis to be excluded.3,4 If the pH was also >7.35 and oxygen saturation (measured by pulse oximetry) was normal, additional blood gases were unnecessary but, if not, arterial blood gas data were obtained.)
gen def6 = (paco2 >= 45 & !missing(paco2)) | (vbg_co2 >= 45 & vbg_ph >= 7.35 & !missing(vbg_co2, vbg_ph))
replace def6 = 0 if missing(def6)
label variable def6 "Vonderbank"
label define vonderbank_lab 1 "Vonderbank" 
label values def6 vonderbank_lab

//Wilson et al: PaCO2 over 45 and pH 7.35-45
gen def7 = (paco2 >= 45 & abg_ph >= 7.35 & abg_ph <= 7.45) if !missing(paco2, abg_ph)
replace def7 = 0 if missing(def7)
label variable def7 "Wilson"
label define wilson_lab 1 "Wilson" 
label values def7 wilson_lab

//ED
//Cavalot* pH <7.35 paCO2 >45 or VBG ph <7.34 and CO2 50; excluding cystic fibrosis, neuromuscular disease, ILD, lung cancer, drug overdose, or tracheostomy.
gen def8 = ( (paco2 >= 45 & abg_ph < 7.35) | (vbg_co2 >= 50 & vbg_ph < 7.34) ) if !missing(paco2, abg_ph) | !missing(vbg_co2, vbg_ph)
replace def8 = 0 if missing(def8)
label variable def8 "Cavalot"
label define cavalot_lab 1 "Cavalot"
label values def8 cavalot_lab 

//Chung et al 2021: PaCO2 over 45 - exclude Iatrogenic causes, trauma, post-arrest. 
gen def9 = hypercap_on_abg
replace def9 = 0 if missing(def9)
label variable def9 "Chung"
label define chung_lab 1 "Chung"
label values def9 chung_lab


list hypercap_resp_failure paco2 vbg_co2 vbg_po2 vbg_ph niv_proc imv_proc acidemia def1 def2 def3 def4 def5 def6 def7 def8 def9 in 1/200

/* Ones I can't do */ 
//Domaradzki et al: Diagnostic code for COPD or respiratory failure & VBG (they did not separate into dichotomous groups)
//Vonderbank: Arterialized capillary blood gas CO2 > 45 mmHg
//Bulbul - paco2 at two time points over 45 (can't get this)
// Fox et al - referral on discharge (we can't get this)
// Brandao - NIV outside the ICU (can't get this)
/*
//Gursel Hypercapnia treated with NIV // Remove - as this is only acute outcomes. 
gen def4 = (paco2 >= 45 & vent_proc == 1)
label variable def4 "Gursel"
label define gursel_lab 1 "Gursel"
label values def4 gursel_lab
*/ 
/*
//Marik 2016: BMI over 40, daytime co2 over 45 and hco3 over 28. No ILD, NMD, smoking, COPD .--- to restrictive
gen def10 = (paco2 >= 45 & bmi >= 40 & serum_hco3 > 28 & !(nmd == 1 | nic == 1 | op_nrt == 1 | copd == 1))
label variable def10 "Marik"
label define marik_lab 1 "Marik"
label values def10 marik_lab
*/ 
//Can't do nowbar given prospective assessment 


//10 = the number of definitions
matrix rel_sens = J(9, 9, .)
forval i = 1/9 {
    forval j = 1/9 {
        if `i' == `j' {
            matrix rel_sens[`i', `j'] = .
        }
        else {
            // Generate a temporary variable for combined criteria
            gen temp_both_defi_defj = def`i' == 1 & def`j' == 1 //of course it's 100% 

            // Summarize the temporary variable
            sum temp_both_defi_defj if def`i' == 1, meanonly
            scalar num_both_defi_defj = r(mean)

            // Calculate the percentage and store it in the matrix
            matrix rel_sens[`i', `j'] = num_both_defi_defj * 100

            // Drop the temporary variable
            drop temp_both_defi_defj
        }
    }
}
matrix list rel_sens
matrix rownames rel_sens = "Adler" "Thille" "Ouanes-Besbes, Bülbül" "Calvo" "Meservey*" "Vonderbank" "Wilson" "Cavalot" "Chung**"
matrix colnames rel_sens = "Adler" "Thille" "Ouanes-Besbes, Bülbül" "Calvo" "Meservey*" "Vonderbank" "Wilson" "Cavalot" "Chung**"

heatplot rel_sens, ///
 aspectratio(0.8) ///
 cuts(0(5)100) ///
 xlabel(,angle(45) labsize(medsmall)) ///
 ylabel(,angle(45) labsize(medsmall)) ///
 legend(off) ///
 p(lcolor(black%10) lwidth(*0.15)) ///
 values(format(%4.0f) size(small) color(white)) ///
 ytitle("Reference Standard Definition", size(medsmall)) ///
 xtitle("Comparison Definition", size(medsmall)) ///
 title("Relative Sensitivity of Case Definitions", size(medsmall)) ///
 color(RdYlGn, intensify(1.25)) ///
 xsize(5) ysize(5)
graph export "Results and Figures/$S_DATE/Definition Overlap HeatPlot.png", as(png) name("Graph") replace

//should I flip the side of either the Y or the X-axis? (and then change the angle to -45)

matrix kappa_results = J(9, 9, .)
forval i = 1/9 {
    forval j = 1/9 {
        if `i' <= `j' {
            if `i' == `j' { 
				matrix kappa_results[`i', `j'] = . // Perfect agreement with itself
			} //Otherwise skip, so we're just doing the lower  corners
        }
        else {
            // Calculate kappa for def`i' and def`j'
            kap def`i' def`j'
            
            // Retrieve the kappa statistic and store it in the matrix
            matrix kappa_results[`i', `j'] = r(kappa)
        }
    }
}

// Now, you can display the matrix or use it as needed.
//To fix the hanging labels, might need to drop down to 8 x 8 table ... would take a bit more work

matrix list kappa_results
matrix rownames kappa_results = "Adler" "Thille" "Ouanes-Besbes, Bülbül" "Calvo" "Meservey*" "Vonderbank" "Wilson" "Cavalot" "Chung**"
matrix colnames kappa_results = "Adler" "Thille" "Ouanes-Besbes, Bülbül" "Calvo" "Meservey*" "Vonderbank" "Wilson" "Cavalot" "Chung**"

heatplot kappa_results, ///
 aspectratio(0.65) ///
 cuts(0(0.05)1) ///
 xlabel(,angle(45) labsize(medsmall)) ///
 ylabel(,angle(45) labsize(medsmall)) ///
 legend(off) ///
 p(lcolor(black%10) lwidth(*0.15)) ///
 values(format(%4.2f) size(small) color(white)) ///
 title("Agreement Beyond Chance of Case Definitions", size(medsmall)) ///
 color(RdYlGn, intensify(1.25)) ///
 xsize(5) ysize(5)
graph export "Results and Figures/$S_DATE/Definition Overlap HeatPlot - Kappa.png", as(png) name("Graph") replace

preserve
//Calculate Median and IQR range
svmat kappa_results, name(kappa_value) //makes separate column for each
drop if missing(kappa_value1) //has the most ; this is to make this go fast (or it will try to reshape the whole dataset)
gen id = _n
reshape long kappa_value, i(id) j(column)
summarize kappa_value, detail
sort kappa_value
list kappa_value 
restore

//Generate tables with the characteristics selected in the emulated data-set
*** Must be manually combined? 
//"Adler" "Thille" "Ouanes-Besbes, Bülbül" "Calvo" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"
forval i = 1/9 {
	table1_mc, by(def`i') ///
		vars( ///
		age_at_encounter contn %4.0f \ ///
		female bin %4.0f \ ///
		black_race bin %4.0f \ ///
		hisp_eth bin %4.0f \ ///
		bmi contn %4.0f \ ///
		chf bin %4.0f \ ///
		ckd bin %4.0f \ ///
		copd bin %4.0f \ ///
		nmd bin %4.0f \ ///
		osa bin %4.0f \ ///
		paco2 contn %4.1f \ ///
		serum_hco3 contn %4.1f \ ///
		cc_time bin %4.0f \ ///
		vent_proc bin %4.0f \ ///
		died bin %4.0f \ ///
		) ///
		percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") missing onecol saving("Results and Figures/$S_DATE/Def`i'-Summary.xlsx", replace)
}

bysort died: summ months_death_or_cens, detail

//TODO: figure out how I want to present this info

//Neither, ABG only, VBG only, both -> likelihood of meeting each case definition?


pvenn2 has_abg has_vbg vent_proc
pvenn2 has_abg has_vbg hypercap_resp_failure
pvenn2 hypercap_on_abg hypercap_on_vbg vent_proc
pvenn2 hypercap_on_abg hypercap_on_vbg hypercap_resp_failure
pvenn2 has_abg hypercap_on_abg has_vbg

table1_mc, by(abg_vbg_confusion_matrix) ///
		vars( ///
		is_emer bin %4.0f \ ///
		is_inp bin %4.0f \ ///
		cc_time bin %4.0f \ ///		
		location cat %4.0f \ ///
		def1 bin %4.0f \ ///
		def2 bin %4.0f \ ///
		def3 bin %4.0f \ ///
		def4 bin %4.0f \ ///
		def5 bin %4.0f \ ///
		def6 bin %4.0f \ ///
		def7 bin %4.0f \ ///
		def8 bin %4.0f \ ///
		def9 bin %4.0f \ ///
		) ///
		total(before) percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") missing onecol saving("Results and Figures/$S_DATE/Case Definition by Workup.xlsx", replace)
		
		
/* Regional Variation */ 

table1_mc, by(location) ///
		vars( ///
		def1 bin %4.1f \ ///
		def2 bin %4.1f \ ///
		def3 bin %4.1f \ ///
		def4 bin %4.1f \ ///
		def5 bin %4.1f \ ///
		def6 bin %4.1f \ ///
		def7 bin %4.1f \ ///
		def8 bin %4.1f \ ///
		def9 bin %4.1f \ ///
		) ///
		total(before) percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") missing onecol saving("Results and Figures/$S_DATE/Location by Case Definitions.xlsx", replace)

label list loc_lab

forval z = 0/3 {
	preserve
	keep if location == `z'
	
	matrix loc_kappa_results = J(9, 9, .)
	forval i = 1/9 {
		forval j = 1/9 {
			if `i' <= `j' {
				if `i' == `j' { 
					matrix loc_kappa_results[`i', `j'] = . // Perfect agreement with itself
				} //Otherwise skip, so we're just doing the lower  corners
			}
			else {
				// Calculate kappa for def`i' and def`j'
				quietly kap def`i' def`j'
				
				// Retrieve the kappa statistic and store it in the matrix
				matrix loc_kappa_results[`i', `j'] = r(kappa)
			}
		}
	}

	// Now, you can display the matrix or use it as needed.
	//To fix the hanging labels, might need to drop down to 8 x 8 table ... would take a bit more work
	matrix list loc_kappa_results
	matrix rownames loc_kappa_results = "Adler" "Thille" "Ouanes-Besbes, Bülbül" "Calvo" "Meservey*" "Vonderbank" "Wilson" "Cavalot" "Chung**"
	matrix colnames loc_kappa_results = "Adler" "Thille" "Ouanes-Besbes, Bülbül" "Calvo" "Meservey*" "Vonderbank" "Wilson" "Cavalot" "Chung**"

	heatplot loc_kappa_results, ///
	 aspectratio(0.65) ///
	 cuts(0(0.05)1) ///
	 xlabel(,angle(45) labsize(medsmall)) ///
	 ylabel(,angle(45) labsize(medsmall)) ///
	 legend(off) ///
	 p(lcolor(black%10) lwidth(*0.15)) ///
	 values(format(%4.2f) size(small) color(white)) ///
	 title("Agreement Beyond Chance by Location", size(medsmall)) ///
	 color(RdYlGn, intensify(1.25)) ///
	 xsize(5) ysize(5)
	graph export "Results and Figures/$S_DATE/Loc`z'-Definition Overlap HeatPlot - Kappa.png", as(png) name("Graph") replace


//Calculate Median and IQR range

	svmat loc_kappa_results, name(loc_kappa_value) //makes separate column for each
	drop if missing(loc_kappa_value1) //has the most ; this is to make this go fast (or it will try to reshape the whole dataset)
	gen id = _n
	reshape long loc_kappa_value, i(id) j(column)
	summarize loc_kappa_value, detail
	sort loc_kappa_value
	list loc_kappa_value 
	restore
}

		
		
		
		
		

/* ---------------
PART 2: ICD CODE (Meservey) VS CLINICAL REFERENCE STANDARD (Chung)
-----------------*/

tab has_abg, missing

gen highest_any_flag = 0 
replace highest_any_flag = 1 if (highest_paco2_flag == 1 | highest_vbg_co2_flag == 1)
tab highest_paco2_flag highest_any_flag, missing
tab highest_vbg_co2_flag highest_any_flag, missing

test_char_from_icd_list, ref_std("paco2_flag") icdcodes("ohs_code" "has_j9602" "has_j9612" "has_j9622" "has_j9692" "hypercap_resp_failure")

preserve 
keep if has_abg == 1
tab hypercap_resp_failure if paco2_flag == 1, missing

//How predictive is each hypercapnic respiratory failure diagnostic code? 

test_char_from_icd_list, ref_std("paco2_flag") icdcodes("ohs_code" "has_j9602" "has_j9612" "has_j9622" "has_j9692" "hypercap_resp_failure")
test_char_from_icd_list, ref_std("highest_paco2_flag") icdcodes("ohs_code" "has_j9602" "has_j9612" "has_j9622" "has_j9692" "hypercap_resp_failure")
test_char_from_icd_list, ref_std("highest_any_flag") icdcodes("ohs_code" "has_j9602" "has_j9612" "has_j9622" "has_j9692" "hypercap_resp_failure")
test_char_from_icd_list, ref_std("paco2_flag") icdcodes("has_j9600" "has_j9601" "has_j9610" "has_j9611" "has_j9620" "has_j9621" "has_j9690" "has_j9691") 
test_char_from_icd_list, ref_std("paco2_flag") icdcodes("resp_acid_dx" "sleep_hypovent_dx" "cchs_dx" "other_sleep_hypovent_dx" "other_abn_of_br") 
restore 

/* ------------
//Degree of elevation modeled with a spline. 

//Unadjusted for other characteristics

------------ */ 

//Full Cohort
preserve 
keep if has_abg == 1
gen paco2_rounded = round(paco2, 0.1)
mkspline2 rc = paco2_rounded, cubic nknots(5) displayknots
assert float(paco2_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE

logistic hypercap_resp_failure rc* 
lroc, nograph

levelsof paco2_rounded if paco2_rounded >= 20 & paco2_rounded <= 100, local(levels)
xblc rc*, covname(paco2_rounded) at(`r(levels)') reference(45) eform generate(pa odds lb ub) pr

gen log_odds = ln(odds)
gen log_lb = ln(lb) 
gen log_ub = ln(ub)

gen prob_hypercap = odds / (1+odds)
gen pr_lb = lb / (1+lb)
gen pr_ub = ub / (1+ub)
// note: getting the predicted options to work here (for just predicted odds is hard)

//Title: Probability of a Diagnostic Code for Hypercapnic Respiratory Failure; Unadjusted
twoway (line pr_lb pr_ub pa, sort lc(black black) lp(longdash longdash)) (line prob_hypercap pa, sort lc(black) lp(l)) if inrange(paco2_rounded,20,100), xscale(range(20 100)) yscale(range(0 0.6)) scheme(cleanplots) legend(off) xlabel(20(10)100) xmtick(20(10)100) ylabel(0(0.2)0.6) ytitle(" " " " ) xtitle(" ") title("All Encounters") yline(0, lp("shortdash") lc(gs10)) xline(45, lp("shortdash_dot") lc(gs10)) note(" ") xsize(7) ysize(3)
graph export "Results and Figures/$S_DATE/Unadjusted Prob of Dx Hypercapnia Splines .png", as(png) name("Graph") replace
graph save "All_Encounters_Prob_Dx_spline.gph", replace
restore		


// Emergency Room Encounter
preserve 
keep if has_abg == 1
keep if is_emer == 1
gen paco2_rounded = round(paco2, 0.1)
mkspline2 rc = paco2_rounded, cubic nknots(5) displayknots
assert float(paco2_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE

logistic hypercap_resp_failure rc* 
lroc, nograph

levelsof paco2_rounded if paco2_rounded >= 20 & paco2_rounded <= 100, local(levels)
xblc rc*, covname(paco2_rounded) at(`r(levels)') reference(45) eform generate(pa odds lb ub) pr

gen log_odds = ln(odds)
gen log_lb = ln(lb) 
gen log_ub = ln(ub)

gen prob_hypercap = odds / (1+odds)
gen pr_lb = lb / (1+lb)
gen pr_ub = ub / (1+ub)
// note: getting the predicted options to work here (for just predicted odds is hard)

//Title: Probability of a Diagnostic Code for Hypercapnic Respiratory Failure; Unadjusted
twoway (line pr_lb pr_ub pa, sort lc(black black) lp(longdash longdash)) (line prob_hypercap pa, sort lc(black) lp(l)) if inrange(paco2_rounded,20,100), xscale(range(20 100)) yscale(range(0 0.6)) scheme(cleanplots) legend(off) xlabel(20(10)100) xmtick(20(10)100) ylabel(0(0.2)0.6) ytitle("Probability of a Hypercapnic" "Respiratory Failure Diagnosis Code") xtitle(" ") title("Emergency Encounters") yline(0, lp("shortdash") lc(gs10)) xline(45, lp("shortdash_dot") lc(gs10)) note(" ") xsize(7) ysize(3)
graph export "Results and Figures/$S_DATE/Unadjusted Prob of Dx Hypercapnia Splines .png", as(png) name("Graph") replace
graph save "Emer_Encounters_Prob_Dx_spline.gph", replace
restore		

// Inpatient Encounter
preserve 
keep if has_abg == 1
keep if is_inp == 1
gen paco2_rounded = round(paco2, 0.1)
mkspline2 rc = paco2_rounded, cubic nknots(5) displayknots
assert float(paco2_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE

logistic hypercap_resp_failure rc* 
lroc, nograph

levelsof paco2_rounded if paco2_rounded >= 20 & paco2_rounded <= 100, local(levels)
xblc rc*, covname(paco2_rounded) at(`r(levels)') reference(45) eform generate(pa odds lb ub) pr

gen log_odds = ln(odds)
gen log_lb = ln(lb) 
gen log_ub = ln(ub)

gen prob_hypercap = odds / (1+odds)
gen pr_lb = lb / (1+lb)
gen pr_ub = ub / (1+ub)
// note: getting the predicted options to work here (for just predicted odds is hard)

//Title: Probability of a Diagnostic Code for Hypercapnic Respiratory Failure; Unadjusted
twoway (line pr_lb pr_ub pa, sort lc(black black) lp(longdash longdash)) (line prob_hypercap pa, sort lc(black) lp(l)) if inrange(paco2_rounded,20,100), xscale(range(20 100)) yscale(range(0 0.6)) scheme(cleanplots) legend(off) xlabel(20(10)100) xmtick(20(10)100) ylabel(0(0.2)0.6) ytitle(" " " ") xtitle("Day 1 PaCO{subscript:2}") title("Inpatient Encounters") yline(0, lp("shortdash") lc(gs10)) xline(45, lp("shortdash_dot") lc(gs10)) note(" ") xsize(7) ysize(3)
graph export "Results and Figures/$S_DATE/Unadjusted Prob of Dx Hypercapnia Splines .png", as(png) name("Graph") replace
graph save "Inp_Encounters_Prob_Dx_spline.gph", replace
restore		

graph combine All_Encounters_Prob_Dx_spline.gph Emer_Encounters_Prob_Dx_spline.gph Inp_Encounters_Prob_Dx_spline.gph, ///
	cols(1) /// 
	xcommon ///
	xsize(7) ysize(9)
graph export "Results and Figures/$S_DATE/Figure 2 Prob Hypercap ICD.png", name("Graph") width(3600) replace
//graph export "Results and Figures/$S_DATE/IPW Figure 2.svg", name("Graph") replace //huge




// Location 0 == South
preserve 
keep if has_abg == 1
keep if location == 0
gen paco2_rounded = round(paco2, 0.1)
mkspline2 rc = paco2_rounded, cubic nknots(5) displayknots
assert float(paco2_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE

logistic hypercap_resp_failure rc* 
lroc, nograph

levelsof paco2_rounded if paco2_rounded >= 20 & paco2_rounded <= 80, local(levels)
xblc rc*, covname(paco2_rounded) at(`r(levels)') reference(45) eform generate(pa odds lb ub) pr

gen log_odds = ln(odds)
gen log_lb = ln(lb) 
gen log_ub = ln(ub)

gen prob_hypercap = odds / (1+odds)
gen pr_lb = lb / (1+lb)
gen pr_ub = ub / (1+ub)
// note: getting the predicted options to work here (for just predicted odds is hard)

//Title: Probability of a Diagnostic Code for Hypercapnic Respiratory Failure; Unadjusted
twoway (line pr_lb pr_ub pa, sort lc(black black) lp(longdash longdash)) (line prob_hypercap pa, sort lc(black) lp(l)) if inrange(paco2_rounded,20,80), xscale(range(20 80)) yscale(range(0 0.8)) scheme(cleanplots) legend(off) xlabel(20(10)80) xmtick(20(10)80) ylabel(0(0.2)0.8) ytitle(" " " ") xtitle("Day 1 PaCO{subscript:2}") title("Region: South") yline(0, lp("shortdash") lc(gs10)) xline(45, lp("shortdash_dot") lc(gs10)) note(" ") xsize(4) ysize(4)
graph export "Results and Figures/$S_DATE/South - Unadjusted Prob of Dx Hypercapnia Splines .png", as(png) name("Graph") replace
graph save "Loc0_Encounters_Prob_Dx_spline.gph", replace
restore		


/* 
REGIONAL VARIATION in the application of ICD codes at various levels of PaCO2 elevation
*/ 


// Location 1 = Northeast
preserve 
keep if has_abg == 1
keep if location == 1
gen paco2_rounded = round(paco2, 0.1)
mkspline2 rc = paco2_rounded, cubic nknots(5) displayknots
assert float(paco2_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE

logistic hypercap_resp_failure rc* 
lroc, nograph

levelsof paco2_rounded if paco2_rounded >= 20 & paco2_rounded <= 80, local(levels)
xblc rc*, covname(paco2_rounded) at(`r(levels)') reference(45) eform generate(pa odds lb ub) pr

gen log_odds = ln(odds)
gen log_lb = ln(lb) 
gen log_ub = ln(ub)

gen prob_hypercap = odds / (1+odds)
gen pr_lb = lb / (1+lb)
gen pr_ub = ub / (1+ub)
// note: getting the predicted options to work here (for just predicted odds is hard)

//Title: Probability of a Diagnostic Code for Hypercapnic Respiratory Failure; Unadjusted
twoway (line pr_lb pr_ub pa, sort lc(black black) lp(longdash longdash)) (line prob_hypercap pa, sort lc(black) lp(l)) if inrange(paco2_rounded,20,80), xscale(range(20 80)) yscale(range(0 0.8)) scheme(cleanplots) legend(off) xlabel(20(10)80) xmtick(20(10)80) ylabel(0(0.2)0.8) ytitle(" " " ") xtitle(" ") title("Region: Northeast") yline(0, lp("shortdash") lc(gs10)) xline(45, lp("shortdash_dot") lc(gs10)) note(" ") xsize(4) ysize(4)
graph export "Results and Figures/$S_DATE/Northeast - Unadjusted Prob of Dx Hypercapnia Splines .png", as(png) name("Graph") replace
graph save "Loc1_Encounters_Prob_Dx_spline.gph", replace
restore		

// Location 2 = Midwest
preserve 
keep if has_abg == 1
keep if location == 3
gen paco2_rounded = round(paco2, 0.1)
mkspline2 rc = paco2_rounded, cubic nknots(5) displayknots
assert float(paco2_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE

logistic hypercap_resp_failure rc* 
lroc, nograph

levelsof paco2_rounded if paco2_rounded >= 20 & paco2_rounded <= 80, local(levels)
xblc rc*, covname(paco2_rounded) at(`r(levels)') reference(45) eform generate(pa odds lb ub) pr

gen log_odds = ln(odds)
gen log_lb = ln(lb) 
gen log_ub = ln(ub)

gen prob_hypercap = odds / (1+odds)
gen pr_lb = lb / (1+lb)
gen pr_ub = ub / (1+ub)
// note: getting the predicted options to work here (for just predicted odds is hard)

//Title: Probability of a Diagnostic Code for Hypercapnic Respiratory Failure; Unadjusted
twoway (line pr_lb pr_ub pa, sort lc(black black) lp(longdash longdash)) (line prob_hypercap pa, sort lc(black) lp(l)) if inrange(paco2_rounded,20,80), xscale(range(20 80)) yscale(range(0 0.8)) scheme(cleanplots) legend(off) xlabel(20(10)80) xmtick(20(10)80) ylabel(0(0.2)0.8) ytitle("Probability of" "Diagnosis Code") xtitle(" ") title("Region: Midwest") yline(0, lp("shortdash") lc(gs10)) xline(45, lp("shortdash_dot") lc(gs10)) note(" ") xsize(4) ysize(4)
graph export "Results and Figures/$S_DATE/Midwest - Unadjusted Prob of Dx Hypercapnia Splines .png", as(png) name("Graph") replace
graph save "Loc2_Encounters_Prob_Dx_spline.gph", replace
restore		


// Location 3 = West
preserve 
keep if has_abg == 1
keep if location == 3
gen paco2_rounded = round(paco2, 0.1)
mkspline2 rc = paco2_rounded, cubic nknots(5) displayknots
assert float(paco2_rounded) == float(rc1) // VERIFY FIRST SPLINE VARIABLE IS THE ORIGINAL VARIABLE

logistic hypercap_resp_failure rc* 
lroc, nograph

levelsof paco2_rounded if paco2_rounded >= 20 & paco2_rounded <= 80, local(levels)
xblc rc*, covname(paco2_rounded) at(`r(levels)') reference(45) eform generate(pa odds lb ub) pr

gen log_odds = ln(odds)
gen log_lb = ln(lb) 
gen log_ub = ln(ub)

gen prob_hypercap = odds / (1+odds)
gen pr_lb = lb / (1+lb)
gen pr_ub = ub / (1+ub)
// note: getting the predicted options to work here (for just predicted odds is hard)

//Title: Probability of a Diagnostic Code for Hypercapnic Respiratory Failure; Unadjusted
twoway (line pr_lb pr_ub pa, sort lc(black black) lp(longdash longdash)) (line prob_hypercap pa, sort lc(black) lp(l)) if inrange(paco2_rounded,20,80), xscale(range(20 80)) yscale(range(0 0.8)) scheme(cleanplots) legend(off) xlabel(20(10)80) xmtick(20(10)80) ylabel(0(0.2)0.8) ytitle("Probability of" "Diagnosis Code") xtitle("Day 1 PaCO{subscript:2}") title("Region: West") yline(0, lp("shortdash") lc(gs10)) xline(45, lp("shortdash_dot") lc(gs10)) note(" ") xsize(4) ysize(4)
graph export "Results and Figures/$S_DATE/West - Unadjusted Prob of Dx Hypercapnia Splines.png", as(png) name("Graph") replace
graph save "Loc3_Encounters_Prob_Dx_spline.gph", replace
restore		



graph combine Loc2_Encounters_Prob_Dx_spline.gph Loc1_Encounters_Prob_Dx_spline.gph Loc3_Encounters_Prob_Dx_spline.gph Loc0_Encounters_Prob_Dx_spline.gph, ///
	cols(2) /// 
	xcommon ///
	ycommon ///
	xsize(8) ysize(8)
graph export "Results and Figures/$S_DATE/Location - Figure S3 Prob Hypercap ICD.png", name("Graph") width(3200) replace

