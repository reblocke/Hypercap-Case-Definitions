* Hypercapnia TriNetX Computable Phenotype Analysis

/* ---------
OVERALL TODO LIST
---
-----------*/

capture log close

* Load data
clear

/* Specify working direction */ 
//cd "C:\Users\reblo\Box\Residency Personal Files\Scholarly Work\Locke Research Projects\TriNetX Code" 
cd "/Users/blocke/Box Sync/Residency Personal Files/Scholarly Work/Locke Research Projects/TriNetX Code"
//cd "/Users/reblocke/Research/trinetx-hypercapnia-code"

/* Create logging / output directories */ 
capture mkdir "Results and Figures"
capture mkdir "Results and Figures/$S_DATE/" //make new folder for figure output if needed
capture mkdir "Results and Figures/$S_DATE/Logs/" //new folder for stata logs
local a1=substr(c(current_time),1,2)
local a2=substr(c(current_time),4,2)
local a3=substr(c(current_time),7,2)
local b = "Hypercapnia Computable Phenotype.do" // do file name
copy "`b'" "Results and Figures/$S_DATE/Logs/(`a1'_`a2'_`a3')`b'"

set scheme cleanplots
graph set window fontface "Helvetica"
log using temp.log, replace

clear
cd "Data"
use full_db
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

replace is_inp = 0 if missing(is_inp) 
// Restrict to ONLY inpatient and emergency
keep if is_inp == 1 | is_emer == 1

// Responses to Reviewers: 
/*
1.	How often primary metabolic alkalosis with respiratory compensation exclusively explains the elevated PaCO2
2.	How often alkalemia is present among patients with PaCO2 ≥ 45 mmHg, suggesting a metabolic alkalosis is present and may at least partially contribute to the elevated PaCO2. 
*/ 

//Number of primary metabolic conditions. 
tab prim_met_alk
tab combo_met_alk

//Testing Strategy label
gen abg_vbg_confusion_matrix = 0  // Default to "Has neither"
replace abg_vbg_confusion_matrix = 1 if has_abg == 1 & has_vbg == 0  // "Has ABG"
replace abg_vbg_confusion_matrix = 2 if has_abg == 0 & has_vbg == 1  // "Has VBG"
replace abg_vbg_confusion_matrix = 3 if has_abg == 1 & has_vbg == 1  // "Has Both"
label define abg_vbg_labels 0 "Has neither" 1 "Only ABG obtained (first day)" 2 "Only VBG obtained (first day)" 3 "Both ABG and VBG obtained (first day)"
label values abg_vbg_confusion_matrix abg_vbg_labels

missings table paco2
missings table vbg_co2

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

/* ICU (not simulated)*/
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

//Ouanes-Besbes pH < 7.35 and PaCO2 > 45 mmHg (& no OSA workup - non-simulated)
gen def3 = (paco2 >= 45 & abg_ph < 7.35) if !missing(paco2, abg_ph)
replace def3 = 0 if missing(def3)
label variable def3 "Ouanes-Besbes"
label define ouanes_lab 1 "Ouanes-Besbes"
label values def3 ouanes_lab

/* Hospitalized / ED  (No ICU requirement) */ 
//Segreelles-Calvo PaCO2 over 45 and pH <7.35  & Received NIV
gen def4 = (def3 == 1 & niv_proc == 1) if !missing(paco2, abg_ph, niv_proc)
replace def4 = 0 if missing(def4)
label variable def4 "Calvo"
label define segreelles_lab 1 "Calvo"
label values def4 segreelles_lab

//Reviewers Bulbul is NOT the same as Ouanes-Besbes **
//Bulbul - just PaCO2 >= 45 mmHg 
gen def5 = hypercap_on_abg
replace def5 = 0 if missing(def5)
label variable def5 "Bülbül"
label define bulbul_lab 1 "Bülbül"
label values def5 bulbul_lab


//Meservey et al 2020:  Admit with code for hypercapnic respiratory failure, exclude Advanced cancer, trauma, stroke, seizure, cardiac arrest, advanced neurologic disease, serious non-pulmonary illness. 
gen def6 = hypercap_resp_failure
replace def6 = 0 if missing(def6)
label variable def6 "Meservey"
label define meservey_lab 1 "Meservey" 
label values def6 meservey_lab

//Vonderbank We preferred capillary blood gas analysis but also included patients with arterial blood gas analysis and some patients with only venous blood gas analysis. (Arterial blood gas analysis is the gold standard in the measurement of blood gases. However, the procedure to obtain arterial blood gas data is painful. Arterialized capillary gases sampled at the ear lobe give similar results for pH and pCO2.2 The interpretation of venous blood gas data is more difficult. The pH is slightly lower (0.02–0.04 pH units) and the pCO2 is slightly higher (3–8 mmHg). However, differences can be greater in patients with hypotension and they depend on local metabolism. We only used venous blood gases if the pCO2 was <45 mmHg and pH was >7.35, which allowed hypercapnia and acidosis to be excluded.3,4 If the pH was also >7.35 and oxygen saturation (measured by pulse oximetry) was normal, additional blood gases were unnecessary but, if not, arterial blood gas data were obtained.)
gen def7 = (paco2 >= 45 & !missing(paco2)) | (vbg_co2 >= 45 & vbg_ph >= 7.35 & !missing(vbg_co2, vbg_ph))
replace def7 = 0 if missing(def7)
label variable def7 "Vonderbank"
label define vonderbank_lab 1 "Vonderbank" 
label values def7 vonderbank_lab

//Wilson et al: PaCO2 over 45 and pH 7.35-45
gen def8 = (paco2 >= 45 & abg_ph >= 7.35 & abg_ph <= 7.45) if !missing(paco2, abg_ph)
replace def8 = 0 if missing(def8)
label variable def8 "Wilson"
label define wilson_lab 1 "Wilson" 
label values def8 wilson_lab

//ED
//Cavalot* pH <7.35 paCO2 >45 or VBG ph <7.34 and CO2 50; excluding cystic fibrosis, neuromuscular disease, ILD, lung cancer, drug overdose, or tracheostomy.
gen def9 = ((paco2 >= 45 & abg_ph <= 7.35) | (vbg_co2 >= 50 & vbg_ph <= 7.34) ) if !missing(paco2, abg_ph) | !missing(vbg_co2, vbg_ph)
replace def9 = 0 if missing(def9)
label variable def9 "Cavalot"
label define cavalot_lab 1 "Cavalot"
label values def8 cavalot_lab 

//Chung et al 2021: PaCO2 over 45 - exclude Iatrogenic causes, trauma, post-arrest. 
gen def10 = (paco2 >= 45 & abg_ph < 7.45) if !missing(paco2, abg_ph)
replace def10 = 0 if missing(def10)
label variable def10 "Chung"
label define chung_lab 1 "Chung"
label values def10 chung_lab


list hypercap_resp_failure paco2 vbg_co2 vbg_po2 vbg_ph niv_proc imv_proc acidemia def1 def2 def3 def4 def5 def6 def7 def8 def9 def10 in 1/200

/* Ones I can't do */ 
//Domaradzki et al: Diagnostic code for COPD or respiratory failure & VBG (they did not separate into dichotomous groups) - not actually include as not general. 
//Vonderbank: Arterialized capillary blood gas CO2 > 45 mmHg
//Bulbul - paco2 at two time points over 45 (can't get this)
// Fox et al - referral on discharge (we can't get this)
// Brandao - NIV outside the ICU (can't get this)
//Can't do nowbar given prospective assessment 
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

//10 = the number of definitions
matrix rel_sens = J(10, 10, .)
forval i = 1/10 {
    forval j = 1/10 {
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
matrix rownames rel_sens = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"
matrix colnames rel_sens = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"

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

matrix agreement_results = J(10, 10, .)
matrix kappa_results = J(10, 10, .)
matrix pabak_results = J(10, 10, .)

forval i = 1/10 {
    forval j = 1/10 {
        if `i' <= `j' {
            if `i' == `j' { 
				matrix agreement_results[`i', `j'] = . // Perfect agreement with itself
				matrix kappa_results[`i', `j'] = . // Perfect agreement with itself
				matrix pabak_results[`i', `j'] = . // Perfect agreement with itself
			} //Otherwise skip, so we're just doing the lower  corners
        }
        else {
            // Calculate kappa for def`i' and def`j'
			kap def`i' def`j'
            // Retrieve the kappa statistic and store it in the matrix
            matrix kappa_results[`i', `j'] = r(kappa)
			matrix agreement_results[`i', `j'] = r(prop_o)
			
			//PABAK
			kappaetc def`i' def`j'
			matrix pabak_results[`i', `j'] = r(b)[1,2]   // Brennan-Prediger aka PABAK is the 2nd coefficient
        }
    }
}

/* Kappa */
matrix list kappa_results
matrix rownames kappa_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"
matrix colnames kappa_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"

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

/* Raw Agreement */ 
matrix list agreement_results
matrix rownames agreement_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"
matrix colnames agreement_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"

heatplot agreement_results, ///
 aspectratio(0.65) ///
 cuts(0(0.05)1) ///
 xlabel(,angle(45) labsize(medsmall)) ///
 ylabel(,angle(45) labsize(medsmall)) ///
 legend(off) ///
 p(lcolor(black%10) lwidth(*0.15)) ///
 values(format(%4.2f) size(small) color(white)) ///
 title("Raw Agreement Between of Case Definitions", size(medsmall)) ///
 color(RdYlGn, intensify(1.25)) ///
 xsize(5) ysize(5)
graph export "Results and Figures/$S_DATE/Definition Overlap HeatPlot - Agreement.png", as(png) name("Graph") replace

/* PABAK */ 
matrix list pabak_results
matrix rownames pabak_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"
matrix colnames pabak_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"

heatplot pabak_results, ///
 aspectratio(0.65) ///
 cuts(0(0.05)1) ///
 xlabel(,angle(45) labsize(medsmall)) ///
 ylabel(,angle(45) labsize(medsmall)) ///
 legend(off) ///
 p(lcolor(black%10) lwidth(*0.15)) ///
 values(format(%4.2f) size(small) color(white)) ///
 title("Prevalence and Bias-Adjusted Kappa between Case Definitions", size(small)) ///
 color(RdYlGn, intensify(1.25)) ///
 xsize(5) ysize(5)
graph export "Results and Figures/$S_DATE/Definition Overlap HeatPlot - PABAK.png", as(png) name("Graph") replace

/* 
Alternative method of multi-rater agreement evaluation 
*/ 
kappaetc def1 def2 def3 def4 def5 def6 def7 def8 def9 def10 // Cohen would have been an alterantive way to Median Kappa

preserve
//Calculate Median and IQR range
svmat kappa_results, name(kappa_value) //makes separate column for each
drop if missing(kappa_value1) //has the most ; this is to make this go fast (or it will try to reshape the whole dataset)
gen id = _n
reshape long kappa_value, i(id) j(column)
summarize kappa_value, detail // main result
sort kappa_value
list kappa_value 
summarize kappa_value if kappa_value > 0, detail // sensitivity excluding 4 structurally conflicting case definitions.
restore

preserve
//Calculate Median and IQR range
svmat agreement_results, name(agreement_value) //makes separate column for each
drop if missing(agreement_value1) //has the most ; this is to make this go fast (or it will try to reshape the whole dataset)
gen id = _n
reshape long agreement_value, i(id) j(column)
summarize agreement_value, detail // raw result
sort agreement_value
list agreement_value 
summarize agreement_value if agreement_value > 0, detail // sensitivity excluding 4 structurally conflicting case definitions.
restore

//PABAK sensitivity analyses: 
preserve
svmat pabak_results, name(pabak_value) //makes separate column for each
drop if missing(pabak_value1) //has the most ; this is to make this go fast (or it will try to reshape the whole dataset)
gen id = _n
reshape long pabak_value, i(id) j(column)
summarize pabak_value, detail // main sensitivity result
restore

/* Sensitivity Analyses */ 

/* --------------------------
Analysis By testing Strategy: 
ABG only vs VBG only vs ABG and VBG 
abg_vbg_labels 0 "Has neither" 1 "Only ABG obtained (first day)" 2 "Only VBG obtained (first day)" 3 "Both ABG and VBG obtained (first day)"
-------------------------- */ 
//ABG only

/* 
Has ABG (whether or not also had VBG)
Had VBG (whether or not also had ABG)
Had either ABG or VBG
*/ 

//Had an ABG
preserve
keep if has_abg == 1
count 
matrix abg_kappa_results = J(10, 10, .)
forval i = 1/10 {
    forval j = 1/10 {
        if `i' <= `j' {
            if `i' == `j' { 
				matrix abg_kappa_results[`i', `j'] = . // Perfect agreement with itself
			} //Otherwise skip, so we're just doing the lower  corners
        }
        else {
            // Calculate kappa for def`i' and def`j'
			kap def`i' def`j'
            // Retrieve the kappa statistic and store it in the matrix
            matrix abg_kappa_results[`i', `j'] = r(kappa)
        }
    }
}

// Now, you can display the matrix or use it as needed.
//To fix the hanging labels, might need to drop down to 8 x 8 table ... would take a bit more work

matrix list abg_kappa_results
matrix rownames abg_kappa_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"
matrix colnames abg_kappa_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"

heatplot abg_kappa_results, ///
 aspectratio(0.65) ///
 cuts(0(0.05)1) ///
 xlabel(,angle(45) labsize(medsmall)) ///
 ylabel(,angle(45) labsize(medsmall)) ///
 legend(off) ///
 p(lcolor(black%10) lwidth(*0.15)) ///
 values(format(%4.2f) size(small) color(white)) ///
 title("Agreement Beyond Chance of Case Definitions, ABG-only", size(medsmall)) ///
 color(RdYlGn, intensify(1.25)) ///
 xsize(5) ysize(5)
graph export "Results and Figures/$S_DATE/ABG-only Definition Overlap HeatPlot - Kappa.png", as(png) name("Graph") replace

kappaetc def1 def2 def3 def4 def5 def6 def7 def8 def9 def10 // Cohen would have been an alterantive way to Median Kappa

//Calculate Median and IQR range
svmat abg_kappa_results, name(abg_kappa_value) //makes separate column for each
drop if missing(abg_kappa_value1) //has the most ; this is to make this go fast (or it will try to reshape the whole dataset)
gen id = _n
reshape long abg_kappa_value, i(id) j(column)
summarize abg_kappa_value, detail // main result
sort abg_kappa_value
list abg_kappa_value 
restore

//had a vbg
preserve
keep if has_vbg == 1
count 
matrix vbg_kappa_results = J(10, 10, .)
forval i = 1/10 {
    forval j = 1/10 {
        if `i' <= `j' {
            if `i' == `j' { 
				matrix vbg_kappa_results[`i', `j'] = . // Perfect agreement with itself
			} //Otherwise skip, so we're just doing the lower  corners
        }
        else {
            // Calculate kappa for def`i' and def`j'
			kap def`i' def`j'
            // Retrieve the kappa statistic and store it in the matrix
            matrix vbg_kappa_results[`i', `j'] = r(kappa)
        }
    }
}

matrix list vbg_kappa_results
matrix rownames vbg_kappa_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"
matrix colnames vbg_kappa_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"

heatplot vbg_kappa_results, ///
 aspectratio(0.65) ///
 cuts(0(0.05)1) ///
 xlabel(,angle(45) labsize(medsmall)) ///
 ylabel(,angle(45) labsize(medsmall)) ///
 legend(off) ///
 p(lcolor(black%10) lwidth(*0.15)) ///
 values(format(%4.2f) size(small) color(white)) ///
 title("Agreement Beyond Chance of Case Definitions, VBG-only", size(medsmall)) ///
 color(RdYlGn, intensify(1.25)) ///
 xsize(5) ysize(5)
graph export "Results and Figures/$S_DATE/VBG-only Definition Overlap HeatPlot - Kappa.png", as(png) name("Graph") replace

//Calculate Median and IQR range
svmat vbg_kappa_results, name(vbg_kappa_value) //makes separate column for each
drop if missing(vbg_kappa_value1) //has the most ; this is to make this go fast (or it will try to reshape the whole dataset)
gen id = _n
reshape long vbg_kappa_value, i(id) j(column)
summarize vbg_kappa_value, detail // main result
sort vbg_kappa_value
list vbg_kappa_value 
restore

//either ABG or VBG
preserve
keep if abg_vbg_confusion_matrix !=0
count
matrix abg_vbg_kappa_results = J(10, 10, .)
forval i = 1/10 {
    forval j = 1/10 {
        if `i' <= `j' {
            if `i' == `j' { 
				matrix abg_vbg_kappa_results[`i', `j'] = . // Perfect agreement with itself
			} //Otherwise skip, so we're just doing the lower  corners
        }
        else {
            // Calculate kappa for def`i' and def`j'
			kap def`i' def`j'
            // Retrieve the kappa statistic and store it in the matrix
            matrix abg_vbg_kappa_results[`i', `j'] = r(kappa)
        }
    }
}

matrix list abg_vbg_kappa_results
matrix rownames abg_vbg_kappa_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"
matrix colnames abg_vbg_kappa_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"

heatplot abg_vbg_kappa_results, ///
 aspectratio(0.65) ///
 cuts(0(0.05)1) ///
 xlabel(,angle(45) labsize(medsmall)) ///
 ylabel(,angle(45) labsize(medsmall)) ///
 legend(off) ///
 p(lcolor(black%10) lwidth(*0.15)) ///
 values(format(%4.2f) size(small) color(white)) ///
 title("Agreement Beyond Chance of Case Definitions, ABG-VBG", size(medsmall)) ///
 color(RdYlGn, intensify(1.25)) ///
 xsize(5) ysize(5)
graph export "Results and Figures/$S_DATE/ABG-VBG Definition Overlap HeatPlot - Kappa.png", as(png) name("Graph") replace

kappaetc def1 def2 def3 def4 def5 def6 def7 def8 def9 def10 // Cohen would have been an alterantive way to Median Kappa

//Calculate Median and IQR range
svmat abg_vbg_kappa_results, name(abg_vbg_kappa_value) //makes separate column for each
drop if missing(abg_vbg_kappa_value1) //has the most ; this is to make this go fast (or it will try to reshape the whole dataset)
gen id = _n
reshape long abg_vbg_kappa_value, i(id) j(column)
summarize abg_vbg_kappa_value, detail // main result
sort abg_vbg_kappa_value
list abg_vbg_kappa_value 
restore


/* Evaluation of the composition of each group */ 

recode months_death_or_cens (min/-1=.) (17/max=.)  //remove impossible values
summ months_death_or_cens, detail

gen died_1mo = 0
replace died_1mo = 1 if died == 1 & months_death_or_cens < 1 //died in month 0 
label variable died_1mo "Died (1 month)"

gen died_2mo = 0
replace died_2mo = 1 if died == 1 & months_death_or_cens <= 1 //died in month 0 
label variable died_2mo "Died (2 months)"

gen death_time = . 
replace death_time = months_death_or_cens if died == 1
label variable death_time "Month of Death"

stset months_death_or_cens, origin(min) failure(died == 1)

//Generate tables with the characteristics selected in the emulated data-set
*** Must be manually combined? 
//"Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey" "Vonderbank" "Wilson" "Cavalot" "Chung"
forval i = 1/10 {
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
		died_1mo bin %4.0f \ ///
		died_2mo bin %4.0f \ ///
		death_time conts %4.0f \ ///
		months_death_or_cens conts %4.0f \ ///
		died bin %4.0f \ ///
		) ///
		percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") missing onecol saving("Results and Figures/$S_DATE/Def`i'-Summary.xlsx", replace)
	
	stcox i.def`i'
	sts test def`i', logrank 
	sts graph, by(def`i') ///
		risktable(0(2.5)15) ///
		xlabel(0(2.5)15) ///
		xtitle("Follow-up (months)") ///
		ytitle("Probability of Survival") ///
		legend(ring(0) position(6) rows(2)) ///
		xsize(10) ysize(5)
}

bysort died: summ months_death_or_cens, detail // median follow-up is potentially misleadsing 

/* Diagnostic Strategy Influence on Inclusion */ 

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
		def10 bin %4.0f \ ///
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
		def10 bin %4.1f \ ///
		) ///
		total(before) percent_n percsign("%") iqrmiddle(",") sdleft(" (±") sdright(")") missing onecol saving("Results and Figures/$S_DATE/Location by Case Definitions.xlsx", replace)

label list loc_lab

forval z = 0/3 {
	preserve
	keep if location == `z'
	
	matrix loc_kappa_results = J(10, 10, .)
	forval i = 1/10 {
		forval j = 1/10 {
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
	matrix rownames loc_kappa_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey*" "Vonderbank" "Wilson" "Cavalot" "Chung**"
	matrix colnames loc_kappa_results = "Adler" "Thille" "Ouanes-Besbes" "Calvo" "Bülbül" "Meservey*" "Vonderbank" "Wilson" "Cavalot" "Chung**"

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

/* Figure */ 
graph combine All_Encounters_Prob_Dx_spline.gph Emer_Encounters_Prob_Dx_spline.gph Inp_Encounters_Prob_Dx_spline.gph, ///
	cols(1) /// 
	xcommon ///
	xsize(7) ysize(9)
graph export "Results and Figures/$S_DATE/Figure 2 Prob Hypercap ICD.png", name("Graph") width(3600) replace
//graph export "Results and Figures/$S_DATE/IPW Figure 2.svg", name("Graph") replace //huge

/* 
REGIONAL VARIATION in the application of ICD codes at various levels of PaCO2 elevation
*/ 


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



