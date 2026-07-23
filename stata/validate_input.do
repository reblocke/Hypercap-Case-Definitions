version 17.0
args contract_file report_file

capture confirm file "`contract_file'"
if _rc {
    di as error "Input contract not found."
    exit 601
}

capture frame drop hcd_contract
frame create hcd_contract
capture frame hcd_contract: import delimited using "`contract_file'", ///
    varnames(1) stringcols(_all) clear
local import_rc = _rc
if `import_rc' {
    capture frame drop hcd_contract
    di as error "Could not read the input contract."
    exit `import_rc'
}

frame hcd_contract: levelsof variable_name if workflow_role == "runtime_input", ///
    local(required_vars) clean
frame hcd_contract: levelsof variable_name if workflow_role == "runtime_input" & ///
    type == "binary", local(binary_vars) clean
frame hcd_contract: levelsof variable_name if workflow_role == "runtime_input" & ///
    type == "continuous", local(continuous_vars) clean

tempname report
file open `report' using "`report_file'", write text replace
file write `report' "check" _tab "variable" _tab "severity" _tab ///
    "count" _tab "detail" _n

local hard_failures = 0
foreach variable of local required_vars {
    capture confirm variable `variable'
    if _rc {
        file write `report' "required_variable" _tab "`variable'" _tab ///
            "error" _tab "1" _tab "missing" _n
        local ++hard_failures
    }
    else {
        capture confirm numeric variable `variable'
        if _rc {
            file write `report' "numeric_type" _tab "`variable'" _tab ///
                "error" _tab "1" _tab "nonnumeric" _n
            local ++hard_failures
        }
        else {
            file write `report' "required_variable" _tab "`variable'" _tab ///
                "pass" _tab "0" _tab "present_numeric" _n
        }
    }
}

if `hard_failures' == 0 {
    foreach variable of local binary_vars {
        quietly count if !missing(`variable') & !inlist(`variable', 0, 1)
        local violations = r(N)
        local severity = cond(`violations' == 0, "pass", "error")
        file write `report' "binary_domain" _tab "`variable'" _tab ///
            "`severity'" _tab "`violations'" _tab "allowed_0_1_or_missing" _n
        if `violations' > 0 local ++hard_failures
    }

    quietly count if !missing(location) & ///
        (location < 0 | location > 3 | location != floor(location))
    local violations = r(N)
    local severity = cond(`violations' == 0, "pass", "error")
    file write `report' "location_domain" _tab "location" _tab ///
        "`severity'" _tab "`violations'" _tab "integer_0_to_3_or_missing" _n
    if `violations' > 0 local ++hard_failures

    quietly count if is_inp == 1 | is_emer == 1
    local retained = r(N)
    local severity = cond(`retained' > 0, "pass", "error")
    file write `report' "analytic_cohort" _tab "_rows" _tab ///
        "`severity'" _tab "`retained'" _tab "ed_or_inpatient_before_mutation" _n
    if `retained' == 0 local ++hard_failures

    foreach variable of local required_vars {
        quietly count if missing(`variable')
        if r(N) > 0 {
            file write `report' "missingness" _tab "`variable'" _tab ///
                "warning" _tab "`r(N)'" _tab "current_code_may_treat_as_negative_or_exclude" _n
        }
    }

    quietly count if has_abg == 0 & !missing(paco2)
    file write `report' "availability_discordance" _tab "has_abg_paco2" _tab ///
        "warning" _tab "`r(N)'" _tab "has_abg_0_with_paco2_present" _n
    quietly count if has_abg == 1 & missing(paco2)
    file write `report' "availability_discordance" _tab "has_abg_paco2" _tab ///
        "warning" _tab "`r(N)'" _tab "has_abg_1_with_paco2_missing" _n
    quietly count if has_vbg == 0 & (!missing(vbg_co2) | !missing(vbg_ph))
    file write `report' "availability_discordance" _tab "has_vbg_values" _tab ///
        "warning" _tab "`r(N)'" _tab "has_vbg_0_with_vbg_value_present" _n
    quietly count if has_vbg == 1 & missing(vbg_co2) & missing(vbg_ph)
    file write `report' "availability_discordance" _tab "has_vbg_values" _tab ///
        "warning" _tab "`r(N)'" _tab "has_vbg_1_with_both_values_missing" _n

    foreach variable of local continuous_vars {
        quietly summarize `variable', meanonly
        local observed_min = cond(r(N) == 0, "missing", string(r(min), "%21.0g"))
        local observed_max = cond(r(N) == 0, "missing", string(r(max), "%21.0g"))
        file write `report' "continuous_range" _tab "`variable'" _tab ///
            "warning" _tab "0" _tab ///
            `"range_unresolved_min_`observed_min'_max_`observed_max'"' _n
    }
}

file close `report'
capture frame drop hcd_contract

if `hard_failures' > 0 {
    di as error "Input validation failed with `hard_failures' hard error(s)."
    exit 459
}

di as result "Input contract validation passed."
