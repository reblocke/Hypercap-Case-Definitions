version 17.0
args harness_root analysis_root input_root output_root run_id status_file ///
    dependency_report run_mode

capture log close _all
clear all

local dependency_rc = 0
local analysis_rc = .
local marker_found = 0
local status "failed_driver"

capture noisily do "`harness_root'/stata/preflight_dependencies.do" ///
    "`dependency_report'"
local dependency_rc = _rc

if `dependency_rc' == 0 {
    capture cd "`analysis_root'"
    local cd_rc = _rc
    if `cd_rc' == 0 {
        if "`run_mode'" == "legacy" {
            capture noisily do "Hypercapnia Case Definitions.do" ///
                "`input_root'" "`output_root'"
        }
        else {
            capture noisily do "Hypercapnia Case Definitions.do" ///
                "`input_root'" "`output_root'" "`run_id'"
        }
        local analysis_rc = _rc
    }
    else {
        local analysis_rc = `cd_rc'
    }
}

capture log close _all

if `dependency_rc' != 0 {
    local status "failed_dependency"
}
else if `analysis_rc' != 0 {
    local status "failed_analysis"
}
else if "`run_mode'" == "candidate" {
    capture confirm file "`output_root'/`run_id'/ANALYSIS_COMPLETE"
    local marker_found = (_rc == 0)
    if `marker_found' local status "success"
    else local status "failed_incomplete"
}
else {
    local status "success"
    local marker_found = 1
}

local status_tmp "`status_file'.tmp"
tempname status_handle
file open `status_handle' using "`status_tmp'", write text replace
file write `status_handle' "status" _tab "`status'" _n
file write `status_handle' "dependency_rc" _tab "`dependency_rc'" _n
file write `status_handle' "analysis_rc" _tab "`analysis_rc'" _n
file write `status_handle' "completion_marker_found" _tab "`marker_found'" _n
file write `status_handle' "stata_version" _tab "`c(stata_version)'" _n
file write `status_handle' "stata_flavor" _tab "`c(flavor)'" _n
file write `status_handle' "operating_system" _tab "`c(os)'" _n
file write `status_handle' "machine_type" _tab "`c(machine_type)'" _n
file close `status_handle'
copy "`status_tmp'" "`status_file'", replace
erase "`status_tmp'"

exit 0
