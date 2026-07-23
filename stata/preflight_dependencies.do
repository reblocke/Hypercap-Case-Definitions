version 17.0
args report_file

tempname report
file open `report' using "`report_file'", write text replace
file write `report' "name" _tab "kind" _tab "required" _tab "found" _tab "resolved_path" _n

local missing_required

foreach command in missings table1_mc heatplot kappaetc diagt mkspline2 xblc colorpalette {
    capture which `command'
    local found = (_rc == 0)
    local resolved
    if `found' {
        capture findfile `command'.ado
        if !_rc local resolved `"`r(fn)'"'
    }
    file write `report' "`command'" _tab "command" _tab "true" _tab ///
        "`=cond(`found', "true", "false")'" _tab `"`resolved'"' _n
    if !`found' local missing_required "`missing_required' `command'"
}

capture findfile scheme-cleanplots.scheme
local found = (_rc == 0)
local resolved
if `found' local resolved `"`r(fn)'"'
file write `report' "cleanplots" _tab "scheme" _tab "true" _tab ///
    "`=cond(`found', "true", "false")'" _tab `"`resolved'"' _n
if !`found' local missing_required "`missing_required' cleanplots"

capture findfile lcolrspace.mlib
local found = (_rc == 0)
local resolved
if `found' local resolved `"`r(fn)'"'
file write `report' "colrspace" _tab "mata_library" _tab "true" _tab ///
    "`=cond(`found', "true", "false")'" _tab `"`resolved'"' _n
if !`found' local missing_required "`missing_required' colrspace"

capture findfile lmoremata.mlib
local found = (_rc == 0)
local resolved
if `found' local resolved `"`r(fn)'"'
file write `report' "moremata" _tab "mata_library" _tab "true" _tab ///
    "`=cond(`found', "true", "false")'" _tab `"`resolved'"' _n
if !`found' local missing_required "`missing_required' moremata"

capture which gcollapse
local found = (_rc == 0)
local resolved
if `found' {
    capture findfile gcollapse.ado
    if !_rc local resolved `"`r(fn)'"'
}
file write `report' "gtools" _tab "optional_command" _tab "false" _tab ///
    "`=cond(`found', "true", "false")'" _tab `"`resolved'"' _n

file close `report'

if `"`missing_required'"' != "" {
    di as error "Missing required Stata dependencies:`missing_required'"
    exit 499
}

di as result "Required Stata dependency preflight passed."
