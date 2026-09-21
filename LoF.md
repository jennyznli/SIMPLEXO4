## LOFTEE Scoring
- Evaluated in this order:  
	- level 4: 
		- AutoGVP = `likely_benign, benign`  OR
		- ClinVar = `likely_benign, benign` OR
		- not ``protein-coding`` 
	- level 1 
		- requires: 
			- ``protein-coding`` 
			- transcript in (`MANE_SELECT` OR `MANE_PLUS_CLINICAL`)
		-  keep if:
			- AutoGVP is `pathogenic, likely_pathogenic`
			- AutoGVP empty AND LOFTEE `HC`
	- level 2 
		- requires:
			- ``protein-coding`` 
			- transcript in (`MANE_SELECT` OR `MANE_PLUS_CLINICAL`
			- AutoGVP = `uncertain_significance` or `empty` 
			-  LOFTEE != `LC` 
			-   NOT (`inframe` AND  polyglutamine/polyglutamate repeat indels) 
		- meets at least one: 
			- if LOFTEE is `HC` 
			- if AlphaMissense`likely_pathogenic` 
			- if spliceAI at least one score > `0.5`
	- level 3
		- everything else, fallback


