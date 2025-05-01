Written 03/03/2025

Key programs:

RWD_extract: gets data from rwd and neurotar, does some processing, creates traces. each session is saved seperately, but as of date of writing multiple channels from
same session are saved together
	-nt_ITImovement function to get data from neurotar setup during ITI

combineData: combines all trace/movement data (generated from RWD_extract) into one file

bstrpCombine: analyzes fiber photometrey signal using bootstrapped confidence intervals and permutation tests