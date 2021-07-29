Data from preliminary limited proteolysis (LiP) experiment with DDM-purified AHA2.

50ug of protein per sample, duplicate digestions with proteinase K for 0, 30, or 60 seconds, followed by SDOC solubilization, 2% FA precipitation, and
overnight digestionwith trypsin and lys-C.

Standard solid phase extraction on C18 Omix tips. LC-MS/MS, 60 and 90 minute gradients to 30% B for each sample.

To quantify semitryptic IDs, look at each terminus of the peptide ID sequence. If the C-terminus is neither R nor K, it is semitryptic;
if the N-terminus is not preceded by R or K in the master protein sequence, it is again semitryptic. All else, tryptic.

LipCounter script (Python 3) will be developed around the 12 result files to quantify semitryptics and calculate percentages, followed by ranking of top-detected (PSM#) semitryptics and tryptics.

May automatically have it produce graphs of these results too.