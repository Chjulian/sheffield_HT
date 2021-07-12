## Characterising within-hospital SARS-CoV-2 transmission events: a retrospective analysis integrating epidemiological and viral genomic data from a UK tertiary care setting across two pandemic waves

This repository contains the code for the analysis of the model presented in the manuscript. The model is created and run with the script `iterate-wave(n).R`, which loads functions from the directory `./src`, read data from the directory `./data` and store the results in the directory `./output`.

### Data spec
---
The dataset should be located in `data/[mydata].csv` and contains the following fields:
```
Barcode (chr) unique identifier per person
Sex (chr) either Male (M) or Female (F)
Age (int) age at diagnosis
DateOfCollection (chr) date of sample collection, DD/MM/YY
DayOfSymptoms (int) Date of onset of symptoms, DD/MM/YY
Lineage (chr) SARS-COV-2 pango-lineage
Category (chr) either STAFF or INPATIENT 
DateOfAdmission (chr) date when the person was admitted to the hospital
HealthcareAssociation (chr) any of the 5 categories (Community Onset-Community Associated, Community Onset-Suspected Healthcare Associated, Hospital Onset-Healthcare Associated, Hospital Onset-Indetermite Healthcare Associated, Hospital Onset-Suspected Healthcare Associated) defined in file S1056 from SAGE (Contribution of nosocomial infections to the first wave) or STAFF
SwabLocation (chr) identifier of the ward where the sample was collected
Loc1...8 (chr) identifier(s) of the ward(s) occupied by the patient or visited by staff
DateOfOnset (chr) date of onset of symptoms, DD/MM/YY
```

The sequence data `.aln` should contain one sequence per person in FASTA format.

The `.csv` files with the ward occupation data contains the following fields:
```
Barcode (chr) unique identifier per person
Ward (chr) identifier of ward
In (chr) date when the person was admitted to the ward,YYYY-MM-DD HH:MM:SS
Out (chr) date when the person was discharged from the ward, YYYY-MM-DD HH:MM:SS
```
The model requires some extra arguments: 
```
scale (int, default=6): scale the occupation data in hours
reporting.probability (probability, default=0.5): the reporting probability
number.import.chains (int, default=5): number of import randomizations
incubation_period: a discrete distribution of the incubation period, can be specified via distcrete::distcrete
generation_time: a discrete distribution of the generation time, can be specified via distcrete::distcrete
```
The information is pass to outbreaker via [`create_config`](https://cran.r-project.org/web/packages/outbreaker2/outbreaker2.pdf).



