# 2023_otero_et_al
R script for:<br/>
*Organ-Chips Enhance the Maturation of Human iPSC-Derived Dopamine Neurons*<br/>
<br/>
Authors:<br/>
&ensp;&ensp;Maria G. Otero<sup>1</sup>, Shaughn Bell<sup>1\*</sup>, Alexander H. Laperle<sup>1</sup>, George Lawless<sup>1</sup>, Zachary Myers<sup>1</sup>, Marian A. Castro<sup>1</sup>, Jaquelyn M. Villalba<sup>1</sup>, Clive N. Svendsen<sup>1*</sup><br/>
<br/>
Affiliations:<br/>
&ensp;&ensp;<sup>1</sup> Cedars-Sinai Board of Governors Regenerative Medicine Institute, Cedars-Sinai Medical Center, Los Angeles CA<br/>
&ensp;&ensp;<sup>\*</sup> Corresponding author, Email:&ensp;Clive.Svendsen@cshs.org<br/>
<br/>
R version 4.2.1<br/>
R Script Title:&ensp;otero_et_al_2023.R<br/>
R Script Author:&ensp;Shaughn Bell<br/>
R Script Corresponding Email:&ensp;Shaughn.Bell@cshs.org<br/>
<br/>
Notes: <br/>
&ensp;&ensp;A)&ensp;Script makes use of the variables set up under "project information" as<br/>
&ensp;&ensp;&ensp;&ensp;&ensp;well as additional "prefixes" throughout the script for ease of saving<br/>
&ensp;&ensp;&ensp;&ensp;&ensp;files with a similar path and naming structure.  When reloading data <br/>
&ensp;&ensp;&ensp;&ensp;&ensp;(see note "B"), you must either use the full path or reload the prefixes.<br/>
&ensp;&ensp;B)&ensp;Script saves intermediate steps at each major manipulation of the seurat<br/>
&ensp;&ensp;&ensp;&ensp;&ensp;object via the "saveRDS" function.  If needed, these RDS objects can then<br/>
&ensp;&ensp;&ensp;&ensp;&ensp;be reloaded to easily restart at one of these save points without needing<br/> 
&ensp;&ensp;&ensp;&ensp;&ensp;to start from scratch.  However, these are not required for analysis, and<br/>
&ensp;&ensp;&ensp;&ensp;&ensp;they can be skipped to save time and disk space.
