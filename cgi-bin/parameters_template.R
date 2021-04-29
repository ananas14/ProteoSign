# parameters_template.R is a template that is used to feed experiment specific parameters to a PS run. The file is copied and renamed to MSDiffexp_parameter.R
# when a PS run is being prepared. all REPLACE tags are replaced with the user's preferences - notice that not all parameters are available for change

# The colour of the histogram boxes
ratios.hist.colour<-"cyan"

# The color of the linear model line
reps.scatter.lmline.colour<-"red"

# PDdata should be true if the data origin from PD and false if they origin from MQ
PDdata<-REPLACE1

# The least number of Bioreps wwhere a protein should be detected to be considered valid
nRequiredLeastBioreps<-REPLACE14

# The least number of peptides that should define a protein o that it may be considered valid
nRequiredLeastPeps<-REPLACE15

# The pThreshold for statistical analysis
pThreshold<-REPLACE16

# The desirable export format for the diagrams
exportFormat<-"pdf"
# working_directory<-"E:\xampp\htdocs\ProteoSign_v2\uploads19728424832\msdiffexp_wd"

# The name of the experimental structure file
experimental_structure_file = "exp_struct.txt"

# The name of the file that matches each MS run (raw file) to each condition
LFQ_conditions_file = "LFQ_conditions.txt"

# The text file that matches all conditions to new names (can also name two or more conditions to a common one to merge them)
Rename_Array_file = "Rename_array.txt"

# The text file that is used to perform label swap
LS_Array_file = "LS_array.txt"

# The proteingroups file (the MultiConsensus file from PD or the proteingroups file from MQ)
pgroups_fname<-"msdiffexp_protein.txt"

# The evidence file from MQ
evidence_fname<-"msdiffexp_peptide.txt"

# Replication Multiplexing information for raw files
RMrawfilesdata_file<-"RMrawfiles.txt"

# Replication multiplexing infor for tags
RMtagsdata_file<-"RMtags.txt"

# The prefix that should be added to all PS output files
outputFigsPrefix<-REPLACE2

# A boolean variable stating if there will be background species filtering
filterL<-REPLACE5

# The filter condition name for background filtering
filterL_lbl<-REPLACE7

# The filter level (peptide or protein)
filterL_lvl<-REPLACE6

# Is the experiment an LFQ one?
LabelFree<-REPLACE8

# Is the experiment an IsobaricLabel one?
IsobaricLabel<-REPLACE10

# All_MQ_labels store all the condition names that are displayed in the headers of the protein groups file of MQ
All_MQ_Labels<-REPLACE11

# exp_desc contains a description of the experiment
exp_desc<-REPLACE3

# A deprecated variable always set to true
ProteinQuantitation<-REPLACE4

# A boolean variable defining if there will be any label renaming
AllowLabelRename<-REPLACE12

# A bollean variable defining if there will be label swapping
AllowLS<-REPLACE13

# A boolean variable defining if there will be replication multiplexing
RMisused<-REPLACE17

# A boolean defining if the breps information is represented in raw files (or in tags)
RMbrepsinrawfiles<-REPLACE18

# A boolean defining if the treps information is represented in raw files (or in tags)
RMtrepsinrawfiles<-REPLACE19

# A boolean defining if the cpnditions' information is represented in raw files (or in tags)
RMconditionsinrawfiles<-REPLACE20

# A boolean defining if Functional enrichment will be done
FNenrichment<-REPLACE21

# The organism to run functional enrichment on
FNorganism<-REPLACE22
addLabel("H",c(""))
addLabel("L",c(""))
