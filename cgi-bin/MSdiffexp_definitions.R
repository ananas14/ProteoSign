ratios.hist.colour<-"cyan"
reps.scatter.lmline.colour<-"red"
nRequiredLeastBioreps<-2
GUI<-F
paramssetfromGUI<-F
keepEvidenceIDs<-F
exportFormat<-"pdf"
rep_order<-NA
mqValidation<-F
#working_directory<-"/var/www/html/ProteoSign/uploads//1489269061141\msdiffexp_wd"
experimental_structure_file = "exp_struct.txt"
LFQ_conditions_file = "LFQ_conditions.txt"
Rename_Array_file = "Rename_array.txt"
LS_Array_file = "LS_array.txt"
pgroups_fname<-"msdiffexp_protein.txt"
evidence_fname<-"msdiffexp_peptide.txt"
outputFigsPrefix<-"wrewerwr"
filterL_lbl<-"H"
filterL_lvl<-F
LabelFree<-F
PDdata<-F
IsobaricLabel<-F
All_MQ_Labels<-c("L", "M", "H")
time.point<-"erwrwer"
# bioreps<-3
# techreps<-1
ProteinQuantitation<-T
filterL<-F
AllowLabelRename<-T
AllowLS<-T
addLabel("L",c(""))
addLabel("H",c(""))
