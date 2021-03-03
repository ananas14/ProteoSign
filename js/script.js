//
//
// PROTEOSIGN - script.js
// Main Javascript for Proteosign: An end-user online differential proteomics statistical analysis platform.
// Reference:
// Efstathiou G., Antonakis A. N., Theodosiou T., Pavlopoulos G. A., Divanach P., Trudgian D. C., Thomas B., Papanikolaou N., Aivaliotis M., Acuto O. and Iliopoulos I.
// https://www.ncbi.nlm.nih.gov/pubmed/28520987
//
//

//--------------------------------------
// GLOBAL VARIABLES __ START
//--------------------------------------
//{
	
	// === Configuration variables ===
	
	var cgi_bin_path = 'cgi-bin/'; // Relative path to cgi-bin folder in the server
	var RSS_prefix = "NAR Advance articles:"; // Prefix before RSS elements displayed while waiting for results
	var RSS_link = "https://academic.oup.com/rss/site_5127/3091.xml"; // RSS link to display contents while waiting for the results
	var server_Feedback_file = "/data/ProteoSign/Feedbacks.txt"; // The absolute path of the server text file that will store all feedbacks
	
	// === DEBUG flags ===
	
	var AllowMergeLabels = true;
	var AllowLS = true;
	
	// === Version variables ===
	
	var sessionid = new Date().getTime();
	var softversion = ''; // Software version as reported from the server
	
	// === Upload relative variables ===
	
	var uploaded_files = 0; // No of files that were uploaded
	// types_of_files_uploaded: Types of files that were uploaded:
	//
	// Possible values:
	// MQ are maxquant evidence files
	// MQP maxquant proteinGroups
	// PD proteome discoverer psm files
	var types_of_files_uploaded = []; 
	var nToUpload = 0; // No of files pending uploading
	
	// === Experiment type variables ===
	
	// isLabelFree and isIsobaricLabel are set to indicate if the experiment is labe; free or isobarically tagged
	// if both are off it is a labelled experiment
	var isIsobaricLabel = false;
	var isLabelFree = false;
	//procProgram refers to the program that is believed to have processed the raw data (MQ or PD) based on the uploaded files' format
	var procProgram = "";
	
	// === Label free conditions variables ===
	
	var LFQconditions = []; // Contains all assigned Label Free conditions
	var RawFileConditions = []; // Array of objects. Each object has two properties (name and condition) that match each rawfile ('name') with the corresponding assigned condition
	
	// === Renaming - merging variables ===
	
	// Label renaming is based on matching 2 or more of the labels to the same condition - this information
	// is stored in RenameArray
	// for example if an experiment contains 3 labels (L, M, H) and L and H both reflect Lung_Ca 
	// then RenameArray should be:
	//
	// L	Lung_Ca
	// M	M
	// H	Lung_Ca
	var RenameArray = [];
	// AuthenticItemsinRename contain the conditions as they are read now in conditions list
	var AuthenticItemsinRename = []; 
	var RenameFromTestData = false; // true if the test dataset chosen forces label merging
	var itemsToRename = []; // Temporary array: conditions to be renamed
	
	
	// === Label Swap variables ===
	
	// LabelSwapArray is in fact a matrix (an array of arrays) that has 4 columns: rawfile, firstlabel, secondlabel and counter (a key)
	// Each raw corresponds to a different raw file where a label swap of firstlabel and secondlabel occur
	// e.g.:
	//
	//
	// +-------------------------------+-------------+--------------+---------------+
	// |		   Raw File			| First Label | Second Label | Counter (key) |
	// +-------------------------------+-------------+--------------+---------------+
	// | OT2_Terhune(...)DMC-HLNF01_02 | L		   | H			|			 0 |
	// | OT2_Terhune(...)DMC-HLNF02_02 | L		   | H			|			 1 |
	// | OT2_Terhune(...)DMC-HLNF03_02 | L		   | H			|			 2 |
	// | OT2_Terhune(...)DMC-HLNF04_02 | L		   | H			|			 3 |
	// | (...)						 |			 |			  |			   |
	// | OT2_Terhune(...)DMC-HLNF12_02 | L		   | H			|			11 |
	// +-------------------------------+-------------+--------------+---------------+
	//
	var LabelSwapArray = [];
	var LS_array_counter = 0; // Stores the highest key that has been assigned to an entry of LabelSwapArray
	//LS_counters_per_Add - each row of this array corresponds to one assigned label swap
	// (and to one row in Label Swaps list in LS Dialog box) and has three items
	// the first is Label Swap's description as shown on the list
	// the second is an array of keys to match the label swap to Rows in LabelSwapArray
	// the last one is an array of the two labels taht are swapped in this label swap
	//e.g.
	//
	// +---------------+-------+--------+
	// | Description   | Keys  | Labels |
	// +---------------+-------+--------+
	// | L - H in b1t2 | 1	   | L	    |
	// |			   +-------+--------+
	// |			   | 2	   | H	    |
	// |			   +-------+--------+
	// |			   | 3	   |		|
	// |			   +-------+		|
	// |			   | 4	   |		|
	// |			   +-------+		|
	// |			   | 5	   |		|
	// |			   +-------+		|
	// |			   | (...) |		|
	// |			   +-------+		|
	// |			   | 11	   |		|
	// +---------------+-------+--------+
	//
	var LS_counters_per_Add = [];
	var LSselected_raws = []; // An array that stores all raw files that will be affected by a selected label swap
	
	// === DOM elements ===
	
	var rawfiles_tbl_allfiles_DT;
	var main_context_menu;
	var label_context_menu;
	
	// === Turning gear variables ===
	
	var sectionSpinnerOn = false; // Is the turning gear shown?
	var sectionSpinnerText = ""; // What text lies next to the turning gear?
	
	// === Experimental structure variables ===
	
	var bioreps; // No of bioreps in current experimental structure
	var techreps; // No of techreps in current experimental structure
	var fractions; // No of fractions in current experimental structure
	var rawfiles; // No of rawfiles in current experimental structure
	// rawfiles_structure is an array of objects having the following structure:
	// {rawfile: ..., biorep: ..., techrep: ..., fraction: ..., used: ...}
	// as an example rawfiles_structure might seem like
	// 0: {rawfile: "D131106_014.raw", biorep: 1, techrep: 1, fraction: 1, used: true}
	// 1: {rawfile: "D131106_016.raw", biorep: 2, techrep: 1, fraction: 1, used: true}
	// 2: {rawfile: "D131106_017.raw", biorep: 3, techrep: 1, fraction: 1, used: true}
	var rawfiles_structure;
	
	// === RSS variables ===
	var rss_i; // Currently displayed RSS element's index
	
	// === Assistant global variables ===
	
	// AddedLabels is a flag to make sure that two different files will not both
	// attempt to populate the conditions list
	var AddedLabels = false;
	var analysis_finished = false; // True if the analysis was successfully completed by at least one server procedure
	var datatestOK_clicked = false; // True if the user chose a dataset from the corresponding dialog box
	var my_lbls_toselect = []; // Valid if test dataset was chosen: labels to select for comparison as defined by test dataset's parameter
	// Since MQ data analysis requires MQ labels to be sent to the backend in the order they appear in file's headers
	// a special var called my_all_mq_labels will be sent in such a case containing an R-structured array [c(... , ..., ...)]
	// to be sent to the backend
	var my_all_mq_labels = "";
	var peptideLabelsNamesFromFile = []; // The conditions as retrieved from a standard label experiment
	var rep_counts; // Assisting variable to auto-fill some experimental structure table entries - see set_reps for more information
	var lastclicked_rawfiles_tbl_tr = null; // Last clicked element on exp structure table
	var lastclicked_RMrawfiles_tbl_tr = null; // Last clicked element on RM table
	var tempAllowMergeLabels = true;
	var tempAllowLS = true;
	var RreturnedErrorsInCurSession = false; // a temporary variable to set if R returned errors in the current R session
	var log_test_dataset = false;
	//The following var is used to ensure that the feedback button will be animated just once
	var feedbackAnimated = false;
	var sessions_that_run = []; // array of the sessionids that run
	// === Advanced options variables ===
	
	var LeastBRep = 2;// least breps where a protein must be found in order not to be disqualified
	var LeastPep = 2; //least peptides where a protein must be found in order not to be disqualified
	var Pthreshold = 0.05; //max p where a protein is considered as differentially expressed
	var my_step = 0; // The step from which parameters are loaded
	
	// === Star-based rating variables ===
	
	var star_hovered = 1; // The last star the user hovered on
	var star_chosen = 1; // Amount of stars the user chose
	var stars_enabled = true; // If star-based rating is enabled
	
	// === Replication multiplexing (RM) variables ===
	
	var RM_RMstep = 1; // The current step in replication multiplexing
	
	// RMrawfilesdata contains all the information for the rawfiles in case Replication Multiplexing is chosen
	// RMrawfilesdata is a two dimension array that has 8 columns:
	// counter, rawfilename, brep, trep, frac, cond, used, selected
	// An example is:
	//
	// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
	// | counter | rawfilename                                               | brep | trep | frac | cond | used | selected |
	// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
	// | 1       | M456-E04-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
	// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
	// | 2       | M456-D10-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
	// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
	// | 3       | M456-E08-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
	// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
	// | 4       | M456-C09-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
	// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
	// | 5       | M456-D05-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
	// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
	//
	var RMrawfilesdata = []; 
	// RMtagsdata contains all the information for the tags in case Replication Multiplexing is chosen, the structure is similar
	// to RMrawfilesdata
	var RMtagsdata = []; 
	// In RM breps treps  and conditions bight be represented as different rawfiles or tags
	// the following three boolean variables if these three attributes are represented as different rawfiles or tags
	var RMbrepsRepInRawFiles = true;
	var RMtrepsRepInRawFiles = true;
	var RMconditionsRepInRawFiles = true;
	var RMtags = []; // The tags of the experiment
	var RMisused = false; //this will tell R if Replication Multiplexing will be used in the current session of PS
	var stepRM2initialized = false;
	var stepRM3initialized = false;
	
	
//}
//--------------------------------------
// GLOBAL VARIABLES __ END
//--------------------------------------



// MAIN ROUTINES __ START
//----------------------------------------
//{
	
	var postFile = function (idx, file, serverSide, postSuccess) {
		
		// postFile handles the asynchronous uploading of a single file calling the upload_files.php on the server
		// postFile sends the sessionid and the filename to the server
		// it also sends server-side that is set to false if the user uploads a file from client
		// and true if the files are simply copied from a test dataset of the server to the current work folder
		//	
		// upload_files.php transfers the data to the current work folder and detects
		// the file type (MQ for MQ_evidence files, MQP for MQ_proteinGroups files and PD for PD psm files) and
		// the experiment type (Metabolically labelled (L), Isobarically labelled (IL) or Label Free (LF))
		// Also, it detects all labels/tags that may be contained in the file and
		// all raw files.
		// The variables that store the aforementioned information are returned in file_type, exp_type, peptide_labels_names and raw_filesnames respectively
		
		var thedata = new FormData();
		thedata.append('thefile', file);
		thedata.append('session_id', sessionid);
		thedata.append('server_side', serverSide);
		
		// progresstditm (progress td item) corresponds to the table item that shows the progress bar of the
		// current uploading process, the following XMLHttpRequest adds a listener so that
		// the progress bar is refreshed regularly.
		
		var progresstditm = $("#s2uluploaders table tr:nth-child(" + (idx + 1) + ") td").get(1);
		
		$.ajax({
			url: cgi_bin_path + 'upload_files.php',
			type: 'POST',
			xhr: function () {  // Custom XMLHttpRequest
				var myXhr = $.ajaxSettings.xhr();
				if (myXhr.upload) { // Check if upload property exists
					myXhr.upload.addEventListener('progress', function (e) {
						if (e.lengthComputable) {
							if (e.loaded / e.total == 1.0) {
								// when the whole file gets uploaded replace the progressbar with a message 
								$(progresstditm).html("<span class='uploadSuccessMsg'><strong><em>Processing, please wait ...<em><strong></span>");
								// display the turning gear next to heading #2 to show that an asynchronous process is running
								toggleCurrentSectionSpinner();
							}
							else
							{
								// update the progress bar value
								helper_setItemAttr("#uploadfile" + idx, {value: e.loaded, max: e.total});
							}
						}
					}, false); // END addEventListener - XHR
				}
				return myXhr;
			},
			// Form data
			data: thedata,
			//Options to tell jQuery not to process data or worry about content-type.
			cache: false,
			contentType: false,
			processData: false,
			//Ajax events
			beforeSend: function (jqXHR, settings) {
				
				// before starting the AJAX request make sure to set the progress bar to 0
				if (serverSide) {
					helper_setItemAttr("#uploadfile" + idx, {value: 0, max: 100});
				}
			}}).done(function (data, textStatus, jqXHR) {
			
			//In case the sessionid changed in session1 do not use the file in ProteoSign
			if (data.ret_session != sessionid)
			{
				return;
			}
			
			// Immediately after uploading a file successfully:
			// turn the turning gear off if this was the last file to be uploaded
			
			if(sectionSpinnerOn & nToUpload == 1) toggleCurrentSectionSpinner();
			
			// Refresh the nToUpload and uploadedFiles integers showing how many files are to be uploaded and how many have been uploaded
			nToUpload = nToUpload - 1;
			uploaded_files++;		
			
			// If everything went fine enable the button for next stage and print OK to the progress item.
			// Print the error message otherwise.
			
			if (data.success) {
				// types_of_files_uploaded is an array that contains all types of files that are uploaded without
				// matching them with the corresponding file
				
				types_of_files_uploaded.push(data.file_type);
				if (data.msg == "")
				{
					// if php did not return a message this states the procedure's success
					$(progresstditm).html("<span class='uploadSuccessMsg'><strong><em>OK<em><strong></span>");
				}
				else
				{
					// Print the warning message in case there is one (e.g. the same file has been uploaded etc.)
					$(progresstditm).html("<span class='uploadWarningMsg'><strong><em>" + data.msg + "<em><strong></span>");
				}
				
				// upload_files.php also detects the RawFiles that were used by MQ or PD to create the
				// uploaded multiconsensus files (e.g. "D131106_014.raw", "D131106_016.raw", "D131106_018.raw"...)
				
				
				// If we have info regarding rawdata file names, save it to global array rawfiles				
				if (data.raw_filesnames.length > 0)
				{
					rawfiles = data.raw_filesnames;
					//Reset the replicate data so that they can be safely defined by the user in Stage 2
					reset_reps();
				}
				
				// Now get the experiment type (MQ and PD are the only filetypes that will return the exp type correctly)
				
				if (data.file_type == "MQ" || data.file_type == "PD")
				{
					switch(data.exp_type)
					{
						case "L":
							isIsobaricLabel = false;
							isLabelFree = false;
							break;
						case "IL":
							isIsobaricLabel = true;
							isLabelFree = false;
							break;
						case "LF":
							isIsobaricLabel = false;
							isLabelFree = true;
							break;
					}
				}
				
				// The last task upload_files.php is assigned to is to find all labels (Light, Heavy etc.) or tags (126, 127 etc.)
				// that are contained in labeled experiments
				// These are stored in data.peptide_labels_names
				// If no labels were detected in previous download and if the current filetype is expected to give us info about labels
				// fetch the labels from data.peptide_labels_names
				
				if (peptideLabelsNamesFromFile.length == 0 && (data.file_type == "MQ" || data.file_type == "PD"))
				{
					
					peptideLabelsNamesFromFile = data.peptide_labels_names.sort();
					
					// some HTML items' visibility depends on whether the experiment is label-free or not get a jquery pointer for them:
					
					// 1. Non-labelled "background" species present? (mainly for pSILAC experiments)
					var item1 = $("#s3advparams input[name='explbl00']").closest("tr").children().first();
					// 2. Quantitation filtering?
					var item2 = $("#s3advparams input[name='expquantfilt']").closest("tr").children().first();
					
					
					// For metabolically and isobarically labelled data:
					if (!isLabelFree)
					{
						
						// Set states of parameters relevant to labelled experiments accordingly
						$(item1).prop('disabled', false);
						$(item2).prop('disabled', false);
						
						// Since this is detected to be a labelled experiment hide the assign condition from right click context menu
						// on the table of Stage 2 - this is only helpful in Label Free experiments
						
						$(document).contextmenu("showEntry", "assign_condition", false);
						
						// Hide a tip for assigning conditions from the user since it is not relevant in labelled experiments
						
						$("#ExtraInfoForLFQ").hide();
						
						// Hide the last column in Stage 2 table that corresponds to assigned conditions
						
						rawfiles_tbl_allfiles_DT.column("4").visible(false);

						// Quantitation filtering is applicable only in labelled and isobarically tagged experiments
						$("#s3QuantitationFiltering").show();
						
						// The conditions list in Stage 3 should contain all conditions that can be cross compared
						// The following block populates this list
						//
						// Fill the conditions list
						//
						// Label merging or Label Renaming is a feature where the user is able to
						// match 2 or more labels to the same condition. See our documentation on "Label Merging"
						// for more information
						//
						// if merge labels is disabled then simply populate the conditions list with the conditions:
						// AddedLabels is a flag to make sure that two different files will not both attempt to populate
						// the conditions list
						
						if(!AddedLabels)
						{
							// First erase everything from the lists
						
							$('#conditions_list').empty();
							
							// The folowing block populates the conditions list and the dropdown menu expquantfiltlbl
							// used to define as background label a label in pSILAC experiments
							
							$("#s3advparams select[name='expquantfiltlbl']").empty();
							
							// If merging is allowed populate conditions list and initialize merging-renaming
							// In case a test dataset was loaded that contained label merging 
							// labels are merged automatically while loading the test dataset
							
							if(!AllowMergeLabels)
							{
								// Get all loaded labels
								$.each(peptideLabelsNamesFromFile, function (idx, lblname_i) {
									// If the label is number parse it into a string
									if (StringisNumber(lblname_i))
									{
										peptideLabelsNamesFromFile[idx] = parseInt(lblname_i);
										lblname_i = parseInt(lblname_i);
									}
									// Populate the lists
									$("#conditions_list").append("<option value='" + lblname_i + "' style='padding: 2px; font-size: 125%;' selected='true'>" + lblname_i + "</option>");
									$("#s3advparams select[name='expquantfiltlbl']").append("<option oncontextmenu='return false;' value='" + lblname_i + "'>" + lblname_i + "</option>");
								});// END for each peptideLabelsNamesFromFile
							}
							else
							{
								if (!RenameFromTestData)
								{
									// Get loaded labels
									$.each(peptideLabelsNamesFromFile, function (idx, lblname_i) {
										// If the label is number parse it into string
										if (StringisNumber(lblname_i))
										{
											peptideLabelsNamesFromFile[idx] = parseInt(lblname_i);
											lblname_i = parseInt(lblname_i);
										}
										// Populate the lists
										$("#conditions_list").append("<option value='" + lblname_i + "' style='padding: 2px; font-size: 125%;' selected='true'>" + lblname_i + "</option>");
										$("#s3advparams select[name='expquantfiltlbl']").append("<option oncontextmenu='return false;' value='" + lblname_i + "'>" + lblname_i + "</option>");
									}); // END for each peptideLabelsNamesFromFile
									InitializeRename();
								}
							}
							// Select labels to compare according to test datasets options (if the dataest is not a test dataset
							// select_labels_according_to_test_dataset will simply do nothing)
							select_labels_according_to_test_dataset();
							
							// Sort MQ labels as needed
							create_my_all_mq_labels();
							
							AddedLabels = true;
						}
				
						// INFO: Add here more labeled data - specific functionality
						
					}
					
					if (isLabelFree)
					{
						
						// Label-free case (disable non-applicable parameters and rename others accordingly)
						// see above comments for more information
						$(item1).prop('disabled', true);
						$(item2).prop('disabled', true);
						$("#ExtraInfoForLFQ").show();
						$(document).contextmenu("showEntry", "assign_condition", true);
						rawfiles_tbl_allfiles_DT.column("4").visible(true);
						$("#s3QuantitationFiltering").hide();
						
						//INFO: Add here more functionality in case of LFQ data
						
					}
					
					if (isIsobaricLabel)
					{
						
						// Enable Replication Multiplexing
						$(document).contextmenu("showEntry", "rep_mult", true);
						
						//INFO: Add here more functionality in case of Isobarically labelled data
						
					}
					else
					{
						$(document).contextmenu("showEntry", "rep_mult", false);
					}
				}
				
			}
			else // if the php did not return success
			{
				$(progresstditm).html("<span class='uploadErrorMsg'><strong><em>A server-side error occurred: " + data.msg + "<em><strong></span>");
			}
			if (typeof (postSuccess) == "function")
			{
				// Run the postSuccess function (the postsuccess function runs if the asynchronous process was completed
				// with no unexpected issues)
				
				postSuccess();
			}
			}).fail(function (jqXHR, textStatus, errorThrown) {
			
			// if there were unexpected issues with asynchronous ajax running
			if(sectionSpinnerOn) toggleCurrentSectionSpinner();
			$(progresstditm).empty();
			$(progresstditm).html("<span class='uploadErrorMsg'><strong><em>An AJAX error occurred: " + errorThrown + "<em><strong></span>");
			nToUpload = nToUpload - 1;
			uploaded_files++;
		});
	}
	
	
	var postFireUpAnalysisAndWait = function () {
		
		// postFireUpAnalysisAndWait initiates the analysis server-sided in an asynchronous manner
		// perform_analysis.php is responsible to setup and run the R script on the server
		
		// Add the current sessions to sessions that run
		if (sessions_that_run.indexOf(sessionid) == -1)
		{
			sessions_that_run.push(sessionid);
		}
		else
		{
			return;
		}  
		// Before initiating the analysis create a save parameters file to the output folder
		SaveParams(1);
		// Also make sure that the analysis_finished var is false since we are initiating a new analysis
		analysis_finished = false;
		RreturnedErrorsInCurSession = false;
		
		var thedata = new FormData();
		thedata.append('session_id', sessionid);
		thedata.append("AllowMergeLabels", AllowMergeLabels ? "T" : "F");
		thedata.append("AllowLS", AllowLS ? "T" : "F");
		
		//get the processing program (MQ or PD)
		var quantlabel = $("#quantitation_prog_lbl").text();
		var pattern = new RegExp("MaxQuant");
		var res = pattern.test(quantlabel);
		var provenprocprogram = "";
		if (res == true)
		{
			provenprocprogram = "MQ";
		}
		else
		{
			provenprocprogram = "PD";
		}
		
		//Append this info to thedata
		thedata.append('proc_program', provenprocprogram);
		
		
		if(!sectionSpinnerOn)
		{
			toggleCurrentSectionSpinner();
		}
		$("#s4btnf").prop('disabled', true);
		$.ajax({
			url: cgi_bin_path + 'perform_analysis.php', //Server script to fire up data analysis
			type: 'POST',
			// Form data
			data: thedata,
			//Options to tell jQuery not to worry about content-type.
			processData: false,
			cache: false,
			contentType: false,
			beforeSend: function (jqXHR, settings) {
				
				// get RSS data to display while waiting for the results
				getRSS(RSS_link, "#server_feedback");
			},
			}).done(function (data, textStatus, jqXHR) {
			
			// When the analysis is finished display the results back to the user
			
			// First make sure that analysis_finished flag is not set (this may occur if erroneously more than
			// one backend procedures ran at the same time)
			//
			// Also make sure that data.ret_session (the session id reported back by server's backend) is the
			// same with the curren sessionid. a session id mismatch may happen if the user opted to hit the Back button in stage 4
			// that would fire a session change
			
			if (data.ret_session != sessionid)
			{
				// Abort results
				return;
			}
			
			// Remove the current session from sessions that run
			if (sessions_that_run.indexOf(sessionid) != -1)
			{
				sessions_that_run.splice(sessions_that_run.indexOf(sessionid), 1);
			}
			
			var R_returned_error = false;
			
			// Hide the spinning gear
			if(sectionSpinnerOn)
			{
				toggleCurrentSectionSpinner();
			}
			
			// Results will be displayed on server_feedback
			// First empty the division:
			$("#server_feedback").empty();
			$("#s4btnf").prop('disabled', !data.success);
			
			// Backend may have reported information warnings or errors as messages back to the user
			// This is sent via data.UserInfo that contain each message per line. If it begins with
			// Warn User: it is a warning, Error User: an error and Info User simple information
			var my_UserInfo = data.UserInfo;
			var UserInfoDisplay = "";
			if(typeof(my_UserInfo) != "undefined" && my_UserInfo != "")
			{
				var lines = my_UserInfo.split("\t\t");
				$.each(lines, function(idx, my_line){
					UserInfoDisplay = UserInfoDisplay + my_line.replace(/\[.*?\]/i, "</p>&#13");
				});
			}
			if (UserInfoDisplay != "")
			{
				// Update the text in user info dialog box and show it
				UserInfoDisplay = UserInfoDisplay.replace(/Warn User: /g, '<p style="color: #BA4A00">')
				UserInfoDisplay = UserInfoDisplay.replace(/Error User: /g, '<p style="color: #E60000">')
				UserInfoDisplay = UserInfoDisplay.replace(/Info User: /g, '<p>')
				$("#s4UserInfoText").empty();
				$("#s4UserInfoText").append(UserInfoDisplay);
				$(".expparamsDlg").css({"left": ($("body").width() / 2) - ($("#s4UserInfo").width() / 2)});
				$('body').append('<div id="mask"></div>');
				$("#s4UserInfo").fadeIn(300);
				$('#mask').fadeIn(300);
			}
			
			if (data.success)
			{
				// Now display the result images to the user and provide him with
				// a download link for the results in zip format
				//
				// WARNING: the following lines consider a "/" file seperator. Might not work on Windows based systems
				$("#results_p").html("Now you can inspect your results. When ready, click <em>Next</em>.");
				$("#dndres").html("<span><a href=" + data.results_url + "><strong>" + data.results_url.substr(data.results_url.lastIndexOf("/") + 1) + "</strong></a></span>");
				patt = new RegExp("limma\-graphs");
				
				// Populate server_feedback division with thumbnails of the result images, data.results_preview store their server paths
				$.each(data.results_preview, function (idx, path_to_img_i)
					{
						var img_i = path_to_img_i.substr(data.results_url.lastIndexOf("/") + 1);
						if (!patt.test(img_i)) {
							$("#server_feedback").append("<div class='resimg'><a href='" + path_to_img_i + "' target='_blank'><img src='" + path_to_img_i + "' width='120'></img></a></div>");
						}
					}); // END for each data.results_preview
					$("#server_feedback").css("box-shadow", "0 4px 8px 0 rgba(0, 0, 0, 0.2), 0 6px 20px 0 rgba(0, 0, 0, 0.19)")
			}
			else
			{
				// Otherwise, display an error message
				$("#results_p").html("");
				$("#server_feedback").html("<span class='uploadErrorMsg'><strong><em>The analysis could not be completed: " + data.msg + "<em><strong></span>");
				if (data.R_dump.length > 0)
				{
					$("#server_feedback").append("<br><br><span style='font-family: Georgia; font-size: 95%; text-align: left; display: inline-block; width: 90%'><p>Please ensure that input parameters (such as number of replicates, number of biological conditions/labels etc) are correctly defined and input data format is valid. The statistical analysis routine relies heavily on the validity of the input parameters.</p><p>If the above does not apply, then the statistical analysis may have failed due to numerical problems (e.g. there were too many missing values/data points).</p><p> If you feel that none of the above is the case, please click <a href='mailto:msdiffexp@gmail.com,iliopj@med.uoc.gr?Subject=Session%20" + sessionid + "' target='_blank'><u>here</u></a> to notify via e-mail (do not delete the session id in the subject) the ProteoSign team for investigation of your analysis issue.</p></span>")
				}
				R_returned_error = true;
				RreturnedErrorsInCurSession = true;
				$("#s4btnb").prop("disabled", false);
			}
			if (!R_returned_error) analysis_finished = true;
			}).fail(function (jqXHR, textStatus, errorThrown) {
			if(sectionSpinnerOn)
			{
				toggleCurrentSectionSpinner();
			}
			$("#server_feedback").empty();
			$("#server_feedback").html("<span class='uploadErrorMsg'><strong><em>An AJAX error occurred: " + errorThrown + "<em><strong></span>");
		});
	}
	
	
	var postParameters = function (params) {
		
		// postParameters is the first function executed when running the analysis (immediately after)
		// validating the input
		//
		// postParameters executes upload_parameters.php and sends the parameters of the analysis
		// Some of them are the values of several DOM elements, the argument params is a collection of these elements
		// Their values will be parsed to upload_parameters.php. The rest of the parameters
		// are global variables
		
		var thedata = new FormData();
		
		// Send the session_id
		thedata.append('session_id', sessionid);
		
		// For each DOM element store its value to thedata
		$.each(params, function (idx, param_i)
			{
				var theval = $(param_i).val();
				switch ($(param_i).attr('type'))
				{
					case "checkbox":
					// If the element is a checkbox
					theval = ($(param_i).prop("checked") ? 1 : 0);
					// Here on/off is transformed to T/F so that it easily read by R
					thedata.append($(param_i).attr('name'), (theval ? "T" : "F"));
					break;
					default:
					if (theval == null)
					{
						theval = '';
					}
					if ($(param_i).attr('name') == "conditions_list")
					{
						// Especially for conditions list get only the selected items and map
						// each one to an anoynymous function in order to create the conditions
						var mytmp = "#" + $(param_i).attr('id') +  " option:selected";
						var tmp_i = 1;
						$(mytmp).map(function ()
							{
								if (!AllowMergeLabels)
								{
									// If label merging is disabled simply append the selected items text to thedata
									thedata.append("explbl" + (tmp_i) + "name", $(this).text());
									tmp_i++;	
								}
								else
								{
									var my_temp_text = $(this).text(); // Selected conditions_list item's text
									// If label merging is enabled append the selected items text to thedata if the condition is authentic
									// otherwise the authentic ones that correspond to the merged one
									if($.inArray(my_temp_text, AuthenticItemsinRename) != -1)
									{
										thedata.append("explbl" + (tmp_i) + "name", my_temp_text);
										tmp_i++;
									}
									else
									{
										$.each(RenameArray, function(idx, my_row)
											{
												if(my_row[1] == my_temp_text)
												{
													// If we found an entry in RenameArray that has my_temp_text as merged condition...
													// append its corresponding authentic condition to thedata
													// see RenameArray's structure for more information
													thedata.append("explbl" + (tmp_i) + "name", my_row[0]);
													tmp_i++;
												}
											}); // END for each RenameArray					
									}
								}
							});
					}
					else
					{
						// Most of the elements hold their value in $(param_i).val() and need no further transformation
						thedata.append($(param_i).attr('name'), theval);
					}
					break;
				}
			}); // END for each parameter
			
			
			// More parameters to append:
			
			thedata.append("labelfree", isLabelFree ? 'T' : 'F');
			thedata.append("AllowMergeLabels", AllowMergeLabels ? "T" : "F");
			thedata.append("AllowLS", AllowLS ? "T" : "F");
			if (RMisused)
			{
				// If Replication Multiplexing is used get a copy of rawfiles structure and 
				// push 2 example entries so that gen_expdesign will send them to the backend
				var temp_rf_structure = rawfiles_structure.slice();
				
				rawfiles_structure.push({rawfile: "RM_exmpl1", biorep: "1", techrep: "1", fraction: "1", used: "1"});
				rawfiles_structure.push({rawfile: "RM_exmpl2", biorep: "2", techrep: "1", fraction: "1", used: "1"});
				thedata.append("exp_struct", gen_expdesign(rawfiles_structure));
				
				// Restore rawfile_structure to its original state
				rawfiles_structure = temp_rf_structure.slice();
			}
			else
			{
				thedata.append("exp_struct", gen_expdesign(rawfiles_structure));
			}
			thedata.append("LFQ_conds", gen_lfqdesign(RawFileConditions));
			thedata.append("LeastBreps", LeastBRep);
			thedata.append("LeastPeps", LeastPep);
			thedata.append("PThreshold", Pthreshold);
			if(AllowMergeLabels) thedata.append("Rename_Array", gen_RenameFile(RenameArray));
			if(AllowLS) thedata.append("LabelSwapArray", gen_LSArray(LabelSwapArray));
			
			// Deal with the RM vars:
			thedata.append("RMRawFiles", generate_tab_file(RMrawfilesdata));
			thedata.append("RMTags", generate_tab_file(RMtagsdata));
			thedata.append("RMbrepsRepInRawFiles", RMbrepsRepInRawFiles ? "T" : "F");
			thedata.append("RMtrepsRepInRawFiles", RMtrepsRepInRawFiles ? "T" : "F");
			thedata.append("RMconditionsRepInRawFiles", RMconditionsRepInRawFiles ? "T" : "F");
			thedata.append("RMisused", RMisused ? "T" : "F");
			
			thedata.append("IsIsobaricLabel", isIsobaricLabel ? "T" : "F");
			thedata.append("All_MQ_Labels", my_all_mq_labels);
			
			// upload_parameters.php is responsible to set up the parameters so that R backend can work flawlessly
			// it produces a file called MSdiffexp_definitions.R that is executed as a script and imports most of the
			// parameters in R. It also creates some files (i.e. exp_struct.txt LFQ_conditions.txt LS_array.txt and Rename_array.txt)
			// that contain relevant information in tabular format and are then parsed by R
			
			$.ajax({
				url: cgi_bin_path + 'upload_parameters.php', //Server script to receive parameters
				type: 'POST',
				// Form data
				data: thedata,
				//Options to tell jQuery not to worry about content-type.
				processData: false,
				cache: false,
				contentType: false
				}).done(function (data, textStatus, jqXHR) {
				//if there was a server-side error alert.
				if (!data.success) {
					msgbox("ERROR on SERVER: " + data.msg);
				}
				else
				{
					// If everything went fine Fire Up the analysis
					postFireUpAnalysisAndWait();
				}
			});
	}
	
	
	var postTestDatasetInfo = function (dataset_desc) {
		
		// When the user selects to display a test dataset, postTestDatasetInfo
		// executes load_test_data.php to retrieve all necessary info to display the test dataset
		// load_test_data.php returns information regarding the exp structure, the conditions and
		// the analysis options retrieved from an SQLITE3 database called testdatadb stored in test data folder
		//
		// For more information about testdatadb and how to edit it see our documentation
		//
		// Data returned by load_test_data.php contain the following arrays:
		// file: an array of the names of the files in the server to be loaded (e.g. [MQ_SILAC_1909_2plex_evidence.txt, MQ_SILAC_1909_2plex_proteinGroups.txt]) 
		//	  they are stored under test data folder
		// raw_file: the names of the rawfiles from top to bottom as displayed in table in stage 2
		// brep: the bio reps from top to bottom as displayed in the table in stage 2
		// trep: the tech reps from top to bottom as displayed in the table in stage 2
		// frac: the fractions from top to bottom as displayed in the table in stage 2
		// used: 0 for false and 1 for true, information about using the corresponding raw file in the analysis as displayed in the tab le in stage 2 from top to bottom
		// condition: for Label free data the condition of the corresponding raw file as displayed in table in stage 2 from top to bottom
		// selector: a jQuery selector
		// value: the value to be set to the corresponding jQuery element
		// 
		// postTestDatasetInfo will set value [0] to item corresponding to selector[0] value [1] to $(selector[1]) and so on
		//
		// option: there might be special options returning from the dataset. If option is not undefined it denotes that a special
		//	  procedure must be executed. Its argument is denoted bu opt_value, another variable returned in this case, retrieved by the database
		//	  Typically opt_value is a one-line encoded array that is parsed to a global array of PS
		// 
		// | Possible options |		  Description		       |
		// |------------------|--------------------------------|
		// | Rename		      | Renames (merges) conditions	   |
		// |				  | by editing the RenameArray	   |
		// | LS_Array		  | Sets Label Swap by editing the |
		// |				  | LS_Array					   |
		// | LS_c_p_Add	      | When a test dataset has	       |
		// |				  | Label Swaping LS_c_p_Add	   |
		// |				  | defines the structure of	   |
		// |				  | LS_counters_per_Add			   |
		// | Select_Labels	  | Selects a subset of the		   |
		// |				  | available conditions		   |
		//
		// opt_value: the argument of the corresponding option, usually denoting an array, opt_value is transformed to
		// an array by breaking the first dimension every time a | is found, the second dimension every time a || is
		// found etc.
		//
		// e.g. the value 0|Lung_SCC||1|Lung_SCC||2|Lung_SCC||3|Lung_ADC||4|Lung_ADC||5|Lung_ADC will become this array:
		//
		// +---+----------+
		// | 0 | Lung_SCC |
		// +---+----------+
		// | 1 | Lung_SCC |
		// +---+----------+
		// | 2 | Lung_SCC |
		// +---+----------+
		// | 3 | Lung_ADC |
		// +---+----------+
		// | 4 | Lung_ADC |
		// +---+----------+
		// | 5 | Lung_ADC |
		// +---+----------+
		//
		
		
		// Send the necessary info to PHP
		var thedata = new FormData();
		thedata.append('session_id', sessionid);
		thedata.append('descriptions_requested', false);
		thedata.append('dataset_info_requested', dataset_desc);
		
		$.ajax({
			url: cgi_bin_path + 'load_test_data.php',
			type: 'POST',
			// Form data
			data: thedata,
			//Options to tell jQuery not to worry about content-type.
			processData: false,
			cache: false,
			contentType: false
			}).done(function (data, textStatus, jqXHR) {
			// If there was a server-side error alert.
			if (!data.success)
			{
				msgbox("ERROR on SERVER: " + data.msg);
			}
			else
			{
				uploadingFiles = data.queryres.file;
				if (uploadingFiles && uploadingFiles.length > 0)
				{
					
					// Start uploading files (in this case we simply copy the files from test data to
					// our work folder server side)
					//
					// Setting serverside argument of uploadingFiles to true, uploadFiles passes the argument (via postFile)
					// to upload_files.php so the server knows to simply copy the test data file to the workfolder
					// on the server
					
					uploadFiles(true, function ()
						{
							// This anonymous function is the postSuccess function
							//
							// After "pseudo uploading" the test files enable next button in Stage 1 and
							// trigger it to advance automatically to stage 2
							
							if (nToUpload > 0)  return;
							
							$("#s2btnb").prop('disabled', false);
							$("#dlgTestDatasetsBtnOK").prop('disabled', false);
							datatestOK_clicked = false;
							$("#s2btnf").triggerHandler("click");
							
							// Replace any illegal characters in experiment description just in case
							
							$("input[name='expid']").val(dataset_desc.replace(/[^a-zA-Z0-9]+/g, "_"));
							
							
							$.each(data.queryres.selector, function (idx, param_selector)
								{
									// For each item denoted by [selector] set its corresponding attribute to [value]
									switch ($(param_selector).attr('type')) {
										case "checkbox":
											// Here the 0/1 means false/true
											var theval = (data.queryres.value[idx] == "0" ? false : true);
											if ($(param_selector).prop("checked") != theval) {
												$(param_selector).prop("checked", theval);
											}
											break;
										default:
											$(param_selector).val(data.queryres.value[idx]);
											break;
									}
								}); // END for each data.queryres.selector
								
								
								if ($("input[name='exppddata']").is(':checked'))
								{
									procProgram = "PD";
									$("#quantitation_prog_lbl").text("\u2014 Raw data were quantified with Proteome Discoverer\u2122");
								}
								else
								{
									procProgram = "MQ";
									$("#quantitation_prog_lbl").text("\u2014 Raw data were quantified with MaxQuant");
								}
								on_exppddatachange();
								for (var i = 0; i < data.queryres.raw_file.length; i++)
								{
									// Get the rawfiles information from the database, refresh rawfile_structure
									// and then refresh the table on stage 2
									
									rawfiles_structure.push({rawfile: data.queryres.raw_file[i], biorep: data.queryres.brep[i], techrep: data.queryres.trep[i], fraction: data.queryres.frac[i], used: data.queryres.used[i]});
									if (data.queryres.condition[i] != "-")
									{
										RawFileConditions.push({name: data.queryres.raw_file[i], condition: data.queryres.condition[i]});
									}
									
									var tds = $('tr[id="tr_' + data.queryres.raw_file[i] + '"] td');
									$(tds[1]).text(data.queryres.brep[i] == 0 ? '-' : data.queryres.brep[i]);
									$(tds[2]).text(data.queryres.trep[i] == 0 ? '-' : data.queryres.trep[i]);
									$(tds[3]).text(data.queryres.frac[i] == 0 ? '-' : data.queryres.frac[i]);
									if (data.queryres.condition[i] != "-")
									{
										$(tds[4]).text(data.queryres.condition[i]);
									}
									if (data.queryres.used[i] == 0)
									{
										$(tds[0]).css("text-decoration", "line-through");
										$(tds[1]).css("text-decoration", "line-through");
										$(tds[2]).css("text-decoration", "line-through");
										$(tds[3]).css("text-decoration", "line-through");
										if (isLabelFree == true) $(tds[4]).css("text-decoration", "line-through");
									}
								}
								
								if (typeof data.queryres.option !== "undefined")
								{
									for (var i = 0; i< data.queryres.option.length; i++)
									{
										switch (data.queryres.option[i])
										{
											//For the extra options:
											case "Rename":
											//The following automates the load of a test dataset that contains merged labels:
											if(!AllowMergeLabels) break;
											RenameArray = [];
											AuthenticItemsinRename = [];
											
											// Transform one line encoded array to an array using the rules detailed above
											var my_val = data.queryres.opt_value[i];
											var lines = my_val.split("||");
											$.each(lines, function(idx, my_line){
												var my_values = my_line.split("|");
												RenameArray.push(my_values);
												AuthenticItemsinRename.push(my_values[0]);
											});
											RenameFromTestData = true;
											Refresh_conds_list();
											break;
											case "LS_Array":
											//The following block automates the loading of a dataset that contains label swap, it builds the LabelSwapArray
											if(!AllowLS) break;
											LabelSwapArray = [];
											// Transform one line encoded array to an array using the rules detailed above
											var my_val = data.queryres.opt_value[i];
											var lines = my_val.split("||");
											$.each(lines, function(idx, my_line){
												var my_values = my_line.split("|");
												var to_add = [];
												to_add.push(my_values[0]);
												to_add.push(my_values[1]);
												to_add.push(my_values[2]);
												to_add.push(parseInt(my_values[3]));
												LabelSwapArray.push(to_add);
											});
											LS_array_counter = lines.length;
											break;
											case "LS_c_p_Add":
											//The following builds the LS_counters_per_Add array:
											if(!AllowLS) break;
											LS_counters_per_Add = [];
											var my_val = data.queryres.opt_value[i];
											var lines = my_val.split("||");
											
											// Transform one line encoded array to an array using the rules detailed above
											$.each(lines, function(idx, my_line){
												var my_values = my_line.split("|");
												//my_values loads 5 values: the first one is the LS assignment's description, the second and third are the first and last indices of rows in LabelSwapArray that correspond to the specific assignment. The 4th and 5th are the two labels that are swapped
												var to_add = [];
												to_add[0] = my_values[0];
												to_add[1] = [];
												for(var i = parseInt(my_values[1]); i<= my_values[2]; i++)
												{
													to_add[1].push(i);
												}
												to_add[2] = [];
												to_add[2].push(my_values[3]);
												to_add[2].push(my_values[4]);
												LS_counters_per_Add.push(to_add);
												
												//Add the respective lines in LSLabelSwaps
												$("#LSLabelSwaps").append("<option value='" + my_values[0] + "'>" + my_values[0] + "</option>");
											});
											break;
											case "Select_Labels":
											//Selects only the labels that are written in the database
											var my_val = data.queryres.opt_value[i];
											my_lbls_toselect = my_val.split("|");
											break;
										}
									}
								}
								
								// Refresh the reps the fractions and the conditions
								
								set_reps();
								refresh_fractions();
								if (isLabelFree) refresh_LFQ_conditions_from_test_data();
								$("#s3expparams").animate({scrollTop: 0}, "slow");
						});
				}
				else
				{
					$("#s2btnf").prop('disabled', true);
				}
			}
			}).fail(function (jqXHR, textStatus, errorThrown) {
			msgbox("An AJAX error occurred: " + errorThrown);
		});
	}
	
	
	function ProteoSignInit() {
		
		// Web Page initialization
		
		// If label swap is disabled hide the corresponding button from the UI
		if(!AllowLS)  $("#s3showdivLabelSwapRow").hide();
		
		// Use a regex to gather all forward (next) and backward buttons in an array
		var forward_buttons = getItems("button.main", /s[0-9]+btnf/);
		var backward_buttons = getItems("button.main", /s[0-9]+btnb/);
		
		// Bind the click event to "toggleNextClass" for each "next button"
		forward_buttons.forEach(function (btn) {
			$(btn).on("click", function () {
				
				// idsToStrings(".main_div .main_section") returns an array of all items having main_div or main_section as classnames
				// this matches all divisions in html corresponding to different steps of the process (s1div, s2div etc.)
				// in order of appeareance in HTML source
				// toggleNextClass hides the currently displayed division and shows the next one in order
				toggleNextClass(idsToStrings(".main_div .main_section"), "hidden", true, executeStage);
				
				if ($(backward_buttons[0]).hasClass("hidden")) {
					$(backward_buttons[0]).removeClass("hidden");
				}
				// Always redraw the data tables after displaying or hiding a step
				rawfiles_tbl_allfiles_DT.columns.adjust().draw();
			});
		});
		
		
		// Binds the click event to "toggleNextClass" for each "backward button" - for more info see forward_buttons.forEach above
		// In two cases pressing the back key (Step 2 and step 4) fires the "Careful Back" function to
		// prompt the user that they may lose changes if proceding
		backward_buttons.forEach(function (btn) {
			if ($(btn).attr("id") == "s22btnb")
			{
				$(btn).on("click", function () {
					carefulBackStep = 2;
					CarefulBack();
				});
			}
			else if($(btn).attr("id") == "s4btnb")
			{
				$(btn).on("click", function () {
					carefulBackStep = 4;
					CarefulBack();
				});
			}
			else
			{
				$(btn).on("click", function () {
					toggleNextClass(idsToStrings(".main_div .main_section").reverse(), "hidden", true, rollbackStage);
					rawfiles_tbl_allfiles_DT.columns.adjust().draw();
				});
			}
		});
		
		
		// Bind event for file upload button (this is necessary because upload buttons can not be stylized thus a hidden upload button is bound to the shown stylized one)
		$("#s2btnupld").on("click", function () {
			$("#__s2btnupld").click();
		});
		
		
		// Bind event when file(s) is/are chosen in step 2
		$("#__s2btnupld").change(function()
			{
				uploadButtonStatusChanged(this.files);
			});
			
			// Bind event to upload parameters button
			$("#__upldparams").change(uploadParameters);
			
			// Bind a function to s3showhideadvparams to toggle advance parameter visibility
			$("#s3showhideadvparams").on("click", showHideAdvancedParameters);
			
			// Bind a function to s3showdivLabelSwap to Show the Label; swap dialog box when the relative text is pressed
			$("#s3showdivLabelSwap").on("click", showLabelSwapDialogBox); 
			
			$("#s3showdivAdvancedOptions").on("click", onShowAdvOptionsClick);
			
			$("#s2linkSaveParams").on("click", function () {
				SaveParams(0);
			});
			
			// Bind functions to Save parameter and Load parameter buttons for all steps
			$("#s3linkSaveParams").on("click", function () {
				SaveParams(0);
			}); 
			$("#s4linkSaveParams").on("click", function () {
				SaveParams(0);
			}); 
			$("#s2linkLoadParams").on("click", function () {
				getparamfile();
				my_step = 2;
			}); 
			$("#s3linkLoadParams").on("click", function () {
				getparamfile();
				my_step = 3;
			}); 
			
			
			// Toggle visibility of quantitation filtering options based on whether the Quantitation Filtering check box
			// is checked or not:
			$("#s3advparams input[name='expquantfilt']").on("click", onclickQuantitationFilteringChkBox);
			
			// Reset CSS when losing focus (because some required fields might have been highlighted)
			$("#s3expparams input[data-required]").on("focusout", resetCssStyles(this));
			
			// Display chopped long raw file names on datatable on a tooltip
			$('#rawfiles_tbl_allfiles').delegate('td', 'mouseenter', longRawFileNameDisplay);
			
			// btnAssignExpStructCoord handles the button in Stage 2 that assigns a biorep and a tech rep
			// to the selected raw file
			$('#btnAssignExpStructCoord').on("click", assignExpStructCoord);
			
			// The following button resets the tabe in stage 2 and the relative variables
			$('#btnResetExpStructCoord').on("click", resetTableStage2);
			
			// Bind mouseleave and clock events to star scoring system
			
			for (var i = 1; i <= 5; i++)
			{
				var myid = "#star" + i;
				$(myid).mouseleave(function(){
					if (!stars_enabled) return;
					star_hovered = star_chosen;
					display_star_score();
				});
				$(myid).click(function(){
					if (!stars_enabled) return;
					star_chosen = star_hovered;
					display_star_score();
				});
			}
			
			// Bind to mousenter events to display star scoring system
			$("#star1").mouseenter(function() {starHover(1);});  
			$("#star2").mouseenter(function() {starHover(2);});  
			$("#star3").mouseenter(function() {starHover(3);});  
			$("#star4").mouseenter(function() {starHover(4);});  
			$("#star5").mouseenter(function() {starHover(5);});  
			
			if (!stars_enabled)
			{
				$("#Feedbackstars").css("display", "none");
			}
			
			// In assign condition dialog (Stage 2) in label free data there are 2 dialog buttons (OK and Cancel)
			// here they are both bound to a simple function that closes the dialog box
			// but OK is also boud to ons2LFQConditionsOK that performs the procedure of assigning 
			// anew label to the selected raw files
			
			$("#s2LFQConditionsOK").on("click", dlgFadeout);
			
			$("#s2LFQConditionsCancel").on("click", dlgFadeout);
			
			// For user info dialog OK button:
			
			$("#s4UserInfoOK").on("click", dlgFadeout);
			
			// CSS
			$(".tooltip").hover(function () {
				$(".tooltip span").css({"margin-left": -$(".tooltip span").width() / 2 + 9});
				$(".callout").css({"left": $(".tooltip span").width() / 2});
			});
			
			$(".tooltip2").hover(function () {
				$(".tooltip2 span").css({"margin-left": 11});
				$(".callout2").css({"left": 0});
			});
			
			
			// Bind test datasets dialog buttons' events to functions
			$("#dlgTestDatasetsBtnOK").on("click", datasetDialogOK);
			
			$("#dlgTestDatasetsBtnCancel").on("click", dlgFadeout);
			
			$("#dlgTestDatasetsBtnDownload").on("click", datasetDialogDownload);
			
			// Bind Label Swap dialog buttons' events to functions
			
			$("#LSBack").on("click", dlgFadeout);
			
			// Initialize Raw file DataTable
			rawFileDataTableInit();
			
			// Initialize the rawfile data table context menu
			rawFileDataTableCMenuInit();
			
			// Initialize the data table in Replication Multiplexing dialog
			RMDataTableInit();
			
			// Initialize the context menu in RM
			RMCMenuInit();
			// Initialize the context menu of label merging
			conditionsMergeCMenuInit();
			
			// Some CSS styling:
			$('#rawfiles_tbl_allfiles').css({border: 'none'});
			$('#rawfiles_tbl_allfiles').css({margin: '0px'});
			$('#rawfiles_tbl_allfiles thead th').css({border: 'none'});
			$('#rawfiles_tbl_allfiles').css({'border-bottom': 'none'});
			$('#rawfiles_tbl_allfiles').css({'border-bottom': 'none'});
			
			// Avoid the need to hold Ctrl in order to select multiple rows in rawfiles data table:
			//from http://stackoverflow.com/questions/8641729/how-to-avoid-the-need-for-ctrl-click-in-a-multi-select-box-using-javascript
			$("#conditions_list").mousedown(function(e){
				e.preventDefault();
				if (e.which == 3)
				{
					return;
				}
				var select = this;
				var scroll = select .scrollTop;
				
				e.target.selected = !e.target.selected;
				
				setTimeout(function(){select.scrollTop = scroll;}, 0);
				
				$(select ).focus();
			}).mousemove(function(e){e.preventDefault()});
			
			// Bind function to user feedback text box value change
			$("#userFeedback").on('change keyup paste', onuserFeedbackchange);
			
			// Initialize the test datasets:
			
			postTestDatasetsInfo();
			
			postClientServerClientInfo();
	}
	
	
	function uploadButtonStatusChanged(myFiles) {
		
		// When a file is chosen for upload in Step 2 automatically start uploading it to the server
		
		
		if ($("#__s2btnupld").val() == "") {
			// If no files were chosen simply stop execution
			return;
		}
		// else...
		
		//myFiles corresponds to the files argument of __s2btnupld button that are the local files chosen by the user
		
		if (myFiles.length > 0)
		{
			
			//Start uploading ...
			uploadingFiles = myFiles; // uploadingFiles is a global variable containing all files to be uploaded
			
			if (window.File && window.FileReader && window.FileList && window.Blob) //check that the browser support all necessary APIs
			{
				oversized_files_idxs = []; // this array will keep the oversized files indexes
				$.each(uploadingFiles, function (idx, file_i)
					{
						if(file_i.size > 2147483648)
						{
							oversized_files_idxs.push(idx); 
						}
					});
			}
			else
			{
				msgbox("WARNING: File APIs not supported by your current browser.");
			}
			
			if(oversized_files_idxs.length == 0)
			{
				// Call upload files to upload the files one by one. the arguments passed are:
				// serverside = false (the files are uploaded by the user and not "pseudouploaded" straight from the server)
				// the second argument is a delegate function that will be executed at the end of a successful (asynchronous)
				// file upload
				uploadFiles(false, function () {
					
					// After succesfully uploading a file:
					// First check what type of files has the user uploaded:
					
					// MQ are maxquant evidence files
					// MQP maxquant proteinGroups
					// PD proteome discoverer psm files
					// The possible valid combinations are MQ - MQP or PD alone
					//
					// If data originate from PD check a hidden check box - exppddata - stating that the data originate from PD (the check box's status will be used
					// in stage 4 to determine if data originate from PD or not)
					// and alter quantitation_prog_lbl a label that shows the user if the data were  quantified
					// by MQ or PD
					
					var Files_valid = false;
					
					if (types_of_files_uploaded.indexOf("MQ") != -1 && types_of_files_uploaded.indexOf("MQP") != -1)
					{
						
						if(types_of_files_uploaded.indexOf("PD") == -1)
						{
							//Maxquant valid files:
							Files_valid = true;
							procProgram = "MQ";
							$("#s3expparams input[name='exppddata']").prop('checked', false);
							$("#quantsoftdetectedspan").text("MaxQuant");
							$("#quantitation_prog_lbl").text("\u2014 Raw data were quantified with MaxQuant");
						}
						else
						{
							procProgram = "";
							//Both MQ and PD files were uploaded:
							msgbox("WARNING: Some of the files you uploaded are processed by MaxQuant and some by Proteome Discoverer, we are sorry to inform you that it is impossible to proceed.");
							$("#s2btnf").prop('disabled', true);
							Files_valid = false;
						}
					}
					if(types_of_files_uploaded.indexOf("PD") != -1)
					{
						if (types_of_files_uploaded.indexOf("MQ") == -1 && types_of_files_uploaded.indexOf("MQP") == -1)
						{
							//PD valid files:
							Files_valid = true;
							procProgram = "PD";
							$("#s3expparams input[name='exppddata']").prop('checked', true);
							$("#quantsoftdetectedspan").text("Proteome Discoverer");
							$("#quantitation_prog_lbl").text("\u2014 Raw data were quantified with Proteome Discoverer\u2122");
						}
						// No need to check if both MQ and PD files were uploaded, this has already been done
					}
					//If the user has uploaded more than one file from the same type alert them:
					if ((procProgram == "MQ" && types_of_files_uploaded.filter(filterFunction).length > 2) 	|| (procProgram == "PD" && types_of_files_uploaded.filter(filterFunction).length > 1))
					{
						msgbox("WARNING: You uploaded too many files from " + procProgram + "! Please make sure you upload MultiConsensus files. Please refresh the current page and try again.");
						$("#s2btnf").prop('disabled', true);
						Files_valid = false;
					}
					
					// If no files are to be uploaded and files are of a valid combination enable next button in Step 2
					if (Files_valid == true && (nToUpload == 0 && (isLabelFree == true || peptideLabelsNamesFromFile.length > 0)))
					{
						$("#s2btnf").prop('disabled', false);
					}
					
				}); //END uploadFiles postSuccess function		
			}
			else //there are oversized files
			{
				msgbox("WARNING: At least one of the files chosen is oversized, we are sorry to inform you that it is impossible to proceed.");
				$("#s2btnf").prop('disabled', true);
			}
		}
		
		$("#__s2btnupld").val(""); // Reset __s2btnupld
	}
	
	
	var uploadFiles = function (serverSide, postSuccess) {
		
		// uploadFiles is used in Stage 2 to upload the files to analyse to the server
		// this can happen either through literally uploading from the client computer to the server
		// or by pseudo-uploading i.e. copying a test dataset's files from the test-data folder of the server
		// to the current work-folder.
		// In the first case serverSide == false otherwise it is true
		
		// postSuccess is a function that will be executed after the file is completely uploaded
		// it will be simply be passed as an argument to postFile from here
		
		$("#s2btnf").prop('disabled', true); // Set forward button to false until the files have been uploaded
		
		// Create a status bar for each of the files to be uploaded
		$.each(uploadingFiles, function (idx, file_i) {
			var temp_index = uploaded_files + nToUpload + idx;
			if (serverSide)
			{
				$("#s2uluploaders table").append("<tr><td>" + file_i.toString() + "</td><td><progress max='100' value='0' id=uploadfile" + temp_index + "><div class='progress-bar'><span style='width: 80%;'></span></div></progress></td></tr>");
			}
			else
			{
				$("#s2uluploaders table").append("<tr><td>" + file_i.name + "</td><td><progress max='100' value='0' id=uploadfile" + temp_index + "><div class='progress-bar'><span style='width: 80%;'></span></div></progress></td></tr>");
			}
			
			
		}); // END each(uploadingfiles)
		
		// fire AJAX calls - postFile is the function that will handle the AJAX cal of upload_files.php for each file seperately
		$.each(uploadingFiles, function (idx, file_i)
			{
				var temp_index = uploaded_files + nToUpload;
				postFile(temp_index, file_i, serverSide, postSuccess); // this starts an asynchronous upload procedure calling "upload_files.php" to the server
				nToUpload++;
			}); // END each(uploadingFiles)
	}
	
	
	var validateParameters = function (params) {
		
		// Verifies that all parameters are valid
		// Return false if at least one required parameter (existence of data-required attr) is not set (empty text field)
		// Mark invalid parameters with a red border color (currently tested on input fields)
		// Marked fields are reset through the 'focusout' event (see respective binding in function document.ready)
		// params is an array of all the DOM elements that correspond to parameters
		
		var nValid = 0;
		var nInvalid = 0;
		$.each(params, function (idx, param_i)
			{
				switch ($(param_i).attr('type')) {
					case "text":
					// If the element is a text box
					if ($(param_i).attr('data-required') != null && $(param_i).attr('data-required'))
					{
						// If it is a required field
						if ($(param_i).val() == "" && !$(param_i).hasClass("hidden"))
						{
							// If the required field is empty and not hidden
							nInvalid++;
							// Set its border color to Red
							setItemAttrColor(param_i, "border", "#E60000");
						}
					}
					break;
					default:
					nValid++;
					break;
				}
			}); // END for each params
			var mytmp = "#conditions_list option:selected"; // Get all selected conditions list-items
			// Count how many conditions are selected
			var tmp_i = 0;
			$(mytmp).map(function () {
				tmp_i++;
			});
			
			// If only one condition was selected set conditions list border to Red
			if (tmp_i <= 1)
			{
				setItemAttrColor("#conditions_list", "border", "#E60000");
			}
			
			// In case the user has selected quantitation filtering make sure that the
			// lael he chose for quant filtering is valid (it appears in conditions list)
			var found_cond = false;
			if ($("#expquantfilt").prop("checked") == true)
			{
				var my_text = $("#expquantfiltlbl option:selected").text(); // The label on which to perform quantitation filtering
				$("#conditions_list > option:selected").each(function()
					{
						if (this.text == my_text)
						{
							found_cond = true;
						}
					});// END for each selected option for conditions list
					if (found_cond == false)
					{
						// If the quant filtering label is not valid paint the border of expquantfiltlbl red
						setItemAttrColor("#expquantfiltlbl", "border", "#E60000");
					}
			}
			else
			{
				found_cond = true;
			}
			
			// Inform the user if a condition has invalid characters
			
			if(!isIsobaricLabel)
			{
				var error_to_display = "";
				var re = new RegExp('^(?![0-9_])[a-zA-Z0-9_]+$');
				$("#conditions_list option").each(function(idx, my_opt)
					{
						if (!re.test($(this).val()))
						{
							error_to_display = error_to_display + "Condition name " + $(this).val() + " contains invalid characters or starts with a number or an underscore, please rename the condition to continue<br>"
						}
					});
					if (error_to_display !== "")
					{
						msgbox(error_to_display);
						nInvalid++;
					}
			}
			
			// If all parameters are valid return true
			return (nInvalid == 0 && tmp_i > 1 && found_cond);
	}
	
	
//}
//----------------------------------------
// MAIN ROUTINES __ END



// HELPER ROUTINES __ START
//----------------------------------------
//{
	
	
	var add_LFQ_conditions = function(tab_sep_string) {
		
		// Used by load params to populate the RawFileConditions array
		//
		// RawFilesConditions is represented in one-line string as described in add_raw_file_structure's comments
		// add_LFQ_conditions folds the argument string to an array
		
		RawFileConditions = [];
		// First split the first dimension of the string to an array
		var lines = tab_sep_string.split("\t\t");
		$.each(lines, function(idx, my_line){
			// For each row of the array split the string even more to create the columns of the array
			if (my_line == "") return;
			var my_vals = my_line.split("\t");
			// Get rid of the header:
			if(my_vals[0] == "rawfile")
			{
				return;
			}
			RawFileConditions.push({name: my_vals[0], condition: my_vals[1]});
		}); // END for each line
		
		//Show everything back to the user
		var all_items = $('#rawfiles_tbl_allfiles').find('tr');
		for (var i = 0; i < all_items.length; i++) {
			var items_tds = $(all_items[i]).find('td');
			var items_cond = items_tds[4];
			$.each(RawFileConditions, function (idx, my_raw_file)
				{
					if (my_raw_file.name == $(items_tds[0]).text())
					{
						$(items_cond).text(my_raw_file.condition);
						return false;
					}
				});
		}
		refresh_fractions();
		set_s22btnf();
		refresh_LFQ_conditions_from_test_data();
	}
	
	
	var add_raw_file_structure = function(tab_sep_string, check_validity) {
		
		// This function prints all rawfile names if tab_rep_string == "" or sets the rawfile_structure otherwise
		// based on tab seperated data where each \t\t denotes a line break and each \t a column break
		// e.g OT2_Terhune_2012-10-09_DMC-HLNF01_01\t1\t1\t1\t1\t\tOT2_Terhune_2012-10-09_DMC-HLNF01_02\t1\t2\t1\t1\t1\t\tOT2_Terhune_2012-10-09_DMC-HLNF01_02\t1\t2\t1\t1
		// will be transformed to:
		//
		// +--------------------------------------+---+---+---+
		// | OT2_Terhune_2012-10-09_DMC-HLNF01_01 | 1 | 1 | 1 |
		// +--------------------------------------+---+---+---+
		// | OT2_Terhune_2012-10-09_DMC-HLNF01_02 | 1 | 2 | 1 |
		// +--------------------------------------+---+---+---+
		//
		// with both rawfiles used = true
		//
		// This function can also be used to check if at some point all raw files in rawfiles table
		// exist in the uploaded file. This is useful when loading parameters from a saved file
		// to avoid parameter file - experiment incompatibilities
		//
		
		if(typeof(tab_sep_string) == "undefined" || tab_sep_string == "")
		{
			// tab_sep_string is programmed to be left undefined or set to "" only for debugging purposes
			var all_items = $('#rawfiles_tbl_allfiles').find('tr');
			for (var i = 0; i < all_items.length; i++) {
				var items_tds = $(all_items[i]).find('td');
				console.log($(items_tds[0]).text() + "\n");
			}
			return;
		}
		
		// Get all raw files from the respectve table
		var all_items = $('#rawfiles_tbl_allfiles').find('tr');
		var local_rawfiles = [];
		for (var i = 0; i < all_items.length; i++) {
			var items_tds = $(all_items[i]).find('td');
			local_rawfiles.push($(items_tds[0]).text());
		}
		
		// Refresh rawfile_structure
		// if check_validity == true then we are just making sure that there is no rawfile in the param save file that does not correspond to any files in the dataset
		// also we accept that some rawfiles in te dataset do not have a match in the param save file but not all of them
		if (!check_validity) rawfiles_structure = [];
		
		var lines = tab_sep_string.split("\t\t");
		var validityresponse = true;
		$.each(lines, function(idx, my_line){
			if (my_line == "") return;
			var my_vals = my_line.split("\t");
			if(my_vals[0] == "rawfile")
			{
				return;
			}
			
			if (check_validity == true)
			{
				//in this case the only reason to run the function is to check if all the rawfiles in the file uploaded exist in the current dataset
				if($.inArray(my_vals[0], local_rawfiles) == -1)
				{
					//if one rawfile in the params file is not in the dataset abort
					validityresponse = false;
				}
			}
			else
			{
				// if check validity == false then we should renew the rawfiles structure object
				rawfiles_structure.push({rawfile: my_vals[0], biorep: my_vals[1], techrep: my_vals[2], fraction: my_vals[3], used: my_vals[4]});
			}
		}); // END for each line
		
		// Now, validityresponse is true if all rawfiles in the param save file have a match in the dataset
		if (check_validity)
		{
			return validityresponse;
		}
		
		// Code from now on will run only if check_validity == false
		
		// Here the validation has ended and we know that there is at least one rawfile in the
		// parameter file that matches one in the dataset. Update these and return false if there are some rawfiles in the dataset left unset
		
		var rawfilesset = []; //dataset's raw files that will be set
		var all_items = $('#rawfiles_tbl_allfiles').find('tr'); // get all table items
		for (var i = 0; i < all_items.length; i++) {
			var items_tds = $(all_items[i]).find('td');
			$.each(rawfiles_structure, function (idx, my_raw_file)
				{
					if (my_raw_file.rawfile == $(items_tds[0]).text())
					{
						rawfilesset.push(my_raw_file.rawfile);
						return false;
					}
				}); // END for each rawfile_structure
		}
		
		var rawfilesunset = []; //dataset's raw files that will not be set
		for (var i = 0; i < all_items.length; i++)
		{
			var items_tds = $(all_items[i]).find('td');
			if ($.inArray($(items_tds[0]).text(), rawfilesset) == -1)
			{
				if ($(items_tds[0]).text() != "") rawfilesunset.push($(items_tds[0]).text());
			}
		}
		var found_unset = 0;
		for (var i = 0; i < all_items.length; i++) {
			var items_tds = $(all_items[i]).find('td');
			var items_brep = items_tds[1];
			var items_trep = items_tds[2];
			var items_frac = items_tds[3];
			if ($(items_tds[0]).text() == "") continue;
			if ($.inArray($(items_tds[0]).text(), rawfilesunset) != -1)
			{//if found unset raw in table create an instance of it since it means that this raw was not set at the time of creating the param file. So set the raw back to its original value (-,-,-,used)
				$(items_brep).text("-");
				$(items_trep).text("-");
				$(items_frac).text("-");
				$(items_tds[0]).css("text-decoration", "none");
				$(items_tds[1]).css("text-decoration", "none");
				$(items_tds[2]).css("text-decoration", "none");
				$(items_tds[3]).css("text-decoration", "none");
				if (isLabelFree == true)
				{
					$(items_tds[4]).css("text-decoration", "none");
				}
				found_unset++;
			}
		}
		
		
		//Show everything back to the user
		var all_items = $('#rawfiles_tbl_allfiles').find('tr');
		// console.log(all_items);
		for (var i = 0; i < all_items.length; i++) {
			var items_tds = $(all_items[i]).find('td');
			var items_brep = items_tds[1];
			var items_trep = items_tds[2];
			var items_frac = items_tds[3];
			$.each(rawfiles_structure, function (idx, my_raw_file)
				{
					if (my_raw_file.rawfile == $(items_tds[0]).text())
					{
						$(items_brep).text(my_raw_file.biorep);
						$(items_trep).text(my_raw_file.techrep);
						$(items_frac).text(my_raw_file.fraction);
						if(my_raw_file.used == false)
						{
							$(items_tds[0]).css("text-decoration", "line-through");
							$(items_tds[1]).css("text-decoration", "line-through");
							$(items_tds[2]).css("text-decoration", "line-through");
							$(items_tds[3]).css("text-decoration", "line-through");
							if (isLabelFree == true)
							{
								$(items_tds[4]).css("text-decoration", "line-through");
							}
						}
						else
						{
							$(items_tds[0]).css("text-decoration", "none");
							$(items_tds[1]).css("text-decoration", "none");
							$(items_tds[2]).css("text-decoration", "none");
							$(items_tds[3]).css("text-decoration", "none");
							if (isLabelFree == true)
							{
								$(items_tds[4]).css("text-decoration", "none");
							}
						}
						return false;
					}
				});
		}
		refresh_fractions();
		set_s22btnf();
	}
	
	
	var array_to_string = function(my_array, curdimension) {
		
		// array_to_string is a recursive function that gets an array of curdimension dimensions and
		// transforms it to a one line string
		//
		// Warning! in contrast to other PS one line representations that use as many repetitions of the
		// delimiter as the corresponding dimension, array_to_string (and its "counterpart" string_to_array)
		// 
		// e.g. the table
		//
		// +--------------------------------------+---+---+---+
		// | OT2_Terhune_2012-10-09_DMC-HLNF01_01 | 1 | 1 | 1 |
		// +--------------------------------------+---+---+---+
		// | OT2_Terhune_2012-10-09_DMC-HLNF01_02 | 1 | 2 | 1 |
		// +--------------------------------------+---+---+---+
		//
		// is transformed to
		//
		// OT2_Terhune_2012-10-09_DMC-HLNF01_01\t1\t1\t1\t1\t\tOT2_Terhune_2012-10-09_DMC-HLNF01_02\t1\t2\t1\t1\t1\t\tOT2_Terhune_2012-10-09_DMC-HLNF01_02\t1\t2\t1\t1
		//
		// in other functions such as add_raw_file_structure
		// but in array_to_string it will be transformed to
		//
		// OT2_Terhune_2012-10-09_DMC-HLNF01_01||1||1||1|OT2_Terhune_2012-10-09_DMC-HLNF01_02||1||2||1
		//
		// This difference was maintained due to backward compatibility issues
		//
		// array_to_string also suppports jagged arrays and arrays of up to 15 dimensions
		//
		// Notice that the input array should not contain the | character in any of its elements
		//
		// Also notice that array_to_string should be called from other functions without defining the curdimension var
		
		var retstring = "";
		var curdelimiter = "|";
		if (typeof(curdimension) === 'undefined') curdimension = 1;
		if (curdimension > 15) return;
		
		// First create a repetition of the delimeter representing the highest dimension seperator (e.g. in 2-dim array "||")
		curdelimiter = curdelimiter.repeat(curdimension);
		
		// Unfold the first dimension of the array:
		for (var i = 0 ; i <= my_array.length - 1; i++)
		{
			// For each element of my_array
			if (my_array[i].constructor === Array)
			{
				// If the element is an array unfold the array using array_to_string recursively 
				retstring += array_to_string(my_array[i], curdimension + 1) + curdelimiter;
			}
			else
			{
				// Otherwise simply add the element value to retstring
				retstring += my_array[i] + curdelimiter;
			}
		}
		// Clear the end of retstring from garbage delimiters
		retstring = retstring.substring(0, retstring.length - curdelimiter.length);
		return retstring;
	}
	
	
	function ArrNoDupe(a) {
		//from http://stackoverflow.com/questions/6940103/how-do-i-make-an-array-with-unique-elements-i-e-remove-duplicates
		//ArrNoDupe takes an array and returns anotehr one conatining only the unique elements of the original
		var temp = {};
		for (var i = 0; i < a.length; i++)
		temp[a[i]] = true;
		var r = [];
		for (var k in temp)
		r.push(k);
		return r;
		
		
	}
	
	
	function assignExpStructCoord() {
		
		// A function bound to btnAssignExpStructCoord that assigns a tech rep and a bio rep
		// to the selected raw file on Stage 2 table
		// If the user has solely filled the biorep textbox and left the techrep blank all selected
		// rawfiles are assigned the same biorep and all rawfiles are assigned to a different techrep
		//
		
		// assignExpStructCoord first alters data in rawfiles_structure and refreshes Stage2 table
		// rawfiles_structure is an array of objects having the following structure:
		// {rawfile: ..., biorep: ..., techrep: ..., fraction: ..., used: ...}
		// as an example rawfiles_structure might seem like
		// 0: {rawfile: "D131106_014.raw", biorep: 1, techrep: 1, fraction: 1, used: true}
		// 1: {rawfile: "D131106_016.raw", biorep: 2, techrep: 1, fraction: 1, used: true}
		// 2: {rawfile: "D131106_017.raw", biorep: 3, techrep: 1, fraction: 1, used: true}
		// since rawfiles_structure keeps the rawfiles name as an argument (rawfile) it is difficult
		// to access a specific emelent of rawfile_structure based on its name, therefore code refering to
		// this array may contain many for each loops
		
		// First get all selected raw files lines in the table and store them to items:
		
		
		var items = $('#rawfiles_tbl_allfiles').find('.rawfiles_tbl_td_selected');
		
		// Get the biorep and techrep the user typed:
		var def_biorep = Number($('#expstructcoord_biorep').val());
		var def_techrep = Number($('#expstructcoord_techrep').val());
		var curr_biorep = def_biorep;
		var curr_techrep = def_techrep;
		var curr_fraction = 0;
		var rep_offset = 1;
		$("#btnResetExpStructCoord").prop('disabled', false);
		
		for (var i = 0; i < items.length; i++)
		{
			// for each line in table in Stage 2 selected
			
			var items_tds = $(items[i]).find('td'); // The table line's cells
			var items_biorep = items_tds[1]; // Cell containing the biorep of the selected line
			var items_techrep = items_tds[2]; // Cell containing the techrep of the selected line
			var items_frac = items_tds[3]; // Cell containing the fraction of the selected line
			if (def_biorep == 0)
			{
				return;
			}
			
			//if the user attempted to overwrite the coordinations of a file, delete its rawfiles_structure entry and set it again
			$.each(rawfiles_structure, function (idx, my_raw_file)
				{
					// Find the rawfiles_structure entry of teh selected raw file
					if(my_raw_file.rawfile == $(items_tds[0]).text())
					{
						// Erase the entry
						rawfiles_structure.splice(idx, 1);
						return false;
					}
				}); // END for each rawfile_structure
				
				// Set the rawfiles_structure values and perform auto completion
				// rep_counts is an array that is used to avoid duplicates in replicates autocompletion
				// see set_reps comments for more information
				
				if (def_biorep > 0)
				{ // i.e. non-blank
					if (def_techrep == 0)
					{ // blank
						//techreps auto completion
						if (def_biorep in rep_counts["biorep"] && !isLabelFree)
						{
							curr_techrep = (rep_counts["biorep"][def_biorep].techrep.length - 1) + rep_offset++;
						}
						else
						{
							curr_techrep = rep_offset++;
						}
					}
					else
					{
						// fractions auto completion
						if (def_biorep in rep_counts["biorep"] && def_techrep in rep_counts["biorep"][def_biorep].techrep && !isLabelFree)
						{
							curr_fraction = (rep_counts["biorep"][def_biorep].techrep[def_techrep].fraction.length - 1) + rep_offset++;
						}
						else
						{
							curr_fraction = rep_offset++;
						}
					}
				}
				else
				{// blank defined biorep
					if (def_techrep == 0)
					{ // blank defined tech rep
						// bioreps, from current count to items.length, except if we have label-free data
						if (!isLabelFree)
						{
							curr_biorep = rep_counts["biorep"].length + rep_offset++;
						}
						else
						{
							curr_biorep = rep_offset++;
						}
					}
					else
					{
						// error, a biorep must be specified
						if (rawfiles_structure.length == 0)
						{
							$("#btnResetExpStructCoord").prop('disabled', true);
						}
						return;
					}
				}
				curr_used = true;
				
				// Push the rawfiles_structure entry
				
				rawfiles_structure.push({rawfile: $(items_tds[0]).text(), biorep: curr_biorep, techrep: curr_techrep, fraction: curr_fraction, used : curr_used});
				
				// Since we just assignes breps and treps to the rawfile it is obvious that it will be used in the analysis so
				// erase all text decoration (strikeout) from the corresponding lines in table in Stage 2
				
				$(items_tds[0]).css("text-decoration", "none");
				$(items_tds[1]).css("text-decoration", "none");
				$(items_tds[2]).css("text-decoration", "none");
				$(items_tds[3]).css("text-decoration", "none");
				if (isLabelFree == true) $(items_tds[4]).css("text-decoration", "none");
				
				// Note that the table displays - in case a brep trep or frac is not set when rawfile_structure has a zero value
				$(items_biorep).text(curr_biorep == 0 ? '-' : curr_biorep);
				$(items_techrep).text(curr_techrep == 0 ? '-' : curr_techrep);
				$(items_frac).text(curr_fraction == 0 ? '-' : curr_fraction);
				
		}
		
		set_reps();
		refresh_fractions();
		$('#rawfiles_tbl_allfiles tbody tr').removeClass('rawfiles_tbl_td_selected');
	}
	
	
	var CarefulBack = function() {
		// CarefulBack warns the user that they might lose changes if they move one stepo back in PS
		// In case of step 2 and step 4 - after analysis completion, the user must be warned
		// that if they step back they might lose their progress (e.g. they should upload new files or rerun the analysis)
		// CarefulBack is called when the back button in the aforementioned steps is clicked
		
		if (carefulBackStep == 2)
		{
			// Show the yes/no dialog box
			$("#CarefulBackTitle").empty();
			$("#CarefulBackTitle").append("<p style='font-size: 85%'>Going back to Step 1 requires that you upload files again.<br>Are you sure you want to proceed?</p>");
			$(".expparamsDlg").css({"left" : ($("body").width()/2) - ($("#WarnCarefulBack").width()/2)});
			$('body').append('<div id="mask"></div>');
			$("#WarnCarefulBack").fadeIn(300);
			$('#mask').fadeIn(300);
		}
		else if (carefulBackStep == 4)
		{
			if (analysis_finished)
			{
				// Show the yes/no dialog box
				$("#CarefulBackTitle").empty();
				$("#CarefulBackTitle").append("<p style='font-size: 85%'>Going back to Step 3 will start a new session.<br>Please make sure you downloaded your results.<br>Are you sure you want to proceed?</p>");
				$(".expparamsDlg").css({"left" : ($("body").width()/2) - ($("#WarnCarefulBack").width()/2)});
				$('body').append('<div id="mask"></div>');
				$("#WarnCarefulBack").fadeIn(300);
				$('#mask').fadeIn(300);
			}
			else
			{
				toggleNextClass(idsToStrings(".main_div .main_section").reverse(), "hidden", true, rollbackStage);
				rawfiles_tbl_allfiles_DT.columns.adjust().draw();
			}
		}
	}
	
	
	var CheckParamsValidity = function(myparamsstring) {
		
		// CheckParamsValidity gets the multiline string myparamsstring that represents a parameter file
		// and validates all parameters
		//
		// If a parameter is invalid it sets the error_message string with an html formatted error message
		// CheckParamsValidity returns the error message so if it is blank no errors were detected
		var error_message = "";
		
		// Break the string into lines
		var paramslines = myparamsstring.split("\n");
		var setvar = "";
		
		$.each(paramslines, function(idx, myparamline){
			
			// Compensate with Windows like line breaks
			if(myparamline[myparamline.length - 1] = "\r")
			{
				myparamline = myparamline.substring(0, myparamline.length - 1);
			}
			// Get rid of comments
			if (myparamline.charAt(0) == "#")
			{
				return;
			}
			
			// setvar contains the name of the variable to be set
			// if setvar is blank then the next line should be a variable name
			
			if (setvar == "")
			{
				// In this case we expect to read a var name
				if (myparamline == "isLabelFree" || myparamline == "isIsobaricLabel" || myparamline == "procprogram" || myparamline == "rawfiles_structure" || myparamline == "expid" || myparamline == "exptpoint" || myparamline == "conditions_to_compare" || myparamline == "quantitation_filtering" || myparamline == "filtering_label" || myparamline == "peptide_level_filtering" || myparamline == "LeastBRep" || myparamline == "LeastPep" || myparamline == "Pthreshold" || myparamline == "LFQconditions" || myparamline == "!Rename" || myparamline == "!LS_Array" || myparamline == "!LS_c_p_Add")
				{
					setvar = myparamline;
				}
			}
			else
			{
				// Here we have the var name in set var and we read the variable value
				
				// The loaded parameter file must correspond to the same type of experiment with
				// the currently loaded on PS (Label free and isobaric labelled:)
				if (setvar == "isLabelFree")
				{
					if (myparamline == "true")
					{
						if (isLabelFree == false)
						{
							error_message += "<li style='text-align: left;'>The parameters correspond to a labelled experiment set but the uploaded data set to a label-free<br></li>";
						}
					}
					else if (myparamline == "false")
					{
						if (isLabelFree == true)
						{
							error_message += "<li style='text-align: left;'>The parameters correspond to a label-free experiment set but the uploaded data set to labelled<br></li>";
						}
					}
				}
				else if (setvar == "isIsobaricLabel")
				{
					if (myparamline == "true")
					{
						if (isIsobaricLabel == false)
						{
							error_message += "<li style='text-align: left;'>The parameters correspond to an isobaric labelled experiment set but the uploaded data set does not<br></li>";
						}
					}
					else if (myparamline == "false")
					{
						if (isIsobaricLabel == true)
						{
							error_message += "<li style='text-align: left;'>The parameters correspond to an experiment which did not employ isobaric labeling but the uploaded data set corresponds to such an experiment<br></li>";
						}
					}
				}
				// The processing program should be the same with the one of the loaded experiment
				else if (setvar == "procprogram")
				{
					var myexppddata;
					var provenprocprogram = "";
					if ($("#s3expparams input[name='exppddata']").prop("checked") == false)
					{
						provenprocprogram = "MQ";
					}
					else
					{
						provenprocprogram = "PD";
					}
					if (myparamline == "MQ")
					{
						if (provenprocprogram == "PD")
						{
							error_message += "<li style='text-align: left;'>The parameters correspond to a data set processed by MaxQuant but the uploaded data set was processed by Proteome Discoverer<br></li>";
						}
					}
					else if (myparamline == "PD")
					{
						if (provenprocprogram == "MQ")
						{
							error_message += "<li style='text-align: left;'>The parameters correspond to a data set processed by Proteome Discoverer but the uploaded data set was processed by MaxQuant<br></li>";
						}
					}
				}
				// Make sure that all rawfiles in the dataset are included in the parameter's file
				else if (setvar == "rawfiles_structure")
				{
					// Here add_raw_file_structure is used to make sure that all rawfiles in the dataset are included in the parameter's file
					if (add_raw_file_structure(myparamline, true) == false)
					{
						error_message += "<li style='text-align: left;'>The parameters correspond to a data set with different raw files than the uploaded data set<br></li>";
					}
				}
				setvar = "";
			}
			
		}); // END for each parameter line
		
		return error_message;
		
	}
	
	
	function conditionsMergeCMenuInit() {
		
		// Initializes the context menu bound to conditions list right click
		// that allows the user to merge conditions
		
		if(AllowMergeLabels)
		{
			
			label_context_menu = $('#conds_list_container').contextmenu({
				// Define the items
				delegate: ".conditions_list",
				menu: [
					{title: "Same Condition", cmd: "same_cond"},
					{title: "Restore Conditions", cmd: "restore_conds"}
				],
				select: function(event, ui) {
					switch(ui.cmd){
						case "same_cond":
						// In the event of same condition click
						itemsToRename = [];
						$("#conditions_list option:selected").each(function(){
							// Get all selected conditions to rename them (merge them)
							itemsToRename.push($(this).val());
						});
						
						// Show the merging dialog
						$(".expparamsDlg").css({"left" : ($("body").width()/2) - ($("#condsadvancedDlg").width()/2)});
						$('body').append('<div id="mask"></div>');
						$("#condsadvancedDlg").fadeIn(300);
						$('#mask').fadeIn(300);
						
						$("#s3AdvancedOK").on("click", function () {
							// When the OK button is pressed on the Dialog Match the old conditions to their 
							// new renamed oned in RenameArray (see Refresh_conds_list for more information on renameArray'w structure)
							dlgFadeout();
							if($("#s3AdvNewCondition").val() == "") return;
							$.each(RenameArray, function(idx, my_row)
								{
									if($.inArray(my_row[0], itemsToRename) != -1)
									{
										RenameArray[idx][1] = $("#s3AdvNewCondition").val();
									}
								});
								$("#s3AdvNewCondition").val("");
								RenameFromTestData = false;
								// Refresh the items in conditions list
								Refresh_conds_list();
								//ADD OK RESULT HERE
						});
						
						$("#s3AdvancedCancel").on("click", function () {
							dlgFadeout();
							//ADD cancel RESULT HERE
						});
						break;
						case "restore_conds":
						// Restore the conditions of a merged condition
						itemsToRename = [];
						$("#conditions_list option:selected").each(function(){ itemsToRename.push($(this).val());});
						$.each(RenameArray, function(idx, my_row)
							{
								// Simply if we find the selected condition in the second column of RenameArray
								// set it back to its original value (the corresponding one from the first column)
								// See Refresh_conds_list for more information on RenameArray's structure
								if($.inArray(my_row[1], itemsToRename) != -1)
								{
									RenameArray[idx][1] = RenameArray[idx][0];
								}
							});
							RenameFromTestData = false;
							// Refresh the items in conditions list
							Refresh_conds_list();
							break;
					}
				},
				beforeOpen: function(event, ui) {
					var $menu = ui.menu,
					$target = ui.target,
					extraData = ui.extraData;
					ui.menu.zIndex(9999); // Bring to top
				}
			});
		}
	}
	
	
	var CreateParamsFile = function() {
		// Returns a multiline string with all current parameters encoded
		//
		// Each line starting with # is a comment
		// Each line that does not start with # contains the name of a variable
		// The immediately following line contains the value of this line
		//
		// Arrays are encoded using a tab seperated format in one line (see an example in add_raw_file_structure for more information)
		// Sometimes the seprator is not "\t" but e.g. "|"
		//
		// Some special commands may be set in the parameters files that begin with "!"
		// !Rename
		
		var mytext = "";
		mytext += "#ProteoSign parameters files for " + $('input[name="expid"]').val() + " (session ID: " + sessionid + ")\n#\n#File Information:\n";
		// Boolean variables:
		mytext += "isLabelFree\n";
		mytext += isLabelFree + "\n";
		mytext += "isIsobaricLabel\n";
		mytext += isIsobaricLabel + "\n";
		// Proc program (MW or PD)
		mytext += "procprogram\n";
		var quantlabel = $("#quantitation_prog_lbl").text();
		var pattern = new RegExp("MaxQuant");
		var res = pattern.test(quantlabel);
		var provenprocprogram = "";
		if (res == true)
		{
			provenprocprogram = "MQ";
		}
		else
		{
			provenprocprogram = "PD";
		}
		mytext += provenprocprogram + "\n";
		
		mytext += "#User Parameters and commands:\n";
		
		// Raw files structure:
		mytext += "rawfiles_structure\n";
		var myval = "";
		var mytemp = "";
		// If RM is used rawfilesstructure does not keep the real values of rawfiles
		// Create a rawfilesstructure - like array from the raw files table in this case
		if (RMisused)
		{
			var temp_rf_structure = [];
			var all_items = $('#rawfiles_tbl_allfiles').find('tr');
			for (var i = 1; i < all_items.length; i++)
			{
				var items_tds = $(all_items[i]).find('td');
				temp_rf_structure.push({rawfile: $(items_tds[0]).text(), biorep: "-", techrep: "-", fraction: "-", used: true});
			}
		}
		
		// rawfiles structure is a two dimension array, each raw is separated from the next one with \t\t
		// and each element from the next one with a simple \t, see add_raw_file_structure for
		// more info on this kind of notation
		
		if (!RMisused)
		{
			$.each(rawfiles_structure, function(idx, my_rawfile)
				{
					if (my_rawfile.used == true)
					{
						mytemp = "1";
					}
					else
					{
						mytemp = "0";
					}
					myval += my_rawfile.rawfile + "\t" + my_rawfile.biorep + "\t" + my_rawfile.techrep + "\t" + my_rawfile.fraction + "\t" + mytemp + "\t\t";
				}); // END for each rawfiles structure
		}
		else
		{
			// If RM is used use the rawfilesstructure - like array we built above
			$.each(temp_rf_structure, function(idx, my_rawfile)
				{
					if (my_rawfile.used == true)
					{
						mytemp = "1";
					}
					else
					{
						mytemp = "0";
					}
					myval += my_rawfile.rawfile + "\t" + my_rawfile.biorep + "\t" + my_rawfile.techrep + "\t" + my_rawfile.fraction + "\t" + mytemp + "\t\t";
				}); // END for each temp rf structure
		}
		mytext += myval + "\n";
		
		// Experiment ID information
		mytext += "expid\n";
		mytext += $('input[name="expid"]').val() + "\n";
		
		// Rename command tells Proteosign to rename the conditions when it loads the parameter file
		mytext += "!Rename\n";
		// !Rename's value is simply a one-line representation of RenameArray with | as seperator
		myval = "";
		$.each(RenameArray, function(idx, my_line)
			{
				myval += my_line[0] + "|" + my_line[1] + "||";
			}); // END for each line in renamearray
			mytext += myval + "\n";
			
			// Experiment time point (experiment description)
			mytext += "exptpoint\n";
			mytext += $('input[name="exptpoint"]').val() + "\n";
			
			// Conditions that will be compared
			mytext += "conditions_to_compare\n";//these are the authentic items in the conditions list
			myval = "";
			$.each(AuthenticItemsinRename, function(idx, my_cond)
				{
					myval += my_cond + "\t";
				}); // END for each authenticitemsinrename
				
				// Some Boolean values
				mytext += myval + "\n";
				mytext += "quantitation_filtering\n";
				if($("#expquantfilt").prop("checked") == true)
				{
					myval = "T";
				}
				else
				{
					myval = "F";
				}
				mytext += myval + "\n";
				mytext += "filtering_label\n";
				mytext += $("#expquantfiltlbl option:selected").val() + "\n";
				mytext += "peptide_level_filtering\n";
				if($("#expquantfiltprot").prop("checked") == true)
				{
					myval = "T";
				}
				else
				{
					myval = "F";
				}
				mytext += myval + "\n";
				
				// Advanced options
				mytext += "LeastBRep\n";
				mytext += LeastBRep + "\n";
				mytext += "LeastPep\n";
				mytext += LeastPep + "\n";
				mytext += "Pthreshold\n";
				mytext += Pthreshold + "\n";
				
				// Label free conditions simply saves RawFileConditions as a one-line representation (\t seperator)
				mytext += "LFQconditions\n";
				myval = "";
				$.each(RawFileConditions, function(idx, my_cond)
					{
						myval += my_cond.name + "\t" + my_cond.condition + "\t\t";
					});
					mytext += myval + "\n";
					
					// !LS_array stores the labelswaparray in one line representation with | as seperator
					mytext += "!LS_Array\n";
					myval = "";
					$.each(LabelSwapArray, function(idx, my_line)
						{
							myval += my_line[0] + "|" + my_line[1] + "|" + my_line[2] + "|" + my_line[3] + "||";
						}); // END for each label swap array element
						mytext += myval + "\n";
						
						// LS_c_p_Add stores the LS_counters_per_Add in one line representation with | as seperator
						// LS_c_p_Add is the second array needed for label swapping, note that it is a three dimension
						// variable, see its description in "Global variables" section for more details
						mytext += "!LS_c_p_Add\n";
						myval = "";
						$.each(LS_counters_per_Add, function(idx, my_line)
							{
								myval += my_line[0] + "||"; //add the description of the ls assignment
								$.each(my_line[1], function(idx, my_index) // the array of indices in this assignment
									{
										myval += my_index + "|";
									});
									myval += "|";
									myval += my_line[2][0] + "||";
									myval += my_line[2][1] + "|||";
							}); // END for each LS_counters_per_Add
							mytext += myval + "\n";
							
							// !Select_Labels store all selected conditions in conditions list as a one line representation of an array
							// with | as a seperator
							
							mytext += "!Select_Labels\n";
							myval = "";
							$("#conditions_list option").each(function(idx, my_opt)
								{
									if ($(my_opt).prop("selected") == true)
									{
										myval += $(my_opt).val() + "|";
									}
								}); // END for each selected condition in conditions list
								mytext += myval + "\n";
								
								// RM related variables:
								// exactly as most of the arrays are reresented in one line, RMrawfilesdata and RMtagsdata
								// are encoded with | as seperator using array_to_string
								
								mytext += "!RMrawfilesdata\n";
								mytext += array_to_string(RMrawfilesdata) + "\n";
								mytext += "!RMtagsdata\n";
								mytext += array_to_string(RMtagsdata) + "\n";
								
								// More RM related variables
								mytext += "!RMbrepsRepInRawFiles\n";
								mytext += RMbrepsRepInRawFiles + "\n";
								mytext += "!RMtrepsRepInRawFiles\n";
								mytext += RMtrepsRepInRawFiles + "\n";
								mytext += "!RMconditionsRepInRawFiles\n";
								mytext += RMconditionsRepInRawFiles + "\n";
								mytext += "!RMisused\n";
								mytext += RMisused + "\n";
								return mytext;
	}
	
	
	var create_my_all_mq_labels = function() {
		
		// Since MQ data analysis requires MQ labels to be sent to the backend in the order they appear in file's headers
		// a special var called my_all_mq_labels will be sent in such a case containing an R-structured array [c(... , ..., ...)]
		// to be sent to the backend
		// In case label merging - renaming (see our docimentation for more details) is enabled these labels are the first column
		// of RenameArray (see comments on InitializeRename) fore more information
		// Otherwise they are simply the labels displayed in the conditions list in Stage 3
		
		if (!RenameFromTestData)
		{
			
			// Get all elements in conditions list and populate my_all_mq_labels with them
			
			var myOpts = document.getElementById('conditions_list').options;
			my_all_mq_labels = "c(";
			for(i = 0; i < myOpts.length; i++)
			{
				if(i != myOpts.length - 1)
				{
					my_all_mq_labels = my_all_mq_labels + '"' + myOpts[i].value + '", ';
				}
				else
				{
					my_all_mq_labels = my_all_mq_labels + '"' + myOpts[i].value + '"';
				}
			}
			
			my_all_mq_labels = my_all_mq_labels + ")";
		}
		else
		{
			// Get the firs column in RenameArray and populate my_all_mq_labels with them
			my_all_mq_labels = "c(";
			for(i = 0; i < RenameArray.length; i++)
			{
				if(i != RenameArray.length - 1)
				{
					my_all_mq_labels = my_all_mq_labels + '"' + RenameArray[i][0] + '", ';
				}
				else
				{
					my_all_mq_labels = my_all_mq_labels + '"' + RenameArray[i][0] + '"';
				}
			}
			
			my_all_mq_labels = my_all_mq_labels + ")";
		}
	}
	
	
	function datasetDialogDownload() {
		
		// Bound on dataset dialog Download button (Stage 1), if a test dataset is selected it downloads
		// the corresponding dataset
		
		dlgFadeout();
		var tdname = $("#s1TestDatasetsSelection option:selected").text();
		if(tdname){
			downloadTestDataset(tdname);
		}
	}
	
	
	function datasetDialogOK() {
		
		// Bound to OK button on dialog of test datasets (Stage 1) datasetDialogOK hides the
		// mask under the dialog and loads the test dataset using postTestDatasetInfo
		dlgFadeout();
		postTestDatasetInfo($("#s1TestDatasetsSelection option:selected").text());
		$("#s2btnb").prop('disabled', true);
		$("#dlgTestDatasetsBtnOK").prop('disabled', true);
		datatestOK_clicked = true;
	}
	
	
	var display_star_score = function() {
		
		// Used to display highlighted all stars up to the one that the user hovers	
		if (!stars_enabled) return;
		for (var i = 1; i <= 5; i++)
		{
			var my_id = "#star" + i;
			if (i <= star_hovered)
			{
				$(my_id).css("opacity", "1");
			}
			else
			{
				$(my_id).css("opacity", "0.6");
			}
		}
	}
	
	
	var dlgFadeout = function () {
		
		// Executed when a dialog box is closed
		$(".expparamsDlg").fadeOut(300, function () {
			$('#mask').remove();
		});
	}
	
	
	var dlgFadeoutFeedback = function () {
		
		// Close the Feedback dialog box
		$(".expparamsDlgFeedback").fadeOut(300, function () {
			$('#maskFeedback').remove();
		});
	}
	
	
	var dlgFadeoutInfo = function () {
		
		// Executed when an info dialog box is closed
		$(".expparamsDlgInfo").fadeOut(300, function () {
			$('#maskInfo').remove();
		});
	}
	
	
	var dlgFadeoutRepMult = function () {
		
		// Close the RM dialog bpx
		$(".RepMultDlg").fadeOut(300, function () {
			$('#maskRepMult').remove();
		});
	}
	
	
	var downloadTestDataset = function (dataset_desc) {
		
		// Script to download a selected test dataset
		window.location = cgi_bin_path + 'download_test_data.php?session_id=' + sessionid + '&dataset_info_requested=' + dataset_desc;
	}
	
	
	function eraseInvalidLabelSwaps() {
		
		// Erases invalid label swaps in case the conditions of the expoeriment have been changed
		
		// First get all current valid labels:
		if (!AllowLS) return;
		var my_valid_options = [];
		$("#conditions_list option").each(function(idx, my_valid_opt)
			{
				my_valid_options.push($(my_valid_opt).val());
			}); // END for each conditions_list option
			
			var my_invalid_LSwaps = [];
			if (!RenameFromTestData)
			{
				// The last element of LS_counters_per_Add is an array of the labels that are swapped in a single swap assignment
				$.each(LS_counters_per_Add, function(idx, my_swap_assignment)
					{
						// In case at least one of the labels in an assignment is not valid remove it from the LSwaps and inform the user
						if($.inArray(my_swap_assignment[2][0], my_valid_options) == -1 || $.inArray(my_swap_assignment[2][1], my_valid_options) == -1)
						{
							my_invalid_LSwaps.push(my_swap_assignment[0]);
						}
					});
			}
			
			//Now my_invalid_LSwaps contains the descriptions of all the LSwaps assignments to be removed, select them in LSLabelSwaps and call onLSRemoveclick to remove them safely
			if(my_invalid_LSwaps.length >0 )
			{
				//Note that if an LS assignment has been removed manually or not in the past, it is candidate to be considered invalid again since no function erases any elements of LS_counters_per_Add.
				// So before trying to remove an LSassignment always check if it exists in LSLabelSwaps list
				var real_invalid_LSs = 0;
				var text_to_display = "";
				$("#LSLabelSwaps option").each(function(idx, my_opt)
					{
						if ($.inArray($(my_opt).val(), my_invalid_LSwaps) != -1)
						{
							$(my_opt).prop("selected", true);
							text_to_display = text_to_display + '<i>"' + $(my_opt).val() + '"</i>, ';
							real_invalid_LSs++;
						}
						else
						{
							$(my_opt).prop("selected", false);
						}
					}); // END for each LSLabelSwaps option
					
					if (real_invalid_LSs > 1)
					{
						onLSRemoveclick();
						text_to_display = "Label swaps: " + text_to_display;
						text_to_display = text_to_display + "are no longer valid and were removed";
						msgbox(text_to_display);//Prompt the user
					}
					else if(real_invalid_LSs == 1)
					{
						onLSRemoveclick();
						text_to_display = "Label swap: " + text_to_display;
						text_to_display = text_to_display + "is no longer valid and was removed";
						msgbox(text_to_display);//Prompt the user
					}
			}
	}
	
	
	var executeStage = function (stageIndex) {
		// Executed when a "Next" type of button is clicked execute stage performs a validity check before
		// advancing to the next stage. This applys only in stage 3 where the analysis is about to begin
		// in Case 5 where the anallysis has finished PS reloads itself and is self-initialized
		
		// executeStage is a togglefun for toggleNextClass
		
		var ret = true;
		switch (stageIndex) {
			case 3:
			var parameterObjs = $("#s3expparams input,#s3expparams select");
			if (validateParameters(parameterObjs)) { //if all parameters are valid
				postParameters(parameterObjs);
				// Clear areas where results information appears.
				resetResultsInfoItms();
				} else {
				ret = false;
			}
			break;
			case 5:	
			resetState();
			window.location.reload();
			break;
			default:
		}
		
		return ret;
	}
	
	
	function filterFunction(value) {
		
		// for filtering filetypes
		if (procProgram == "MQ")
		{
			return (value == "MQ" || value == "MQP");
		}
		else if( procProgram == "PD")
		{
			return (value == "PD");
		}
	}
	
	
	var generate_tab_file = function(my_array) {
		
		// This function gets a 2 dimensional array and transforms it to a tabular file
		var ret = "";
		for(var i = 0; i <= my_array.length - 1; i++)
		{
			// For each row
			for (var j = 0; j <= my_array[i].length - 1; j++)
			{
				// And for each column
				// break the columns using a tab seperator
				ret += my_array[i][j] + "\t";
			}
			ret = ret.substr(0, ret.length - 1);
			// Break the rows using a newline character
			ret += "\n";
		}
		ret = ret.substr(0, ret.length - 1);
		return ret;
	}
	
	
	var gen_expdesign = function(struct) {
		
		// Returns rawfiles structure in tab seperated format
		// Used to send the structure to R but can also be used for debugging purposes
		
		var ret = "";
		for (var i = 0; i < struct.length; i++) 
		{
			if(struct[i].used == true)
			{
				ret = ret + struct[i].rawfile + "\t" + struct[i].biorep + "\t" + struct[i].techrep + "\t" + struct[i].fraction + "\n"; 
			}
		}
		return ret;
	}
	
	
	var gen_lfqdesign = function (struct) {
		
		// Transforms RawFileConditions (an array that assigns every raw file to a condition)
		// to tab seperated format	
		// Used to send the structure to R but can also be used for debugging purposes
		
		var ret = "";
		for (var i = 0; i < struct.length; i++)
		{
			ret = ret + struct[i].name + "\t" + struct[i].condition + "\n"; 
		}
		return ret;
	}
	
	
	var gen_LSArray = function (struct) {
		
		// Transforms Label Swap array to a tab seperated format
		// Used to send data to R but can also used for debugging purposes
		
		if(!AllowLS) return;
		var ret = "";
		for (var i = 0; i < struct.length; i++) {
			ret = ret + struct[i][0] + "\t" + struct[i][1] + "\t" + struct[i][2] + "\n"; 
		}
		if (ret == "")
		{
			ret = "This run contains\tno label\tswaps\n"
		}
		return ret;
	}
	
	
	var gen_RenameFile = function (struct) {
		
		// Transforms Rename array to tab seperated format
		// Used to send data to R but can also be used for debugging purposes
		
		var ret = "";
		for (var i = 0; i < struct.length; i++) {
			ret = ret + struct[i][0] + "\t" + struct[i][1] + "\n"; 
		}
		return ret;
	}
	
	
	var getItems = function (className, id_pattern) {
		
		// Get array of items from the DOM using a selector and id pattern
		return $(className).toArray().filter(function (element) {
			var match;
			if($(element).attr("id")){
				match = ($(element).attr("id")).match(id_pattern);
			}
			else
			{
				match = null;
			}
			return (match == null ? false : true);
		});
	}
	
	
	var getparamfile = function() {
		
		// Perform a click to the hidden upload parameters button to load a parameters file
		$("#__upldparams").click();
	}
	
	
	var getRSS = function (rssurl, renderelem) {
		
		// Get RSS elements from RSS feed
		
		var thedata = new FormData();
		thedata.append('session_id', sessionid);
		thedata.append('rssurl', rssurl);
		$.ajax({
			url: cgi_bin_path + 'get_rss.php',
			type: 'POST',
			// Form data
			data: thedata,
			//Options to tell jQuery not to worry about content-type.
			processData: false,
			cache: false,
			contentType: false,
			beforeSend: function (jqXHR, settings) {
			}
			}).done(function (data, textStatus, jqXHR) {
			
			// Now retrieve all links and their descriptions from NAR site to immitate an RSS feed and send them to renderRSSData
			//
			// INFO: in case we alter the RSS source change the global variables RSS_prefix and RSS_link - if the link points to 
			// an RSS feed use get_rss.php (slight changes might be needed in get_rss.php but it will work in most cases)
			//
			// If they point to a simple webpage you can use get_NAR_articles.php - alter the php code as needed, comment in the block below and replace renderRSSData(data, renderelem); with renderRSSData(jsonData, renderelem);:
			// get_NAR_articles.php will retrieve URL links and their descriptions using regexes on the html code
			//
			//
			//{
				// if (analysis_finished || RreturnedErrorsInCurSession) return; //in case the rss delays make sure that the r script has not already returned any errors and that the results are not yet displayed
				// var htmlcode = data.html_code;
				// //get the div where NAR stores their current articles
				// var match =  htmlcode.match(/<strong>([\s\S]+?)<\/a>/gi);
				// Link retrieving START
				// var jsonData = {};
				// $.each(match, function(idx, my_match)
				// {
				// var mytempmatch = my_match.match(/<strong>(.+?)<\/strong>/i);
				// var mytitle = mytempmatch[1];
				// //clear the title from formats
				// mytitle = mytitle.replace(/<.*?>/g, "");
				// var mytempmatch2 = my_match.match(/<a href="(.+?)">/i);
				// var myurl = mytempmatch2[1];
				// jsonData[mytitle] = myurl;
				// });
				// // Link retrieving END
			//}
			//
			//
			
			renderRSSData(data, renderelem);
			}).fail(function (jqXHR, textStatus, errorThrown) {
			// Add failure content here
		});
	}
	
	
	var helper_setItemArrAttr = function (selector, idx, json) {
		
		// From an array of items matched by "selector", choose the one with index "idx" and set some attribute(s)
		$($(selector).get(idx)).attr(json);
	}
	
	
	var helper_setItemAttr = function (selector, json) {
		// For a single item matched by "selector" set some attribute(s)
		
		var x = $(selector).attr(json);
	}
	
	
	var idsToStrings = function (itemsSelector) {
		// Get an array of items IDs' strings (starting with a #)
		
		// Using jQuery syntax $(itemsSelector).toArray() returns an array of all items matching the 
		// itemsSelector string. e.g. if itemsSelector == ".main_div .main_section" (as in PS initialization)
		// $(itemsSelector).toArray() will return an array of all items in main_div or main_section class
		
		//  The ids of these items are then extracted and returned as an array
		
		
		return $.map($(itemsSelector).toArray(), function (itm, i) {
			return "#" + $(itm).attr("id");
		});
	}
	
	
	var InitializeAdvParam = function() {
		
		// Initializes the advanced parameters dialog box
		// the dialog box is used to set the global variables LeastBRep that
		// is the least breps where a protein must be found in order not to be disqualified
		// LeastPep: least peptides where a protein must be found in order not to be disqualified
		// and Pthreshold
		
		// Populate the drop down list for LeastPep and set the current one as selected
		$("#AdvParamLPep").empty();
		for(var i = 1; i < 11; i++)
		{
			if(i == LeastPep)
			{
				$("#AdvParamLPep").append("<option value='" + i + "' selected='true'>" + i + "</option>");
			}
			else{
				$("#AdvParamLPep").append("<option value='" + i + "'>" + i + "</option>");
			}
		}
		
		// Populate the drop down list for LeastBRep and set the current one as selected
		// The possibele values in the drop down list idepends on the Breps of the experiment
		$("#AdvParamLBRep").empty();
		var my_bioreps = [];
		$.each(rawfiles_structure, function (idx, my_raw_file)
			{
				if (my_raw_file.used == false) return;
				my_bioreps.push(my_raw_file.biorep);
			});
			var unique_breps = ArrNoDupe(my_bioreps);
			for(var i = 1; i < unique_breps.length + 1; i++)
			{
				if(i == LeastBRep)
				{
					$("#AdvParamLBRep").append("<option value='" + i + "' selected='true'>" + i + "</option>");
				}
				else{
					$("#AdvParamLBRep").append("<option value='" + i + "'>" + i + "</option>");
				}	
			}
			$("#PThreshold").val(Pthreshold);
	}
	
	
	var InitializeLS = function() {
		
		// Initializes the LS dialog box
		
		if (!AllowLS) return;
		
		// First get all the bioreps and techreps available
		var my_bioreps = [];
		var my_techreps = [];
		$.each(rawfiles_structure, function (idx, my_raw_file)
			{
				// Get all bioreps and techreps to an array fetched from rawfiles_structure
				if (my_raw_file.used == false) return;
				my_bioreps.push(my_raw_file.biorep);
				my_techreps.push(my_raw_file.techrep);
			}); // END for each rawfiles_structure
			
			// Convert the array to a unique array
			var unique_breps = ArrNoDupe(my_bioreps);
			var unique_treps = ArrNoDupe(my_techreps);
			
			// Refresh the dropdown lists that correspond to breps and treps respectively
			$("#LSbrep").empty();
			$("#LStrep").empty();
			// First add a header:
			$("#LSbrep").append("<option value='" + "BioReps" + "' selected='true'>" + "Bio Reps" + "</option>");
			$("#LStrep").append("<option value='" + "TechReps" + "' selected='true'>" + "Tech Reps" + "</option>");
			// Then append the valid breps and treps
			$.each(unique_breps, function (idx, my_brep)
				{
					$("#LSbrep").append("<option value='" + my_brep + "'>" + my_brep + "</option>");
				}); // END for each unique_breps
				$.each(unique_treps, function (idx, my_trep)
					{
						$("#LStrep").append("<option value='" + my_trep + "'>" + my_trep + "</option>");
					});// END for each unique_treps
					
					
					// Now refresh the label drop down lists
					$("#LSfirstlabel").empty();
					$("#LSsecondlabel").empty();
					$("#LSRawFilesList").empty();
					$("#conditions_list option").each(function(idx, my_item){
						
						$("#LSfirstlabel").append("<option value='" + $(this).val() + "'>" + $(this).val() + "</option>");
						$("#LSsecondlabel").append("<option value='" + $(this).val() + "'>" + $(this).val() + "</option>");
					}); // END for eac option in conditions_list
					
					LSselected_raws = [];
	}
	
	
	var InitializeRename = function() {
		// Initializes label renaming (otherwise label merging). See Label Merging in our documentation for more details
		//
		// Label renaming is based on matching 2 or more of the labels to the same condition - this information
		// is stored in RenameArray
		// for example if an experiment contains 3 labels (L, M, H) and L and H both reflect Lung_Ca 
		// then RenameArray should be:
		//
		// L	Lung_Ca
		// M	M
		// H	Lung_Ca
		//
		// AuthenticItemsinRename contain the conditions as they are read now in conditions list
		
		
		RenameArray = [];
		AuthenticItemsinRename = [];
		
		//first get all conditions in conditions list to an array
		
		var conditions_in_list = [];
		$("#conditions_list option").each(function(){ conditions_in_list.push($(this).val());})
		
		// Then for each condition 
		$.each(conditions_in_list, function(idx, my_cond)
			{
				// match each condition to itself so that the user can group the labels later safely
				//
				//e.g. RenameArray may be initialized as:
				//
				// L	L
				// M	M
				// H	H
				//
				
				var value_to_push = [];
				value_to_push[0] = my_cond;
				value_to_push[1] = my_cond;
				RenameArray.push(value_to_push);
				AuthenticItemsinRename.push(value_to_push[0]);
			});
			
			// The user may group labels by right clicking on the conditions list and using the context menu that appears
			// Refresh_conds_list_cmenu_items initializes this menu
			
			Refresh_conds_list_cmenu_items();
	}
	
	
	var inputChCheck = function (e, repatt, maxCharacters) {
		
		// inputChCheck validates that a user-typed text does not
		// contain invalid or too many characters
		// In case it does it prevents Default (i.e. most of the times prevent typing a new character)
		// e: the element, repatt: the regular expression that needs to be matched to consider a text valid
		
		var theEvent = e || window.event;
		var key = theEvent.keyCode || theEvent.which;
		var re = new RegExp(repatt);
		var srcelem = e.target || e.srcElement;
		var carposition = srcelem.selectionStart; // the carret position
		if (srcelem.selectionEnd - srcelem.selectionStart >= 1)
		{
			maxCharacters = 9999;
		}
		var txt = [srcelem.value.slice(0, carposition), String.fromCharCode(e.which), srcelem.value.slice(carposition)].join('');
		if (key === 8) {
			return;
		}
		
		if (!re.test(txt) || txt.length > maxCharacters) {
			e.preventDefault();
		}
	}
	
	
	var inputPasteCheck = function(e, repatt, maxCharacters) {
		// This function makes sure that in case the user pastes data in a text box where character rules apply,
		// the data will contain only valid characters
		// e: the element, repatt: the regular expression that needs to be matched to consider a text valid
		
		var re = new RegExp(repatt);
		var srcelem = e.target || e.srcElement;
		var carposition = srcelem.selectionStart; // the carret position
		var clipboardData, pastedData;
		
		// Prevent default anyway
		e.preventDefault();
		
		// Get the clipboard data
		clipboardData = e.clipboardData || window.clipboardData;
		pastedData = clipboardData.getData('Text'); //the data that the user tries to paste
		
		if (pastedData == "") return;
		
		// Make sure that the data contain only valid characters, in case of invalid chars transform them to underscores
		if (!re.test(pastedData))
		{
			// If it contains only valid characters replace anything that is  to character or number with underscores
			// this prevents errors in R backend
			pastedData = pastedData.replace(/[^a-zA-Z0-9]+/g, "_");
		}
		
		if(pastedData.slice(0,1) == "_")
		{
			pastedData = pastedData.slice(1);
		}
		
		// Get current carret position
		var srcValue = srcelem.value;
		if (srcelem.selectionEnd - srcelem.selectionStart >= 1)
		{
			srcValue = [srcValue.slice(0, srcelem.selectionStart), srcValue.slice(srcelem.selectionEnd)].join('');
		}
		
		// The text that will be produced will be:
		var txt = [srcValue.slice(0, carposition), pastedData, srcValue.slice(carposition)].join('');
		if(txt.length > maxCharacters)
		{
			txt = txt.slice(0, maxCharacters)
		}
		
		// Replace the text immitating the paste procedure
		srcelem.value = txt;
	}
	
	
	var is_RM_ready = function(step) {
		
		// This function returns false if the RM data given by the user are not enough or invalid.
		// if "step" is defined it only validates the corresponding RMstep, Otherwise it validates both RM2 and RM3
		var isready = true;
		if (typeof(step) === 'undefined' || step == 2)
		{
			for (var i = 0; i <= RMrawfilesdata.length - 1; i++)
			{
				// For each RMrawfiledata row
				if (!(RMrawfilesdata[i][6] == 'true')) continue; // If the raw file is not used continue
				// If the breps are represented as different rawfiles and a rawfile is not set return false
				if (RMbrepsRepInRawFiles && RMrawfilesdata[i][2] == "-")
				{
					isready = false;
					break;
				}
				// The same applies for treps and conditions
				if (RMtrepsRepInRawFiles && RMrawfilesdata[i][3] == "-")
				{
					isready = false;
					break;
				}
				if (RMconditionsRepInRawFiles && RMrawfilesdata[i][5] == "-")
				{
					isready = false;
					break;
				}
			}
		}
		if (typeof(step) === 'undefined' || step == 3)
		{
			for (var i = 0; i <= RMtagsdata.length - 1; i++)
			{
				// For each RMtagsdata row
				if (!(RMtagsdata[i][6] == 'true')) continue; //if the tag is not used continue
				// If the breps are represented as different tags and a tag is not set return false
				if (!RMbrepsRepInRawFiles && RMtagsdata[i][2] == "-")
				{
					isready = false;
					break;
				}
				// The same applies for treps and conditions
				if (!RMtrepsRepInRawFiles && RMtagsdata[i][3] == "-")
				{
					isready = false;
					break;
				}
				if (!RMconditionsRepInRawFiles && RMtagsdata[i][5] == "-")
				{
					isready = false;
					break;
				}
			}
		}
		return isready;
	}
	
	
	var LoadParams = function(myparametersstring) {
		
		// Load params gets the parameter file as a multiline string (myparametersstring) and loads ProteoSign's variables
		// one by one (see CreateParamsFile for parameters file structure)
		
		// CheckParamsValidity returns "" if all parameters in the parameters file are valid and an html block
		// describing the errors otherwise
		
		var checkres = CheckParamsValidity(myparametersstring);
		if (checkres == "")
		{
			// If all parameters are valid:
			var additional_info = "";
			
			// Split all lines and check that the first one contains a valid comment
			// (not the best validity check but makse sure that no random file is uploaded as parameter file)
			
			var paramslines = myparametersstring.split("\n");
			var setvar = "";
			var mypatt = new RegExp("#ProteoSign parameters files");
			var myres = mypatt.test(paramslines[0]);
			if (!myres)
			{
				msgbox("The uploaded file is not valid")
				return;
			}
			
			$.each(paramslines, function(idx, myparamline)
				{
					// For each parameter file line:
					// Get rid of \r characters to compensate with windows-edited parameter files
					if(myparamline[myparamline.length - 1] == "\r")
					{
						myparamline = myparamline.substring(0, myparamline.length - 1);
					}
					
					// Get rid of comments
					if (myparamline.charAt(0) == "#")
					{
						return;
					}
					
					// setvar contains the name of the variable to be set
					// if setvar is blank then the next line should be a variable name
					// most of these are self explanatory
					if (setvar == "")
					{
						if (myparamline == "isLabelFree" || myparamline == "isIsobaricLabel" || myparamline == "procprogram" || myparamline == "rawfiles_structure" || myparamline == "expid" || myparamline == "exptpoint" || myparamline == "conditions_to_compare" || myparamline == "quantitation_filtering" || myparamline == "filtering_label" || myparamline == "peptide_level_filtering" || myparamline == "LeastBRep" || myparamline == "LeastPep" || myparamline == "Pthreshold" || myparamline == "LFQconditions" || myparamline == "!Rename" || myparamline == "!LS_Array" || myparamline == "!LS_c_p_Add" || myparamline == "!Select_Labels" || myparamline == "!RMrawfilesdata" || myparamline == "!RMtagsdata" || myparamline == "!RMbrepsRepInRawFiles" || myparamline == "!RMtrepsRepInRawFiles" || myparamline == "!RMconditionsRepInRawFiles" || myparamline == "!RMisused")
						{
							setvar = myparamline;
						}
					}
					else
					{
						// setvar contains the PS variable that has to be set
						// The following lines set the corresponding variable or alter a DOM element value depending on the setvar variable
						if (setvar == "rawfiles_structure")
						{
							// Here add_raw_file_structure is used to make sure that all rawfiles in the dataset are included in the parameter's file
							if (add_raw_file_structure(myparamline, false) == false)
							{
								
								additional_info = "<br>Some raw files where not matched to a biological coordinate"
								if (my_step == 3)
								{ // If we are on step 3 send the user back to step 2 to take care of the unset rawfiles
									toggleNextClass(idsToStrings(".main_div .main_section").reverse(), "hidden", true, rollbackStage);
									rawfiles_tbl_allfiles_DT.columns.adjust().draw();
								}
							}
						}
						// Experiment information
						else if(setvar == "expid")
						{
							$('input[name="expid"]').val(myparamline);
						}
						else if (setvar == "exptpoint")
						{
							$('input[name="exptpoint"]').val(myparamline);
						}
						// Conditions that will be used for comparison
						else if (setvar == "conditions_to_compare")
						{
							// Unfold one line representation of AuthenticItemsinRename and set the variable
							var myvalues = myparamline.split("\t");
							AuthenticItemsinRename = [];
							$.each(myvalues, function(idx, myval){
								if (myval == "") return;
								AuthenticItemsinRename.push(myval);
							});
						}
						else if (setvar == "quantitation_filtering")
						{
							if (myparamline == "F")
							{
								$("#expquantfilt").prop("checked", false);
								//hide the rest
								$("#s3advparams select[name='expquantfiltlbl']").closest("tr").addClass("hidden");
								$("#s3advparams input[name='expquantfiltprot']").closest("tr").addClass("hidden");
							}
							else if (myparamline == "T")
							{
								$("#expquantfilt").prop("checked", true);
								//show them again
								$("#s3advparams select[name='expquantfiltlbl']").closest("tr").removeClass("hidden");
								$("#s3advparams input[name='expquantfiltprot']").closest("tr").removeClass("hidden");
								var parentdiv = $(this).closest("div");
								parentdiv.scrollTop(parentdiv.prop("scrollHeight"));
							}
						}
						// Peptide filtering variables
						else if (setvar == "filtering_label")
						{
							$('#expquantfiltlbl option:contains("' + myparamline + '")').prop('selected', true);
						}
						else if (setvar == "peptide_level_filtering")
						{
							if (myparamline == "F")
							{
								$("#expquantfiltprot").prop("checked", false);
							}
							else if (myparamline == "T")
							{
								$("#expquantfiltprot").prop("checked", true);
							}
						}
						// Advanced parameters
						else if (setvar == "LeastBRep")
						{
							LeastBRep = myparamline;
						}
						else if (setvar == "LeastPep")
						{
							LeastPep = myparamline;
						}
						else if (setvar == "Pthreshold")
						{
							Pthreshold = myparamline;
						}
						// Label free conditions are assigned to all rawfiles in table in step 2
						// and add_LFQ_conditions takes care of setting RawFileConditions array
						else if (setvar == "LFQconditions")
						{
							if (myparamline == "")
							{
								RawFileConditions = [];
								//Show everything back to the user
								var all_items = $('#rawfiles_tbl_allfiles').find('tr');
								for (var i = 0; i < all_items.length; i++) {
									var items_tds = $(all_items[i]).find('td');
									var items_cond = items_tds[4];
									$.each(RawFileConditions, function (idx, my_raw_file)
										{
											if (my_raw_file.name == $(items_tds[0]).text())
											{
												$(items_cond).text(my_raw_file.condition);
												return false;
											}
										});
								}
								setvar = "";
								return;
							}
							add_LFQ_conditions(myparamline);
						}
						// Rename - merge conditions
						else if (setvar == "!Rename")
						{
							// If merging is disabled or the parameter is blank simply erase data in Rename related variables
							if (myparamline == "")
							{
								RenameArray = [];
								AuthenticItemsinRename = [];
								setvar = "";
								return
							}
							if(!AllowMergeLabels)
							{
								RenameArray = [];
								AuthenticItemsinRename = [];
								setvar = "";
								return
							}
							// Initialize the arrays
							RenameArray = [];
							AuthenticItemsinRename = [];
							var my_val = myparamline;
							// Unfold the parameter line and assign the corresponding values to RenameArray
							// The first element in each row in RenameArray is an element in AuthenticItemsinRename
							var lines = my_val.split("||");
							$.each(lines, function(idx, my_line){
								if (my_line == "") return;
								var my_values = my_line.split("|");
								RenameArray.push(my_values);
								AuthenticItemsinRename.push(my_values[0]);
							});
							
							// Run Refresh_conds_list to refresh the conditions list based on Renaming arrays
							// Temporarily set RenameFromTestData to true so that Refresh_conds_list does not worry about
							// expquantfiltlbl contents since this will be automatically loaded from the parameter file
							var temprenamefromtestdata = RenameFromTestData;
							RenameFromTestData = true;
							Refresh_conds_list();
							RenameFromTestData = temprenamefromtestdata;
						}
						else if (setvar == "!LS_Array")
						{
							// Set the label swap arrays
							
							// If the parameter is blank simply clear LabelSwapArray
							if (myparamline == "")
							{
								LabelSwapArray = [];
								setvar = "";
								return
							}
							if(!AllowLS) return;
							// Unfold the parameter and set the LabelSwapArray
							LabelSwapArray = [];
							var my_val = myparamline;
							var lines = my_val.split("||");
							$.each(lines, function(idx, my_line){
								if (my_line == "") return;
								var my_values = my_line.split("|");
								var to_add = [];
								to_add.push(my_values[0]);
								to_add.push(my_values[1]);
								to_add.push(my_values[2]);
								to_add.push(parseInt(my_values[3]));
								LabelSwapArray.push(to_add);
							}); // END for each line
							LS_array_counter = lines.length;
						}
						else if (setvar == "!LS_c_p_Add")
						{
							// Set LS_counters_per_Add
							
							// If the parameter is blank simply clear LS_counters_per_Add
							if (myparamline == "")
							{
								LS_counters_per_Add = [];
								$("#LSLabelSwaps").empty();
								setvar = "";
								return
							}
							if(!AllowLS) return;
							LS_counters_per_Add = [];
							var my_val = myparamline;
							// Unfold the parameter and set LS_counters_per_Add
							var lines = my_val.split("|||");
							$("#LSLabelSwaps").empty();
							$.each(lines, function(idx, my_line){
								if (my_line == "") return;
								var my_values = my_line.split("||");
								var to_add = [];
								to_add[0] = my_values[0];
								to_add[1] = [];
								var seclines = my_values[1].split("|");
								$.each(seclines, function(idx, my_secline){
									to_add[1].push(parseInt(my_secline));
								}); // END for each second line
								to_add[2] = [];
								to_add[2].push(my_values[2]);
								to_add[2].push(my_values[3]);
								LS_counters_per_Add.push(to_add);
								//Add the respective lines in Label Swaps list
								$("#LSLabelSwaps").append("<option value='" + my_values[0] + "'>" + my_values[0] + "</option>");
							}); // END for each line
						}
						else if (setvar == "!Select_Labels")
						{
							// Select specific conditions from the condition list
							if (myparamline == "")
							{
								setvar = "";
								return
							}
							// Emulate a condition where specific labels are selected from a test dataset
							// to make the procedure faster
							my_lbls_toselect = myparamline.split("|");
							select_labels_according_to_test_dataset();
							my_lbls_toselect = [];
						}
						// The following blocks use string_to_array to unfold the respective parameters and set
						// RM related variables
						else if (setvar == "!RMrawfilesdata")
						{
							RMrawfilesdata = string_to_array(myparamline);
						}
						else if (setvar == "!RMtagsdata")
						{
							RMtagsdata = string_to_array(myparamline);
						}
						else if (setvar == "!RMbrepsRepInRawFiles")
						{
							RMbrepsRepInRawFiles = (myparamline == 'true');
						}
						else if (setvar == "!RMtrepsRepInRawFiles")
						{
							RMtrepsRepInRawFiles = (myparamline == 'true');
						}
						else if (setvar == "!RMconditionsRepInRawFiles")
						{
							RMconditionsRepInRawFiles = (myparamline == 'true');
						}
						else if (setvar == "!RMconditionsRepInRawFiles")
						{
							RMconditionsRepInRawFiles = (myparamline == 'true');
						}
						else if (setvar == "!RMisused")
						{
							// Sets RMisused and uses set_RMisused to alter the UI and some RM related variables
							// depending on if RM is used or not
							// RM_init_RMsteps_from_load_data initializes RM if RM is to be used
							RMisused = (myparamline == 'true');
							if (RMisused)
							{
								set_RMisused(true);
								RM_init_RMsteps_from_load_data();
								reset_reps();
							}
							else
							{
								set_RMisused(false, false);
							}
						}
						setvar = "";
					}
				});
				// Display a success message
				msgbox("Loading completed." + additional_info);
		}
		else
		{
			// Display an error message
			msgbox("<p style='border-bottom: solid; text-align: center'><strong>The following errors occured when loading parameters:</strong></p><ol>" + checkres + "</ol>");
		}
	}
	
	
	var log_php = function(texttoappend) {
		
		// This function appends the specified text to php log file of the current session
		// alter log_php.php to redirect these logs to another location
		// Notice that this function has no error handling code
		
		var thedata = new FormData();
		thedata.append('texttoappend', texttoappend);
		thedata.append('session_id', sessionid);
		$.ajax({
			url: cgi_bin_path + 'log_php.php', //Server script to send the feedback
			type: 'POST',
			// Form data
			data: thedata,
			//Options to tell jQuery not to worry about content-type.
			processData: false,
			cache: false,
			contentType: false,
		});
	}
	
	
	function longRawFileNameDisplay() {
		
		// When the mouse hovers over a raw file that does not fit the cell of Stage 2 table
		// longRawFileNameDisplay shws a tooltip with the full name of the rawfile
		
		var $this = $(this);
		if (this.offsetWidth < this.scrollWidth && !$this.attr('title')) {
			$this.attr('title', $this.text());
		}
	}
	
	
	var msgbox = function(displaytext) {
		
		// Displays a simple message in a box
		$(".expparamsDlgInfo").css({"left": ($("body").width() / 2) - ($("#VariousInfoDialog").width() / 2)});
		$('body').append('<div id="maskInfo"></div>');
		$("#VariousInfoDialog").fadeIn(300);
		$('#maskInfo').fadeIn(300);
		$("#InfoDiv").empty();
		$("#InfoDiv").append("<span>" + displaytext + "</span>");
	}
	
	
	var nextStageIsFirstStage = function (lastStageItem) {
		
		// Called by "toggleNextClass". It is used for hiding the "Previous" button when we are at stage 1
		//i.e. there is no previous stage.
		
		$("#s1btnb").addClass("hidden");
	}
	
	
	var onADVoptionsCancel_click = function() {
		
		// Bound to cancel button on advanced parameters dialog
		setItemAttrColor("#PThreshold", "border", "#DDDDDD");
		dlgFadeout();
	}
	
	
	var onADVoptionsOK_click = function() {
		
		// Bound to click event of ADVoptionsOK, the OK button in Advanced parameters dialog
		// sets LeastPep LeastBRep and pThreshold according to the options of the user
		
		// Replace , with . to compensate with decimal seperator issues
		$("#PThreshold").val($("#PThreshold").val().replace(",", "."));
		
		// Always make sure that the value is valid before assigning it
		if (StringisNumber($("#PThreshold").val()))
		{
			if ($("#PThreshold").val() > 0 && $("#PThreshold").val() < 1)
			{
				Pthreshold = $("#PThreshold").val();
			}
			else
			{
				setItemAttrColor("#PThreshold", "border", "#E60000");
				return;	
			}
		}
		else{
			setItemAttrColor("#PThreshold", "border", "#E60000");
			return;
		}
		LeastPep = $("#AdvParamLPep").val();
		LeastBRep = $("#AdvParamLBRep").val();
		setItemAttrColor("#PThreshold", "border", "#DDDDDD");
		dlgFadeout();
	}
	
	
	var onAssignCondition = function(){
		
		// Bound to "Assign Condition" item on rawfiletable's context menu to
		// display the assign condition dialog box
		
		$(".expparamsDlg").css({"left" : ($("body").width()/2) - ($("#s2LFQConditions").width()/2)});
		$('body').append('<div id="mask"></div>');
		$("#s2LFQConditions").fadeIn(300);
		$('#mask').fadeIn(300);
	}
	
	
	var onCarefulBackNoclick = function() {
		
		// Bound to No button on carefulback dialog box, this function simply closes the dialog box
		// leaving the user to the current step
		dlgFadeout();
	}
	
	
	var onCarefulBackYesclick = function() {
		
		// Bound to Yes button on carefulback dialog box, this function shows the previous step if the user decided
		// to move a step back despite the warning message
		
		toggleNextClass(idsToStrings(".main_div .main_section").reverse(), "hidden", true, rollbackStage);
		rawfiles_tbl_allfiles_DT.columns.adjust().draw();
		// Close the popup window
		dlgFadeout();
	}
	
	
	var onChooseFromTestDatasets = function(){
		
		// Bound to "choose" text in Stage 1 to display the test dataset dialog box
		
		if (!datatestOK_clicked)
		{
			$(".expparamsDlg").css({"left" : ($("body").width()/2) - ($("#s1TestDatasets").width()/2)});
			$('body').append('<div id="mask"></div>');
			$("#s1TestDatasets").fadeIn(300);
			$('#mask').fadeIn(300);
		}
	}
	
	
	function onclickQuantitationFilteringChkBox() {
		
		// Toggles visibility of quantitation filtering options based on whether the Quantitation Filtering check box
		// is checked or not (bound to its click event)
		
		if ($(this).prop('checked')) {
			$("#s3advparams select[name='expquantfiltlbl']").closest("tr").removeClass("hidden");
			$("#s3advparams input[name='expquantfiltprot']").closest("tr").removeClass("hidden");
			var parentdiv = $(this).closest("div");
			parentdiv.scrollTop(parentdiv.prop("scrollHeight"));
		}
		else
		{
			$("#s3advparams select[name='expquantfiltlbl']").closest("tr").addClass("hidden");
			$("#s3advparams input[name='expquantfiltprot']").closest("tr").addClass("hidden");
		}
	}
	
	
	var oncondslistclick = function(){
		
		// Refresh conditions list cintext menu items every time it is being (right) clicked
		if (AllowMergeLabels)
		{
			Refresh_conds_list_cmenu_items();	
		}
	}
	
	
	var onexpquantfiltlblclick = function() {
		
		// Erase possible red border style from expquantfiltlbl that might had been added due to validation failure
		
		setItemAttrColor("#expquantfiltlbl", "border", "#DDDDDD");
	}
	
	
	var onFeedbackclick = function() {
		
		// Show the feedback panel and bind a simple cancel function to the cancel button
		
		$(".expparamsDlgFeedback").css({"left": ($("body").width() / 2) - ($("#s4FeedbackPanel").width() / 2)});
		$('body').append('<div id="maskFeedback"></div>');
		$("#s4FeedbackPanel").fadeIn(300);
		$('#maskFeedback').fadeIn(300);
		
		$("#FeedbackPanelCancel").unbind();
		$("#FeedbackPanelCancel").on("click", function () {
			dlgFadeoutFeedback();
			$("#userFeedback").val("");
			var mytext = 600 + " characters left";
			$("#FBcharleft").text(mytext);
			//ADD CANCEL RESULT HERE
		});
	}
	
	
	var onFeedbackPanelSendclick = function() {
		// This function uses send_feedback.php to send the user's feedback to ProteoSign's creators
		
		// First check that the user has written some feedback text or that has chosen a star rating
		if (parseInt($("#userFeedback").val().length) == 0 && !stars_enabled)
		{
			msgbox("<p>Please fill the text box with your feedback.</p>");
			return;
		}
		
		// send feedback text and star rating to send_feedback.php
		var thedata = new FormData();
		thedata.append('texttoappend', $("#userFeedback").val());
		thedata.append('session_id', sessionid);
		if (!stars_enabled)
		{
			thedata.append('stars_score', "NA");
		}
		else
		{
			thedata.append('stars_score', star_chosen);
		}
		thedata.append('feedback_file', server_Feedback_file);
		$.ajax({
			url: cgi_bin_path + 'send_feedback.php', //Server script to send the feedback
			type: 'POST',
			// Form data
			data: thedata,
			//Options to tell jQuery not to worry about content-type.
			processData: false,
			cache: false,
			contentType: false,
			}).done(function (data, textStatus, jqXHR) {
			// If the connection finished successfully displa a thank you message
			if (data.success == true)
			{
				msgbox("<p>Thank you! Your feedback has been submitted successfully!</p>");
				dlgFadeoutFeedback();
				$("#userFeedback").val("");
				var mytext = 600 + " characters left";
				$("#FBcharleft").text(mytext);
			}
			else{
				msgbox("<p>An error occured! Please try again.</p>");
			}
			}).fail(function (jqXHR, textStatus, errorThrown){
			msgbox("<p>An error occured! Please try again.</p>");
		});
	}
	
	
	var onlistclick = function() {
		
		// Erase possible red border style from conditions list that might had been added due to validation failure
		
		setItemAttrColor("#conditions_list", "border", "#DDDDDD");
	}
	
	
	var onLSAddclick = function() {
		
		// onLSAddclick is bound to Add button in LS dialog box
		// it first runs some validation tests to make sure that the label swap is valid
		// and if it is it resets the form, updates the Current label swaps list
		// and refreshes LabelSwapArray and LS_counters_per_Add
		
		if(!AllowLS) return;
		
		// First validate that the labels to be swapped are not the same and 
		var proc_failed = false;
		var myfirstlabel = $("#LSfirstlabel").val();
		var myseclabel = $("#LSsecondlabel").val();
		if (myfirstlabel == myseclabel)
		{
			setItemAttrColor("#LSfirstlabel", "border", "#E60000");
			setItemAttrColor("#LSsecondlabel", "border", "#E60000");
			proc_failed = true;
		}
		
		// Also make sure that at least one raw file will be affected
		if (LSselected_raws.length < 1 )
		{
			setItemAttrColor("#LSbrep", "border", "#E60000");
			setItemAttrColor("#LStrep", "border", "#E60000");
			proc_failed = true;
		}
		if (proc_failed) return;
		
		// Now validate that the user has not lready set the same Label swap
		
		// In label swap, label order does not matter, e.g. if the user has already set a label swap of L vs H
		// then setting a label swap of H vs L (for the same breps and treps) is invalid
		
		// First create a simple description of the label swap
		var my_desc = myfirstlabel + " - " + myseclabel + " in b" + $("#LSbrep").val() + "t" + $("#LStrep").val();
		
		// Also create a shifted description (that is a label swap swapping the second with the first element)
		var shifted_desc = myseclabel + " - " + myfirstlabel + " in b" + $("#LSbrep").val() + "t" + $("#LStrep").val();
		
		// In case Label swaps list contains either current Label Swap's description or the shifted one
		// call it a duplicate and abort the procedure
		
		var found_duplicate = false;
		$("#LSLabelSwaps option").each(function(idx, my_item){	
			if ($(this).val() == my_desc || $(this).val() == shifted_desc)
			{
				$(this).prop("selected", true);
				found_duplicate = true;
			}
			else
			{
				$(this).prop("selected", false);
			}
		}); // END for each option of LSLabelSwaps
		
		if (found_duplicate)
		{
			setItemAttrColor("#LSLabelSwaps", "border", "#E60000");
			return;
		}
		
		
		// Refresh LabelSwapArray and LS_counters_per_Add
		//
		// LabelSwapArray is in fact a matrix (an array of arrays) that has 4 columns: rawfile, firstlabel, secondlabel and counter (a key)
		// Each raw corresponds to a different raw file where a label swap of firstlabel and secondlabel occur
		// e.g.:
		//
		//
		// +-------------------------------+-------------+--------------+---------------+
		// |		   Raw File			   | First Label | Second Label | Counter (key) |
		// +-------------------------------+-------------+--------------+---------------+
		// | OT2_Terhune(...)DMC-HLNF01_02 | L		     | H			|			 0  |
		// | OT2_Terhune(...)DMC-HLNF02_02 | L		     | H			|			 1  |
		// | OT2_Terhune(...)DMC-HLNF03_02 | L		     | H			|			 2  |
		// | OT2_Terhune(...)DMC-HLNF04_02 | L		     | H			|			 3  |
		// | (...)						   |		     |			    |			    |
		// | OT2_Terhune(...)DMC-HLNF12_02 | L		     | H			|			11  |
		// +-------------------------------+-------------+--------------+---------------+
		//
		//
		// This matrix can be sent to the backend of PS to perform Label Swapping
		// For more on Label Swapping see our documentation
		
		var ArrayOfIndices = [];
		$.each(LSselected_raws, function (idx, my_raw_file)
			{
				ArrayOfIndices.push(LS_array_counter);
				var toadd = [];
				toadd[0] = my_raw_file;
				toadd[1] = myfirstlabel;
				toadd[2] = myseclabel;
				toadd[3] = LS_array_counter++;
				LabelSwapArray.push(toadd);
			}); // END for each LSselected_raws
			
			// Also refresh LS_counters_per_Add - each row of this array corresponds to one assigned label swap
			// (and to one row in Label Swaps list in LS Dialog box) and has three items
			// the first is Label Swap's description as shown on the list
			// the second is an array of keys to match the label swap to Rows in LabelSwapArray
			// the last one is an array of the two labels taht are swapped in this label swap
			//e.g.
			//
			// +---------------+-------+--------+
			// | Description   | Keys  | Labels |
			// +---------------+-------+--------+
			// | L - H in b1t2 | 1	   | L	    |
			// |			   +-------+--------+
			// |			   | 2	   | H	    |
			// |			   +-------+--------+
			// |			   | 3	   |		|
			// |			   +-------+		|
			// |			   | 4	   |		|
			// |			   +-------+		|
			// |			   | 5	   |		|
			// |			   +-------+		|
			// |			   | (...) |		|
			// |			   +-------+		|
			// |			   | 11	   |		|
			// +---------------+-------+--------+
			//
			//
			var to_add2 = [];
			to_add2[0] = my_desc;
			to_add2[1] = ArrayOfIndices;
			to_add2[2] = [myfirstlabel, myseclabel];
			LS_counters_per_Add.push(to_add2);
			
			// Finally refresh the LSLabelSwaps list:
			$("#LSLabelSwaps").append("<option value='" + my_desc + "'>" + my_desc + "</option>");
			setItemAttrColor("#LSLabelSwaps", "border", "#DDDDDD");
	}
	
	
	var onLSbackclick = function() {
		// Bound to back button in LS dialog box, onLSbackclick simply hides the dialog box
		// and resets its elements border colors to black in case
		// one was marked red due to input error
		//
		// Reset the border colors back to their initial value
		
		setItemAttrColor("#LSfirstlabel", "border", "#DDDDDD");
		setItemAttrColor("#LSsecondlabel", "border", "#DDDDDD");
		setItemAttrColor("#LSbrep", "border", "#DDDDDD");
		setItemAttrColor("#LStrep", "border", "#DDDDDD");
		setItemAttrColor("#LSLabelSwaps", "border", "#DDDDDD");
	}
	
	
	var onLSlabelchange = function() {	
		
		// Bound to change event of the Label dropdown lists in LS dialog box
		// onLSlabelchange verifies that the two selected labels are not identical.
		// If they are (something invalid) it marks the border with red
		
		if(!AllowLS) return;
		var myfirstlabel = $("#LSfirstlabel").val();
		var myseclabel = $("#LSsecondlabel").val();
		if (myfirstlabel == myseclabel)
		{
			setItemAttrColor("#LSfirstlabel", "border", "#E60000");
			setItemAttrColor("#LSsecondlabel", "border", "#E60000");
		}
		else
		{
			setItemAttrColor("#LSfirstlabel", "border", "#DDDDDD");
			setItemAttrColor("#LSsecondlabel", "border", "#DDDDDD");
		}
	}
	
	
	var onLSLabelSwapsclick = function() {
		
		// Bound to LS slist on LS dialog Box onLSLabelSwapsclick resets its border color to black
		
		if(!AllowLS) return;
		setItemAttrColor("#LSLabelSwaps", "border", "#DDDDDD");
	}
	
	
	var onLSRemoveclick = function() {
		
		// Bound on the remove button in LS Dialog Box, onLSRemoveclick Removes the selected Label Swap
		// by altering the structure of LabelSwapArray and LS_counters_per_Add
		// To learn more about LabelSwapArray and LS_counters_per_Add see the comments on onLSAddclick
		
		if(!AllowLS) return;
		
		var my_keys = []; // The keys of entries in LabelSwapArray and LS_counters_per_Add to be removed
		var options_to_erase = []; // the description of Label Swaps to be removed
		
		
		$("#LSLabelSwaps option:selected").each(function(idx, my_opt){
			$.each(LS_counters_per_Add, function(idx, LS_counters_val)
				{
					if(LS_counters_val[0] == $(my_opt).val()) // If the description in LS_counters_per_Add matches the text in the selected item in LS list
					{
						// Here a Label Swap to be erased is found
						options_to_erase.push($(my_opt).val());
						$.each(LS_counters_val[1], function(idx, my_cval)
							{
								my_keys.push(my_cval);
							}); // END for each key of the selected LS
					}
				}); // END for each LS_counters_per_Add entry
		}); // END for each selected option on LSLabelSwaps
		
		// Remove all the entries in LS_counters_per_Add that contain the key to be removed
		for(var i = LS_counters_per_Add.length - 1; i >= 0;i--)
		{
			if($.inArray(LS_counters_per_Add[i][0], options_to_erase) != -1)
			{
				LS_counters_per_Add.splice(i, 1);
			}
		}
		
		// Remove all the entries in LabelSwapArray that contain the key to be removed
		for(var i = LabelSwapArray.length - 1; i >= 0;i--)
		{
			if($.inArray(LabelSwapArray[i][3], my_keys) != -1)
			{
				LabelSwapArray.splice(i, 1);
			}
		}
		
		
		// Now refresh the lslabelswaps list
		// First find the list options that will remain in display
		options_that_remain = [];
		$("#LSLabelSwaps option").each(function(idx, my_option)
			{
				if($.inArray($(my_option).val(), options_to_erase) == -1)
				{
					options_that_remain.push($(my_option).val());
				}
			}); // END for each LSLabelSwaps option
			
			// Empty the list
			$("#LSLabelSwaps").empty();
			
			// Repopulate it with the options that should still be displayed
			
			$.each(options_that_remain, function(idx, my_option)
				{
					$("#LSLabelSwaps").append("<option value='" + my_option + "'>" + my_option + "</option>");
				}); // END for each options_that_remain
				
				setItemAttrColor("#LSLabelSwaps", "border", "#DDDDDD");
	}
	
	
	var onLSrepchange = function() {
		
		// Bound to both trep and brep dropdown lists value change event in LS dialog box
		// onLSrepchange updates the list of raw files (LSRawFilesList) in the same dialog box
		// Depending on breps and treps that were chosen from the drop down lists, onLSrepchange
		// fetches all rawfiles that are assigned these reps and displays them to the corresponding list
		
		if(!AllowLS) return;
		//get all raw files corresponding the selected brep and trep
		LSselected_raws = [];
		var my_selectedbrep = $("#LSbrep").val();
		var my_selectedtrep = $("#LStrep").val();
		$("#LSRawFilesList").empty();
		$.each(rawfiles_structure, function (idx, my_raw_file)
			{
				if (my_raw_file.used == false) return;
				if (my_raw_file.biorep == my_selectedbrep && my_raw_file.techrep == my_selectedtrep)
				{
					LSselected_raws.push(my_raw_file.rawfile);
				}		 
			}); // END for each rawfiles_structure
			
			// Populate LSRawFilesList with the Label Swap rawfiles that will be affected
			if(!(LSselected_raws.length < 1))
			{
				$.each(LSselected_raws, function (idx, my_raw_file)
					{
						$("#LSRawFilesList").append("<option value='" + my_raw_file + "'>" + my_raw_file + "</option>");
					});  // END for each LSselected_raws
			}
	}
	
	
	var onLSrepclick = function() {
		
		// Bound to brep and trep dropdown list in LS dialog box, onLSrepclick resets
		// their border color to black (see onLSlabelchange for more information)
		
		if(!AllowLS) return;
		setItemAttrColor("#LSbrep", "border", "#DDDDDD");
		setItemAttrColor("#LStrep", "border", "#DDDDDD");
	}
	
	
	var onRMDialogBack_click = function() {
		
		// Simply decrease RM_RMstep by one and run showstepRMDialog to show the corresponding RMstep
		if (RM_RMstep != 1) RM_RMstep--;
		showstepRMDialog();
	}
	
	
	var onRMDialogCancel_click = function () {
		
		// Bound to the cancel button in Replication Multiplexing dialog
		// in case RM is not ready ask the user
		// if they want to discard RM or not
		// if RM is ready simply close the dialog
		
		if (RMisused && is_RM_ready() == false)
		{
			questionbox("The data provided for <i>Replication Multiplexing</i> are not enough. Do you want to discard <i>Replication Multiplexing</i>?", function()
				{
					// postYes
					set_RMisused(false, false);
					dlgFadeoutRepMult();
				},
				function()
				{
					// postNo
				});
		}
		else
		{
			dlgFadeoutRepMult();
		}
	}
	
	
	var onRMDialogNext_click = function() {
		
		// Replication Multiplexing (abbreviated as RM) is a PS feature dedicated to isobarically labelled experiments that allows the user to
		// represent either breps or treps or conditions with different tags or raw files (MS runs)
		// This is a technique used to lower technical errors 
		// For more information on RM please refer to our documentation
		//
		// PS uses the RM dialog to ask the user all necessary information
		// RM is broken into 3 steps (RM1 RM2 and RM3) that are all happening in the RM dialog
		// In RM1 the user chooses how breps, treps and conditions will be reresented (either as different tags or as different raw files - MSruns)
		// In RM2 the user sees a table looking like the rawfile table in step 2 and assigns all attributes represented by different raw files
		// the table is initialized in RMDataTableInit and its right click context menu in RMCMenuInit
		// In RM3 the user assigns all attributes represented by different tags to the corresponding tags using a similar table
		// Notice that even if RM steps are 3, in RM dialog there are simply 2 (stepRM1 stepRM2) since in fact the tables
		// in RM2 and RM3 is the same (RMrawfilesDT)
		//
		// The behaviour of the next button click event depends on whether RM is used for the first
		// time in this PS session and on if the user has used RM again but with different options in RM1
		// if at least one if the above conditions is met (re-)initialize RM2 and RM3:
		
		if ((stepRM2initialized || stepRM3initialized) && RM_RMstep == 1)
		{
			// Here, RM has been used again and the user hit the next button having possibly chosen different options
			// in RM1
			
			// First determine how are different attributes represented:
			var newRMbrepsRepInRawFiles = $('input[name=RMbreprepres]:checked').val() == "rawfiles";
			var newRMtrepsRepInRawFiles = $('input[name=RMtreprepres]:checked').val() == "rawfiles";
			var newRMconditionsRepInRawFiles = $('input[name=RMconditionrepres]:checked').val() == "rawfiles";
			
			// If at least one of the options in RM1 are different from the old ones, ask the user if
			// they want to reset their RM options and revise RM2 and RM3, if so reinitialize steps RM2 and RM3
			
			if (newRMbrepsRepInRawFiles != RMbrepsRepInRawFiles || newRMtrepsRepInRawFiles != RMtrepsRepInRawFiles || newRMconditionsRepInRawFiles != RMconditionsRepInRawFiles)
			{
				if ($('input[name=RMbreprepres]:checked').val() == "rawfiles" && $('input[name=RMtreprepres]:checked').val() == "rawfiles" && $('input[name=RMconditionrepres]:checked').val() == "rawfiles")
				{
					// Always the tags mus represent at least 1 attribute
					msgbox("Something must be represented in your experiment's tags");
					RM_RMstep = 1;
					return;
				}
				// Ask the user if they want to revise RM2 and RM3 with different RM1 options
				questionbox("<p>You changed your options in this step. <strong>ProteoSign</strong> will clear any data pertaining to your experimental structure to setup new forms for you<br>Do you want to proceed?</p>", function(){
					//post-yes function
					stepRM2initialized = false;
					stepRM3initialized = false;
					// Advance one step
					RM_RMstep++;
					showstepRMDialog();
					}, function(){
					// postno function
					// Reset RM1 options
					if (RMbrepsRepInRawFiles)
					{
						$('input[name=RMbreprepres][value=rawfiles]').prop('checked', 'checked');
					}
					else
					{
						$('input[name=RMbreprepres][value=tags]').prop('checked', 'checked');
					}
					if (RMtrepsRepInRawFiles)
					{
						$('input[name=RMtreprepres][value=rawfiles]').prop('checked', 'checked');
					}
					else
					{
						$('input[name=RMtreprepres][value=tags]').prop('checked', 'checked');
					}
					if (RMconditionsRepInRawFiles)
					{
						$('input[name=RMconditionrepres][value=rawfiles]').prop('checked', 'checked');
					}
					else
					{
						$('input[name=RMconditionrepres][value=tags]').prop('checked', 'checked');
					}
				});
				return;
			}
		}
		
		// If the conditions above are not met advance one step
		RM_RMstep++;
		showstepRMDialog();
	}
	
	
	var onRMreview_click = function() {
		
		// Bound on "Review" RM text on the interactove text displayed in Step 2 in case RM is used
		showRepMultDialog();
	}
	
	
	var onRMUndo_click = function() {
		
		// Bound to "Undo" text on the interactove text displayed in Step 2 in case RM is used
		set_RMisused(false);
	}
	
	
	var ons22btnfclick = function() {
		
		// when Next Button in stage 2 is prressed 
		if (isLabelFree)
		{
			// Get a copy of the RawFileConditions array and erase all not used
			RawFileConditions_copy = RawFileConditions.slice();
			$.each(rawfiles_structure, function (idx, my_raw_file)
				{
					if (my_raw_file.used == false)
					{// For each rawfile not used
						$.each(RawFileConditions_copy, function (idxC, my_cond)
							{
								if (my_raw_file.rawfile == my_cond.name)
								{// If there is a row in RawFileConditions_copy corresponding to a not used rawfile, erase this row
									RawFileConditions_copy.splice(idxC, 1);
									// Normally there can not be two rows in RawFileConditions_copy so escape the loop if erased a row
									// to avoid unnecessary loops
									return false;
								}
							});
					}
				}); // END for each rawfiles_structure
				
				// Now RawFileConditions_copy contains only the used conditions
				// Refresh conditions list with these conditions
				
				$('#conditions_list').empty();
				
				// Retrieve them and erase possible dupplicates
				var non_un_condiitions = RawFileConditions_copy.map(function(value, index){return value.condition;});
				var unique_conditions = ArrNoDupe(non_un_condiitions);
				
				// Populate the list
				$.each(unique_conditions, function (idx, my_un_cond)
					{
						$("#conditions_list").append("<option value='" + my_un_cond + "' style='padding: 2px; font-size: 125%;' selected='true'>" + my_un_cond + "</option>");
					}); // END for each unique_conditions
					
					// If the dataset is a test dataset it is possible that specific conditions should be selected for comparison
					// run select_labels_according_to_test_dataset now to compensate this
					
					select_labels_according_to_test_dataset();
					
					// Since new conditions were introduced we should initialize rename (merging) and create mqlabels
					if(AllowMergeLabels) InitializeRename();
					create_my_all_mq_labels();
					
					
					on_exppddatachange();
					
					// Since new conditions were introduced make sure there are no invalid Label Swaps
					eraseInvalidLabelSwaps();
		}
	}
	
	
	var ons2btnfclick = function() {
		
		// Bound to forward clock of the second step:
		// If the user successfuly entered second step information, log the following data to a php log file
		var provenprocprogram = "";
		var typeofdataset = "";
		if (document.getElementsByName("exppddata")[0].checked)
		{
			provenprocprogram = "PD";
		}
		else
		{
			provenprocprogram = "MQ";
		}
		if (isLabelFree)
		{
			typeofdataset = "LabelFree";
		}
		else
		{
			if (isIsobaricLabel)
			{
				typeofdataset = "IsobaricLabel";
			}
			else
			{
				typeofdataset = "PrecursorIon";
			}
		}
		var toappend = "ProcProgram: " + provenprocprogram + "\n" + "ExperimentType: " + typeofdataset + "\n" + "ReachedStep2: T\n" + "Test_dataset: " + (log_test_dataset ? "T" : "F") + "\n";
		log_test_dataset = false;
		log_php(toappend);
		$("#s2btnupld").prop('disabled', false); // When leaving step 1 make sure the upload button comes back enabled (it may have been disabled if the user chose a test dataset)
	}
	
	
	var ons2LFQConditionsOK_click = function() {
		
		// Bound to Ok button in assign new condition dialog box
		// ons2LFQConditionsOK_click
		//
		// If the new text box is filled with a new condition it assigns the conditions to the selected raw files in
		// raw files table and appends this condition to the conditions drop down list in assign conditions
		// dialog box
		var CondToAdd;
		if ($("#s2LFQConditionsNew").prop("value") != "")
		{
			if (typeof($("#s2LFQConditionsNew").prop("value")) == "undefined") return;
			
			// Add new condition to LFQconditions
			LFQconditions.push($("#s2LFQConditionsNew").prop("value"));
			CondToAdd = $("#s2LFQConditionsNew").prop("value");
			
			// Erase new-condition text box's content
			$("#s2LFQConditionsNew").val("");
			
			// refresh the drop down list
			var ListLength = $("#s2LFQConditionsNew").prop("length");
			var optiontext = '<option val="' + ListLength + '">' + CondToAdd + '</option>';
			$("#s2LFQConditionsList").append(optiontext);
			$("#s2LFQConditionsList").attr("disabled", false);
		}
		else
		{
			// If there is a selected condition in the dropdown list simply assign this to the selected raw files
			if ($("#s2LFQConditionsList").prop("value") != "")
			{
				CondToAdd = $("#s2LFQConditionsList").prop("value");
			}
		}
		
		// Here CondToAdd stores the condition to be assigned to the selected raw files
		if (typeof(CondToAdd) == "undefined") return;
		
		// First clean up all RawFileConditions rows that will be overwritten
		// all RawFileConditions entries that correspond to a selected item in RawFiles table will be
		// first deleted to be reset in RawFileConditions
		var items = $('#rawfiles_tbl_allfiles').find('.rawfiles_tbl_td_selected');
		var item_idxs_to_splice = [];
		$.each(RawFileConditions, function (idx, my_cond)
		{
			
			for (var i = 0; i < items.length; i++) {
				var items_tds = $(items[i]).find('td');
				var items_name = items_tds[0];
				var name_txt = $(items_name).text();
				if(my_cond.name == name_txt)
				{
					item_idxs_to_splice.push(idx);
				}
			}
		}); // END for each RawFileConditions
		
		// Itterate the array from last element to first since we delete items in parallel
		// and if we iterate the array we will skip some elements
		for (var i = item_idxs_to_splice.length -1; i >= 0; i--)
		{
			RawFileConditions.splice(item_idxs_to_splice[i], 1);
		}
		
		// Now update the RawFileConditions list and the "condition" cell
		// of the corresponding rows in rawfiles table
		for (var i = 0; i < items.length; i++)
		{
			var items_tds = $(items[i]).find('td');
			var items_cond = items_tds[4];
			var items_name = items_tds[0];
			RawFileConditions.push({name: $(items_name).text(), condition: CondToAdd});
			$(items_cond).text(CondToAdd);
		}
		
		//Deselect all rows
		var my_table = $('#rawfiles_tbl_allfiles').dataTable();
		var my_rows = my_table._('tr', {"filter":"applied"});
		$.each(my_rows, function (idx, cur_row) {
			var my_tmp = cur_row["DT_RowId"];
			$(document.getElementById(my_tmp.toString())).removeClass('rawfiles_tbl_td_selected');
		});
		
		//Also check if the next button should be enabled now
		set_s22btnf();
		//Also refresh the fractions:
		refresh_fractions();
	}
	
	
	var ons4btnbclick = function() {
		
		// Bound to Step 4 back button: display a shadow on server feedback box to make it look nice
		$("#server_feedback").css("box-shadow", "")	
	}
	
	
	var ons4btnfclick = function() {
		
		// Bound to step 4 forward button, ons4btnfclick shows the Feedback button on the side
		// if it has not been already shown.
		if (feedbackAnimated == false)
		{
			$("#feedback_div").animate({
				width: "30px",
				opacity: 1
			}, 1000, function(){});
			feedbackAnimated = true;
		}
	}
	
	
	function onShowAdvOptionsClick() {
		
		// Shows the advanced options dialog to let the user choose p threshold and more
		$(".expparamsDlg").css({"left" : ($("body").width()/2) - ($("#divAdvancedOptions").width()/2)});
		$('body').append('<div id="mask"></div>');
		$("#divAdvancedOptions").fadeIn(300);
		$('#mask').fadeIn(300);
		$("#ADVoptionsOK").on("click", function () {
			//ADD OK RESULT HERE
		});
		InitializeAdvParam();
	}
	
	
	var onShowDialog = function (selector){
		
		// When a dialog box is shown:
		$(".expparamsDlg").css({"left": ($("body").width() / 2) - ($(selector).width() / 2)});
		$('body').append('<div id="mask"></div>');
		$(selector).fadeIn(300);
	}
	
	
	var onuserFeedbackchange = function() {
		
		// When the user alters the feedback text inform him of how many characters they have left to type
		// the text box has maxlength property set to 600 so they are not able to type more than that
		var mynum = 600 - parseInt($("#userFeedback").val().length);
		var mytext = mynum + " characters left";
		$("#FBcharleft").text(mytext);
	}
	
	
	var onVIDDiscard_click = function() {
		
		// When various parameters dialog box is discarded
		dlgFadeoutInfo();
	}
	
	
	var on_exppddatachange = function() {
		
		// Called when ProteoSign realizes that the dataset is from PD or MQ
		// this function serves the only reason to inform
		// the user for the procprogram (MQ or PD) by altering the quantitation_prog_lbl
		// this function will run only in label-free data since in labelled ones,
		// quantitation_prog_lbl is set in another way
		if (isLabelFree)
		{
			if ($("#s3expparams input[name='exppddata']").prop("checked") == false)
			{
				//MaxQuant:
				$("#quantitation_prog_lbl").text("\u2014 Raw data were quantified with MaxQuant");
			}
			else
			{
				//Proteome Discoverer:
				$("#quantitation_prog_lbl").text("\u2014 Raw data were quantified with Proteome Discoverer\u2122");
			}
		}
	}
	
	
	var postClientServerClientInfo = function () {
		
		// Gets ProteoSign's version to be displayed correctly
		var thedata = new FormData();
		thedata.append('session_id', sessionid);
		$.ajax({
			url: cgi_bin_path + 'get_serverclieninfo.php',
			type: 'POST',
			// Form data
			data: thedata,
			//Options to tell jQuery not to worry about content-type.
			processData: false,
			cache: false,
			contentType: false,
			beforeSend: function (jqXHR, settings) {
			}}).done(function (data, textStatus, jqXHR) {
			softversion = data.version;
			$("#scrollingtext").html("Welcome to ProteoSign ver. 2.0");
			$("#proteosignversion").html("ProteoSign version " + softversion);
			}).fail(function (jqXHR, textStatus, errorThrown) {
		});
	}
	
	
	var postTestDatasetsInfo = function () {
		
		// Get all test dataset info so that the dialog in Step 1 is ready to display
		// the test datasets when shown
		//
		// postTestDatasetsInfo runs load_test_data.php with descriptions requested true
		// and returns the descriptions of test datasets (i.e. the names that will be displayed in the dropdown list
		// in the dialog box of Step 1) and populates the dropdown list with them
		
		var thedata = new FormData();
		thedata.append('session_id', sessionid);
		thedata.append('descriptions_requested', true);
		$.ajax({
			url: cgi_bin_path + 'load_test_data.php',
			type: 'POST',
			// Form data
			data: thedata,
			//Options to tell jQuery not to worry about content-type.
			processData: false,
			cache: false,
			contentType: false
			}).done(function (data, textStatus, jqXHR) {
			//if there was a server-side error alert.
			if (!data.success) {
				msgbox("ERROR on SERVER: " + data.msg);
			}
			else if(data.queryres.desc)
			{
				var i = 1;
				$.each(data.queryres.desc, function (idx, dataset_desc)
					{
						// Populate the dropdown list
						$("#s1TestDatasetsSelection").append("<option value='" + (i++) + "'>" + dataset_desc + "</option>");
					}); // END for each data.queryres.desc
			}
			}).fail(function (jqXHR, textStatus, errorThrown) {
			msgbox("An AJAX error occurred: " + errorThrown);
		});
	}
	
	
	var questionbox = function(caption, postYes, postNo) {
		
		
		// This function shows a question box, a dialog with two options for the user, a yes and a no
		
		// First show the box
		
		$(".expparamsDlgInfo").css({"left": ($("body").width() / 2) - ($("#VariousInfoDialog").width() / 2)});
		$('body').append('<div id="maskInfo"></div>');
		$("#QuestionBoxDialog").fadeIn(300);
		$('#maskInfo').fadeIn(300);
		$("#CaptionDiv").empty();
		$("#CaptionDiv").append("<span>" + caption + "</span>");
		
		
		// In a question box, some functions should be executed after clicking on yes and no, these are sent as parameters to question box:
		// to assign them.
		
		// First remove any handlers from QBYes and QBNo and reassign the sent functions
		// also make sure that the dialog closes after hitting yes or no$("#QBYes").off("click");
		$("#QBYes").off("click");
		$("#QBYes").on("click", postYes);
		$("#QBYes").on("click", function(){
			dlgFadeoutInfo();
		});
		
		$("#QBNo").off("click");
		$("#QBNo").on("click", postNo);
		$("#QBNo").on("click", function() {
			dlgFadeoutInfo();
		});
	}
	
	
	function rawFileDataTableCMenuInit() {
		
		// Initializes the context menu that appears when right clicking the raw files table in stage 2
		
		main_context_menu = $(document).contextmenu({
			delegate: ".dataTable td",
			// Define the context menu items
			menu: [
				{title: "Select All", cmd: "slc_all"},
				{title: "Deselect All", cmd: "dslc_all"},
				{title: "Invert Selection", cmd: "inv_slc"},
				{title: "Exclude Selected", cmd: "rmv_slc"},
				{title: "Include Selected", cmd: "add_slc"},
				{title: "Clear filter", cmd: "clr_fltr"},
				{title: "Reset", cmd: "reset"},
				{title: "Select unassigned", cmd: "slc_unassigned"},
				{title: "Assign Condition", cmd: "assign_condition", visible: false},
				{title: "Replication Multiplexing", cmd: "rep_mult", visible: false}
			],
			select: function(event, ui) {
				// Bind functions to the items' click
				switch(ui.cmd){
					case "slc_all":
					
					// Select All
					var my_table = $('#rawfiles_tbl_allfiles').dataTable();
					var my_rows = my_table._('tr', {"filter":"applied"});
					$.each(my_rows, function (idx, cur_row) {
						var my_tmp = cur_row["DT_RowId"];
						$(document.getElementById(my_tmp.toString())).addClass('rawfiles_tbl_td_selected');
					});
					break;
					case "dslc_all":
					
					// Deselect All
					var my_table = $('#rawfiles_tbl_allfiles').dataTable();
					var my_rows = my_table._('tr', {"filter":"applied"});
					$.each(my_rows, function (idx, cur_row) {
						var my_tmp = cur_row["DT_RowId"];
						$(document.getElementById(my_tmp.toString())).removeClass('rawfiles_tbl_td_selected');
					});
					break;
					case "inv_slc":
					
					// Invert the selection to all
					var my_table = $('#rawfiles_tbl_allfiles').dataTable();
					var my_rows = my_table._('tr', {"filter":"applied"});
					$.each(my_rows, function (idx, cur_row) {
						var my_tmp = cur_row["DT_RowId"];
						$(document.getElementById(my_tmp.toString())).toggleClass('rawfiles_tbl_td_selected');
					});	
					break;
					case "clr_fltr":
					
					// Clear the filter
					rawfiles_tbl_allfiles_DT.search('');
					rawfiles_tbl_allfiles_DT.columns().search('');
					rawfiles_tbl_allfiles_DT.draw();
					break;
					case "reset":
					
					// Reset the table and the rawfiles structure
					$("#btnResetExpStructCoord").trigger("click");
					break;
					case "rmv_slc":
					
					// Remove the selected raw files from the analysis
					var my_table = $('#rawfiles_tbl_allfiles').dataTable();
					var items = $('#rawfiles_tbl_allfiles').find('.rawfiles_tbl_td_selected');
					for (var i = 0; i < items.length; i++)
					{
						// For each selected item
						var items_tds = $(items[i]).find('td');
						var rec_found = false;
						$.each(rawfiles_structure, function (idx, my_raw_file)
							{
								// find the corresponding entry in rawfiles_structure
								if(my_raw_file.rawfile == $(items_tds[0])[0].textContent)
								{
									// When the raw file is found set it to used: false
									my_raw_file.used = false;
									rec_found = true;
									return false;
								}
							}); // END for each rawfiles_structure
							if (!rec_found)
							{
								// If an entry was not found simply add it
								rawfiles_structure.push({rawfile: $(items_tds[0]).text(), biorep: "-", techrep: "-", fraction: "-", used : false});
							}
							// Strikeout the raw file in the table
							$(items_tds[0]).css("text-decoration", "line-through");
							$(items_tds[1]).css("text-decoration", "line-through");
							$(items_tds[2]).css("text-decoration", "line-through");
							$(items_tds[3]).css("text-decoration", "line-through");
							if (isLabelFree == true) $(items_tds[4]).css("text-decoration", "line-through");
					}
					// Check if Forward button in Stage 1 should be enabled
					set_s22btnf();
					// Deselect all - disable reset and reset fractions
					$('#rawfiles_tbl_allfiles tbody tr').removeClass('rawfiles_tbl_td_selected');
					$("#btnResetExpStructCoord").prop('disabled', rawfiles_structure.length == 0);
					refresh_fractions();
					break;
					case "add_slc":
					
					// Add selected to the analysis
					var my_table = $('#rawfiles_tbl_allfiles').dataTable();
					// Find the selected items:
					var items = $('#rawfiles_tbl_allfiles').find('.rawfiles_tbl_td_selected');
					for (var i = 0; i < items.length; i++)
					{
						// For each selected item
						var items_tds = $(items[i]).find('td');
						var rec_found = false;
						$.each(rawfiles_structure, function (idx, my_raw_file)
							{
								// Find the corresponding entry in rawfiles_structure
								if(my_raw_file.rawfile == $(items_tds[0])[0].textContent)
								{
									if(!(my_raw_file.biorep == '-' && my_raw_file.techrep == '-' && my_raw_file.fraction == '-'))
									{
										// When the raw file is found set it to used: true
										my_raw_file.used = true;
										rec_found = true;
									}
									else
									{
										// If the brep trep and frac were all set to '-' it means that the user never assigned a structure coordinate and simply removed the rawfile from the analysis but now changed his mind
										// remove the rawfile from the array totally as he/she never removed and re-added the raw file
										rawfiles_structure.splice(idx, 1);
									}
									return false;
								}
							}); // END for each rawfiles_structure
							// Remove strikeout from the corresponding table row
							$(items_tds[0]).css("text-decoration", "none");
							$(items_tds[1]).css("text-decoration", "none");
							$(items_tds[2]).css("text-decoration", "none");
							$(items_tds[3]).css("text-decoration", "none");
							if (isLabelFree == true) $(items_tds[4]).css("text-decoration", "none");
					}
					// Check if Forward button in Stage 1 should be enabled
					set_s22btnf();
					// Deselect all - disable reset and reset fractions
					$('#rawfiles_tbl_allfiles tbody tr').removeClass('rawfiles_tbl_td_selected');
					$("#btnResetExpStructCoord").prop('disabled', rawfiles_structure.length == 0);
					refresh_fractions();
					break;
					case "slc_unassigned":
					
					// Select unassigned rows
					// First get all rows where the filter applies
					var my_table = $('#rawfiles_tbl_allfiles').dataTable();
					var my_rows = my_table._('tr', {"filter":"applied"});
					$.each(my_rows, function (idx, cur_row)
						{
							// For each row
							var my_tmp = cur_row["DT_RowId"];
							var rec_found = false;
							var unused_file = false;
							// First find if the rawfile is in rawfiles_structure
							$.each(rawfiles_structure, function (idx, my_raw_file)
								{
									if(my_raw_file.rawfile == cur_row["fname"])
									{
										rec_found = true;
									}
								}); // END for each rawfiles_structure
								
								if (isLabelFree == false)
								{// In labelled data a row is unassigned if it is not found in the rawfiles_structure array
									if(rec_found == false)
									{
										$(document.getElementById(my_tmp.toString())).addClass('rawfiles_tbl_td_selected');
									}
									else
									{
										$(document.getElementById(my_tmp.toString())).removeClass('rawfiles_tbl_td_selected');
									}
								}
								else
								{// in label free the row is also unassigned if a condition is not assigned to it
									var cond_found = false;
									$.each(RawFileConditions, function (idxC, my_cond)
										{
											if (cur_row["fname"] == my_cond.name)
											{
												cond_found = true;
											}
										});
										if(rec_found == true && cond_found == true)
										{
											$(document.getElementById(my_tmp.toString())).removeClass('rawfiles_tbl_td_selected');
										}
										else
										{
											$(document.getElementById(my_tmp.toString())).addClass('rawfiles_tbl_td_selected');
										}							
								}
								
								
						});
						break;
						case "assign_condition":
						// Assign condition dialog in case of label free data
						onAssignCondition();
						break;
						case "rep_mult":
						showRepMultDialog();
						break;
				}
			},
			beforeOpen: function(event, ui) {
				var $menu = ui.menu,
				$target = ui.target,
				extraData = ui.extraData;
				ui.menu.zIndex(9999); // Bring to top
			}
		});
	}
	
	
	function rawFileDataTableInit() {
		
		// Initializes the raw data table in stage 2
		rawfiles_tbl_allfiles_DT = $('#rawfiles_tbl_allfiles').DataTable({
			paging: false,
			bInfo: false,
			select: true,
			"pageLength": 308,
			scrollY: "142px",
			scrollCollapse: true,
			
			// Define the columns
			"columnDefs": [
				{
					targets: [1, 2, 3],
					width: '4%',
					className: "dt-center"
				},
				{
					targets: 0,
					width: '70%',
					className: "dt-left",
					title: "Raw File"
				},
				{
					targets: 1,
					title: "#B"
				},
				{
					targets: 2,
					title: "#T"
				},
				{
					targets: 3,
					title: "#F"
				},
				{
					targets: 4,
					title: "Condition",
					width: '18%',
					visible: false,
					className: "dt-center"
				}
			],
			"dom": '<"top ibvspace"f>rtT',
			"sDom": 'lfrtip',
			"oLanguage": {
				"sSearch": "Filter: "
			},
			"aoColumns": [
				{"mData": "fname"},
				{"mData": "brep"},
				{"mData": "trep"},
				{"mData": "frac"},
				{"mData": "cond"}
			]
		});
	}
	
	
	var RefreshConditionsList = function(my_conditions) {
		
		// Refresh the conditions list given an array of conditions
		$("#conditions_list").empty();
		$.each(my_conditions, function (idx, lblname_i) {
			$("#conditions_list").append("<option value='" + lblname_i + "' style='padding: 2px; font-size: 125%;' selected='true'>" + lblname_i + "</option>");
		}); // END for each condition
	}
	
	
	var Refresh_conds_list = function() {
		
		// This routine refreshes the conditions list so that in case we have label merging
		// only applicable conditions will be displayed
		// if 2 or more conditions are merged to another one only the last one is available for comparisons
		// and it should be displayed in bold
		// all not merged conditions (authentic conditions) are also available for comparisons
		// and are displayed in normal text
		// If label merging (aka label renaming) is disabled simply return
		// e,g, the following RenameArray merges 3 conditions to Lung_SCC and three others
		// to Lung_ADC
		//
		// +---+----------+
		// | 0 | Lung_SCC |
		// +---+----------+
		// | 1 | Lung_SCC |
		// +---+----------+
		// | 2 | Lung_SCC |
		// +---+----------+
		// | 3 | Lung_ADC |
		// +---+----------+
		// | 4 | Lung_ADC |
		// +---+----------+
		// | 5 | Lung_ADC |
		// +---+----------+
		//
		// Only Lung_SCC and Lung_ADC should be displayed in the conitions list and they should be both in bold
		
		if (!AllowMergeLabels) return;
		
		// Clear the list
		$("#conditions_list").empty();
		
		var temp_array = [];
		$.each(RenameArray, function(idx, my_row)
			{
				// Get the second column of renamearray to temp_array
				temp_array.push(my_row[1]);
			}); // END for each RenameArray
			
			var unique_conditions = ArrNoDupe(temp_array); // Get only unique values of temp_array (e.g. Lung_ADC - Lung_SCC)
			
			// The second column always contains all conditions to be displayes since if a condition is authentic it is simply matched to itself
			// in the RenameArray
			
			$.each(unique_conditions, function (idx, my_un_cond)
				{
					// This gets all conditions to be displayed
					if($.inArray(my_un_cond, AuthenticItemsinRename) == -1)
					{
						// If the condition is authentic display in normal text	
						$("#conditions_list").append("<option value='" + my_un_cond + "' style='padding: 2px; font-size: 125%; font-weight: bold;' selected='true'>" + my_un_cond + "</option>");
					}
					else
					{
						// If the condition is not authentic display in bold
						$("#conditions_list").append("<option value='" + my_un_cond + "' style='padding: 2px; font-size: 125%;' selected='true'>" + my_un_cond + "</option>");
					}
					
				}); // END for each unique_onditions
				
				// Make sure that expquantfiltlbl containing the label used for quantitation filtering contains a valid label after merging
				// in case the user had selected one label that became invalid prompt him/her by making expquantfiltlbl's border red
				
				var my_sel_item = ""; // my selected item will temporary save the item that is currently selected in expquantfiltlbl
				if (!RenameFromTestData)
				{
					my_sel_item = $("#expquantfiltlbl option:selected").val();
				}
				
				// Refresh expquantfiltlbl
				$("#expquantfiltlbl").empty();
				$.each(unique_conditions, function (idx, my_un_cond)
					{
						$("#expquantfiltlbl").append("<option oncontextmenu='return false;' value='" + my_un_cond + "'>" + my_un_cond + "</option>");
					}); // END for each unique condition
					
					// Check if the previously selected item in expquantfiltlbl is now valid
					if (!RenameFromTestData && my_sel_item !== "")
					{
						var prev_selected_item_exists = false;
						$("#expquantfiltlbl option").each(function(idx, my_option)
							{
								if($(this).val() == my_sel_item)
								{
									prev_selected_item_exists = true;
								}
							}); // END  for each expquantfiltlbl option
							
							if(prev_selected_item_exists)
							{
								// If it is still valid simply select it
								$("#expquantfiltlbl").val(my_sel_item);
							}
							else
							{
								if($("#expquantfilt").prop("checked") == true)
								{
									// If it is not valid make expquantfiltlbl's border red
									setItemAttrColor("#expquantfiltlbl", "border", "#E60000");
								}
							}
					}
					
					//Make sure that after altering the valid labels, there are no label swaps assigned that contain invalid labels. In such a case delete the swaps and prompt the user
					eraseInvalidLabelSwaps();
	}
	
	
	var Refresh_conds_list_cmenu_items = function() {
		
		// The user can choose to merge 2 or more labels to the same condition in Stage 3 by right clicking on the
		// conditions list and using the context menu that appears
		// conditions list displays all applicable for comparisons conditions that may be either original labels
		// or merged labels. For example if an experiment had the conditions L H M and L and H are
		// merged to Lung_ca, conditions list will display M and Lung_Ca
		// Notice that all conditions that are result of merging are displayed in bold thus Lung_Ca
		// in this example should be displayed in Bold
		
		//For more information on AuthenticItemsinRename see InitializeRename() comments
		
		// First, check which items on the cmenu of conditions list should be visible Merge button:
		
		// If the user has selected original labels then the context menu should show a merge button to allow him merge
		// the labels to a condition
		
		var my_items = [];
		var show_Merge = true; // Show merge if and only if all selected options are authentic labels
		
		$("#conditions_list option:selected").each(function(idx, my_item){
			// For each selected item in conditions list, check if it exists in AuthenticItemsinRename
			// if at least one is not an authentic label do not show merge
			if ($.inArray($(this).val(), AuthenticItemsinRename) == -1)
			{
				show_Merge = false;
			}
			my_items.push(my_item);
		});
		
		if(my_items.length < 2)
		{
			show_Merge = false;
		}
		
		$("#conds_list_container").contextmenu("showEntry", "same_cond", show_Merge);
		
		//Restore should be visible only if all selected items in conditions list are m,ergerd conditions
		var my_items = [];
		var show_restore = true;
		$("#conditions_list option:selected").each(function(idx, my_item){
			// For each selected item in conditions list, check if it exists in AuthenticItemsinRename
			// if at least one is an authentic label do not show restore
			if ($.inArray($(this).val(), AuthenticItemsinRename) != -1)
			{
				show_restore = false;
			}
			my_items.push(my_item);
		});
		
		if(my_items.length == 0)
		{
			show_restore = false;
		}
		
		$("#conds_list_container").contextmenu("showEntry", "restore_conds", show_restore);
		
		// The context menu might contain either a restore or a merge button
		// if none is to be displayed simply do not show the menu on right click:
		
		if (!show_restore && !show_Merge)
		{
			$('#conds_list_container').contextmenu("option", "autoTrigger", false);
		}
		else
		{
			$('#conds_list_container').contextmenu("option", "autoTrigger", true);
		}
	}
	
	
	function refresh_fractions(){
		
		// Sets all fractions correctly
		// The algorithm is slightly different in label free and labeled cases
		if (isLabelFree == false)
		{
			// First erase all fractions from all used rawfiles
			$.each(rawfiles_structure, function (idx, my_raw_file)
				{
					if (my_raw_file.used == false) return; // Do not work on unused raw files
					my_raw_file.fraction = 0;
				}); // END for each rawfiles_structure
			$.each(rawfiles_structure, function (idx, my_raw_file)
				{
					if (my_raw_file.used == false) return; // Do not work on unused raw files
					if (my_raw_file.fraction != 0) return; // Do not set fractions for second time
					var my_brep = my_raw_file.biorep;
					var my_trep = my_raw_file.techrep;
					var my_cur_fraction = 1;
					$.each(rawfiles_structure, function (idxJ, my_raw_fileJ)
						{
							if (my_raw_fileJ.used == false) return;
							if (my_raw_fileJ.biorep == my_brep && my_raw_fileJ.techrep == my_trep)
							{
								// Add successive numbers to fractions of raw files that have the same brep and trep
								my_raw_fileJ.fraction = my_cur_fraction++;
							}
						}); // END for each rawfiles_structure
				});	// END for each rawfiles_structure
		}
		else
		{//Label free cases:
			$.each(rawfiles_structure, function (idx, my_raw_file)
				{
					if (my_raw_file.used == false) return;
					var my_brep = my_raw_file.biorep;
					var my_trep = my_raw_file.techrep;
					//In label free case the fractions are automatically set to ascending order
					//in groups that have not only biorep and techrep in common but conditions as well
					var my_assigned_condition = "";
					$.each(RawFileConditions, function (idxC, my_cond)
						{
							if (my_raw_file.rawfile == my_cond.name)
							{
								my_assigned_condition = my_cond.condition;
								return false;
							}
						}); // END for each RawFileConditions
						
					//if my_assigned_condition == "" then no condition has been assigned and this raw file can not be a part of any groups
					if (my_assigned_condition == "")
					{
						return; // continue
					}
					
					// If code reaches this line my_assigned_condition contains the assigned condition of a rawfile
					var my_cur_fraction = 1;
					$.each(rawfiles_structure, function (idxJ, my_raw_fileJ)
						{
							if (my_raw_fileJ.used == false) return; // Do not work on unused raw files
							
							//Find out the assigned condition (if any) of the J raw file
							var my_assigned_conditionJ = "";
							$.each(RawFileConditions, function (idxJC, my_condJ)
								{
									if (my_raw_fileJ.rawfile == my_condJ.name)
									{
										my_assigned_conditionJ = my_condJ.condition;
										return false;
									}
								}); // END for each RawFileConditions
								
								//if my_assigned_conditionJ == "" then no condition has been assigned and this raw file can not be a part of any groups
								if (my_assigned_conditionJ == "")
								{
									return;
								}
								// else if both [brep], [trep] and [condition assigned] are the same, these rawfiles should
								// get successive fractions
								if (my_raw_fileJ.biorep == my_brep && my_raw_fileJ.techrep == my_trep && my_assigned_condition == my_assigned_conditionJ)
								{
									my_raw_fileJ.fraction = my_cur_fraction++;
								}
						});// END for each rawfiles_structure
				});// END for each rawfiles_structure
		}
		
		//Show the new structure to the user:
		var all_items = $('#rawfiles_tbl_allfiles').find('tr');
		
		for (var i = 0; i < all_items.length; i++)
		{
			var items_tds = $(all_items[i]).find('td');
			var items_frac = items_tds[3];
			
			$.each(rawfiles_structure, function (idx, my_raw_file)
				{
					
					if (my_raw_file.rawfile == $(items_tds[0]).text())
					{
						$(items_frac).text(my_raw_file.fraction);
						return false;
					}
				});
		}
	}
	
	
	var refresh_LFQ_conditions_from_test_data = function() {
		
		// refresh_LFQ_conditions_from_test_data erases conditions from RawFileConditions
		// that are not used and refreshes the LabelFree conditions list in Assign condition dialog box
		
		RawFileConditions_copy = RawFileConditions.slice();
		$.each(rawfiles_structure, function (idx, my_raw_file)
			{
				if (my_raw_file.used == false)
				{//for each rawfile not used splice the respective row
					$.each(RawFileConditions_copy, function (idxC, my_cond)
						{
							if (my_raw_file.rawfile == my_cond.name)
							{// If there is a row in RawFileConditions_copy corresponding to a not used rawfile, erase this row
								RawFileConditions_copy.splice(idxC, 1);
								// Normally there can not be two rows in RawFileConditions_copy with the same name so escape the loop if erased a row
								// to avoid unnecessary loops
								return false;
							}
						});
				}
			}); // END for each rawfiles_structure
			
			// Refresh LabelFree conditions list in Assign condition dialog box with the unique conditions that are currently usec
			var non_un_condiitions = RawFileConditions_copy.map(function(value, index){return value.condition;});
			var unique_conditions = ArrNoDupe(non_un_condiitions);
			$.each(unique_conditions, function (idx, my_un_cond)
				{
					$("#s2LFQConditionsList").append("<option val='undefined'>" + my_un_cond + "</option>");
				});
	}
	
	
	var renderRSSData = function (data, renderelem) {
		
		// Gets the data returned by get_rss.php retrieved from an RSS feed
		// Using setTimeOut here we define a function called updateFun and
		// asynchronously execute it successively to display the next RSS element
		
		// First get the elements and shuffle them
		var prev_html = "<notset>";
		rss_i = 0;
		var items = Object.keys(data);
		items = shuffle(items);
		
		// Define the updateFun and start executing it
		var updateFun = function ()
		{
			// If the current html content of renderelem is not the same as the last one set by this function, it means some other function has set the html content
			// so we have to terminate this infinite update of renderelem's content by not calling the setTimeout function (essentially calling oneself) and just returning.
			if (prev_html.localeCompare("<notset>") != 0 && $(renderelem).html().localeCompare(prev_html) != 0)
			{
				return;
			}
			
			// If we performed a whole loop get back to the first element
			if (rss_i == (items.length - 1)) {
				rss_i = 0;
			}
			var item = items[rss_i++];
			
			if (!/^http/.test(data[item]) || /gif|png|jpg|jpeg$/.test(data[item]))
			{
				// In case of an invalid item reexecute updateFun to move to the next one
				inter = setTimeout(updateFun, 100);
				return;
			}
			
			// renderelem is the element that will display the data from RSS:
			$(renderelem).hide();
			prev_html = '<p>' + RSS_prefix + '<br><br><a href="' + data[item] + '" target="_blank"><strong>' + item + '</strong></a></p>';
			$(renderelem).html(prev_html).fadeIn(300);
			
			// Compute how much will be needed for the next phrase to be read
			var intertime = (Math.log((item.split(/\s/).length) + 1) / Math.log(1.6)) * 1000 + 1000;
			inter = setTimeout(updateFun, intertime);
		}
		
		// Begin asynchronous procedure
		var inter = setTimeout(updateFun, 100);
	}
	
	
	var resetCssStyles = function (selector) {
		// Reset styles based on source CSS (i.e. discard changes made through JS)
		$(selector).removeAttr("style");
	}
	
	
	var resetResultsInfoItms = function () {
		// Clear areas where results information appears.
		$("#server_feedback").empty();
		$("#dndres").empty();
		$("#results_p").html("Now analysing your data. Please wait for the results.<p style='font-size: 14px; margin-top: 8px; margin-bottom: 0;'>Or press <i>Back</i> to stop the analysis and change your parameters.</p>");
	}
	
	
	var resetSession = function () {
		// When resetSession is executed change_session php is server side executed and copies all session files
		// to a new session folder. the session is reset to a new session id in both server and client sidebar
		// this allows the user to alter the parameters and rerun the same experiment with different params.
		
		
		var thedata = new FormData();
		thedata.append('old_session_id', sessionid);
		sessionid = new Date().getTime();
		thedata.append('session_id', sessionid);
		$.ajax({
			url: cgi_bin_path + 'change_session.php', //Server script to fire up data analysis
			type: 'POST',
			// Form data
			data: thedata,
			//Options to tell jQuery not to worry about content-type.
			processData: false,
			cache: false,
			contentType: false
		}).done(function (data, textStatus, jqXHR) {});
	}
	
	
	var resetState = function (uploading_new_file = false) {
		//resetState resets all counters and data to their initial value as if PS starts over
		uploading_new_file = (typeof uploading_new_file !== 'undefined') ? uploading_new_file : false;
		if(uploading_new_file == false)
		{
			// Reset items that contained information from previous data input files (e.g. label information)
			$("#step2ToolsWrapper").css({"display": "inline-flex"});
			$("#step2RMuseinfo").css({"display": "none"});
			$("#step2InfoHeader").css({"display": "block"});
			$("#step2desc").css({"display": "block"});
			$("#s3expparamsDlgLabelsSelection").empty();
			$("#s3showdivLabelSwapRow").show();
			$("#s3div h2").html("Step 3 ");
			$("#s3advparams select[name='expquantfiltlbl']").empty();
			rawfiles = undefined;
			procProgram = "";
			peptideLabelsNamesFromFile = [];
			$("#s2uluploaders > table").empty();
			sessionid = new Date().getTime();
			uploaded_files = 0;
			$("#s2btnf").prop('disabled', true);
			types_of_files_uploaded = [];
			nToUpload = 0;
			AddedLabels = false;
			RawFileConditions = [];
			isLabelFree = false;
			LFQconditions = [];
			isIsobaricLabel = false;
			RenameFromTestData = false;
			LabelSwapArray = [];
			LSselected_raws = [];
			LS_counters_per_Add = [];
			my_all_mq_labels = "";
			itemsToRename = [];
			AuthenticItemsinRename = [];
			RenameArray = [];
			analysis_finished = false;
			datatestOK_clicked = false;
			my_lbls_toselect = [];
			
			RM_RMstep = 1;
			RMrawfilesdata = [];
			RMtagsdata = [];
			RMbrepsRepInRawFiles = true;
			RMtrepsRepInRawFiles = true;
			RMconditionsRepInRawFiles = true;
			RMtags = [];
			RMisused = false;
			stepRM2initialized = false;
			stepRM3initialized = false;
		}
	}
	
	
	function resetTableStage2() {
		
		// Bound to reset clicking on stage 2 and reset on stage2_table's contextmenu
		// resetTableStage2 resets the tabe in stage 2 and the relative variables
		
		reset_reps();
		$("#btnResetExpStructCoord").prop('disabled', true);
		$("#s22btnf").prop('disabled', true);
	}
	
	
	var reset_reps = function () {
		
		// reset_reps resets the table in Stage 2 and the main variables related to exp structure (repcounts rawfiles_structure and RawFileConditions)
		bioreps = 0;
		techreps = 0;
		fractions = 0;
		rawfiles_structure = [];
		RawFileConditions = [];
		rep_counts = {biorep: []};
		rawfiles_tbl_allfiles_DT.clear();
		
		// Repopulate the table
		
		$.each(rawfiles, function (idx, filename_i) {
			rawfiles_tbl_allfiles_DT.row.add(
				{
					'fname': filename_i,
					'brep': '-',
					'trep': '-',
					'frac': '-',
					'cond': '-',
					'DT_RowClass': "rawfiles_tbl_td_not_selected",
					'DT_RowId': 'tr_' + filename_i
				}
			);
		});
		
		// Redarw the table
		rawfiles_tbl_allfiles_DT.draw();
		
		// Bind a click event so that when Shift is pressed wile clicking 2 consecutive rawfiles on the table
		// all rawfiles between these two toggle their selection - somehow like in Windows explorer when selecting files
		$('#rawfiles_tbl_allfiles tbody tr').click(function (event) {
			if (event.shiftKey) {
				if (lastclicked_rawfiles_tbl_tr !== null) {
					var i1 = $('#rawfiles_tbl_allfiles tbody tr').index(lastclicked_rawfiles_tbl_tr);
					var i2 = $('#rawfiles_tbl_allfiles tbody tr').index(this);
					var trs = $('#rawfiles_tbl_allfiles tbody tr');
					if (i2 > i1) {
						for (var i = (i1 + 1); i <= i2; i++) {
							$(trs[i]).toggleClass('rawfiles_tbl_td_selected');
						}
						} else {
						for (var i = (i1 - 1); i >= i2; i--) {
							$(trs[i]).toggleClass('rawfiles_tbl_td_selected');
						}
					}
				}
				} else {
				$(this).toggleClass('rawfiles_tbl_td_selected');
			}
			lastclicked_rawfiles_tbl_tr = this;
		});
	}
	
	
	function RMCMenuInit() {
		
		// Initializes the context menu displayed on Right click on the data table in RM dialog		RMcontext_menu = $("#RMtablecontainer").contextmenu({
		// See rawFileDataTableCMenuInit comments for more information
		RMcontext_menu = $("#RMtablecontainer").contextmenu({
			delegate: ".RMrawfilestbl td",
			menu: [
				{title: "Select All", cmd: "slc_all"},
				{title: "Deselect All", cmd: "dslc_all"},
				{title: "Invert Selection", cmd: "inv_slc"},
				{title: "Exclude Selected", cmd: "rmv_slc"},
				{title: "Include Selected", cmd: "add_slc"},
				{title: "Clear filter", cmd: "clr_fltr"},
				{title: "Reset", cmd: "reset"},
				{title: "Select unassigned", cmd: "slc_unassigned"}
			],
			select: function(event, ui) {
				switch(ui.cmd){
					case "slc_all":
					var my_table = $('#RMrawfiles').dataTable();
					var my_rows = my_table._('tr', {"filter":"applied"});
					$.each(my_rows, function (idx, cur_row) {
						var my_tmp = cur_row["DT_RowId"];
						$(document.getElementById(my_tmp.toString())).addClass('rawfiles_tbl_td_selected');
						RM_set_row(my_tmp.toString().substr(5), true, 7);
					});
					break;
					case "reset":
					if (RM_RMstep == 2)
					{
						RMrawfilesdata = [];
						var mycounter = 1;
						$.each(rawfiles, function (idx, filename_i) {
							RMrawfilesdata.push([mycounter++, filename_i, '-', '-', '-', '-', true, false]);
						});
						RM_refresh_fractions();
					}
					else if (RM_RMstep == 3)
					{
						RMtagsdata = [];
						RMtags = peptideLabelsNamesFromFile;
						var mycounter = 1;
						$.each(RMtags, function (idx, tag_i) {
							RMtagsdata.push([mycounter++, tag_i, '-', '-', '-', '-', true, false]);
						});
					}
					RM_redraw_table();
					//notice the absence of break, reseting the state also deselects all
					case "dslc_all":
					RM_deselect_all();
					break;
					case "inv_slc":
					var my_table = $('#RMrawfiles').dataTable();
					var my_rows = my_table._('tr', {"filter":"applied"});
					$.each(my_rows, function (idx, cur_row) {
						var my_tmp = cur_row["DT_RowId"];
						$(document.getElementById(my_tmp.toString())).toggleClass('rawfiles_tbl_td_selected');
						RM_set_row(my_tmp.toString().substr(5));
					});	
					break;
					case "clr_fltr":
					RMrawfilesDT.search('');
					RMrawfilesDT.columns().search('');
					RMrawfilesDT.draw();
					break;
					case "rmv_slc":
					var my_table = $('#RMrawfiles').dataTable();
					var items = $('#RMrawfiles').find('.rawfiles_tbl_td_selected');
					for (var i = 0; i < items.length; i++) {
						var items_tds = $(items[i]).find('td');
						RM_set_row($(items_tds[0]).text(), 'false', 6);//set the row to unused
					}
					if (RM_RMstep == 2) RM_refresh_fractions();
					RM_redraw_table();
					RM_deselect_all();
					break;
					case "add_slc":
					var my_table = $('#RMrawfiles').dataTable();
					var items = $('#RMrawfiles').find('.rawfiles_tbl_td_selected');
					for (var i = 0; i < items.length; i++) {
						var items_tds = $(items[i]).find('td');
						RM_set_row($(items_tds[0]).text(), 'true', 6);//set the row to unused
					}
					if (RM_RMstep == 2) RM_refresh_fractions();
					RM_redraw_table();
					RM_deselect_all();
					break;
					case "slc_unassigned":
					var my_table = $('#RMrawfiles').dataTable();
					var my_rows = my_table._('tr', {"filter":"applied"});
					$.each(my_rows, function (idx, cur_row) {
						var mustbeselected = false;
						if (RM_RMstep == 3)
						{
							RMbrepsRepInRawFiles = !RMbrepsRepInRawFiles;
							RMtrepsRepInRawFiles = !RMtrepsRepInRawFiles;
							RMconditionsRepInRawFiles = !RMconditionsRepInRawFiles;
						}
						if (RMbrepsRepInRawFiles && RM_get_row(cur_row["DT_RowId"].toString().substr(5), 2) == "-") mustbeselected = true;
						if (RMtrepsRepInRawFiles && RM_get_row(cur_row["DT_RowId"].toString().substr(5), 3) == "-") mustbeselected = true;
						if (RMconditionsRepInRawFiles && RM_get_row(cur_row["DT_RowId"].toString().substr(5), 5) == "-") mustbeselected = true;
						if (RM_RMstep == 3)
						{
							RMbrepsRepInRawFiles = !RMbrepsRepInRawFiles;
							RMtrepsRepInRawFiles = !RMtrepsRepInRawFiles;
							RMconditionsRepInRawFiles = !RMconditionsRepInRawFiles;
						}
						if (RM_get_row(cur_row["DT_RowId"].toString().substr(5), 6) == 'false') mustbeselected = false;
						if (mustbeselected)
						{
							$(document.getElementById(cur_row["DT_RowId"].toString())).addClass('rawfiles_tbl_td_selected');
							RM_set_row(cur_row["DT_RowId"].toString().substr(5), true, 7);
						}
						else
						{
							$(document.getElementById(cur_row["DT_RowId"].toString())).removeClass('rawfiles_tbl_td_selected');
							RM_set_row(cur_row["DT_RowId"].toString().substr(5), false, 7);
						}
					});
					break;
				}
				RM_check_next_enable(); //check if the next button in RM Dialog must be enabled
			},
			beforeOpen: function(event, ui) {
				var $menu = ui.menu,
				$target = ui.target,
				extraData = ui.extraData;
				ui.menu.zIndex(99999);
			}
		});
	}
	
	
	function RMDataTableInit() {
		
		// Initializes the data table of raw files in Replication Multiplexing dialog
		RMrawfilesDT = $("#RMrawfiles").DataTable({
			paging: false,
			bInfo: false,
			select: true,
			"pageLength": 308,
			scrollY: "200px",
			scrollCollapse: true,
			"columnDefs": [
				{
					targets: [1, 2, 3, 4],
					width: '4%',
					className: "dt-center"
				},
				{
					targets: 0,
					width: '100%',
					className: "dt-left",
					title: "Raw File"
				},
				{
					targets: 1,
					title: "#B",
				},
				{
					targets: 2,
					title: "#T",
				},
				{
					targets: 3,
					title: "#F",
				},
				{
					targets: 4,
					title: "Condition",
					width: '18%',
					className: "dt-center"
				}
			],
			"dom": '<"top ibvspace"f>rtT',
			"sDom": 'lfrtip',
			"oLanguage": {
				"sSearch": "Filter: "
			},
			"aoColumns": [
				{"mData": "fname"},
				{"mData": "brep"},
				{"mData": "trep"},
				{"mData": "frac"},
				{"mData": "cond"}
			]
		});
	}
	
	
	var RM_check_next_enable = function() {
		
		// If all options are set correctly in the current RM step enable the next button
		// Use is_RM_ready to check the validity of user's options
		
		var next_should_be_enabled = true;
		if (RM_RMstep == 2)
		{
			next_should_be_enabled = is_RM_ready(2);
		}
		else if (RM_RMstep == 3)
		{
			next_should_be_enabled = is_RM_ready(3);
		}
		$("#RMDialogNext").attr('disabled', !next_should_be_enabled);
	}
	
	
	var RM_deselect_all = function() {
		
		// Deselects all filtered rows from the RM datatable
		
		var my_table = $('#RMrawfiles').dataTable();
		var my_rows = my_table._('tr', {"filter":"applied"});
		$.each(my_rows, function (idx, cur_row) {
			var my_tmp = cur_row["DT_RowId"];
			$(document.getElementById(my_tmp.toString())).removeClass('rawfiles_tbl_td_selected');
			RM_set_row(my_tmp.toString().substr(5), false, 7);
		}); // END for each filtered row
	}
	
	
	var RM_getRMconditions = function() {
		
		// This routine gets the conditions set in Replication Multiplexing as a one-dimension array
		// used to display the conditions to the conditions list
		
		if (RMconditionsRepInRawFiles)
		{
			var temp_array = [];
			// If the conditions are represented as different rawfiles get all conditions for each rawfile to temparray
			for(var i = 0; i<= RMrawfilesdata.length - 1; i++)
			{
				if (RMrawfilesdata[i][6] == 'true') temp_array.push(RMrawfilesdata[i][5]);
			}
			// Remove duplicates
			temp_array = ArrNoDupe(temp_array);
		}
		else
		{
			var temp_array = [];
			// If the conditions are represented as different tags get all conditions for each tag to temparray
			for(var i = 0; i<= RMtagsdata.length - 1; i++)
			{
				if (RMtagsdata[i][6] == 'true') temp_array.push(RMtagsdata[i][5]);
			}
			// Remove duplicates
			temp_array = ArrNoDupe(temp_array);
		}
		return temp_array;
	}
	
	
	var RM_get_row = function(name, valueidx) {
		
		// Returns the value of the row in RMrawfilesdata (or RMtagsdata if RM_RMstep = 3) of valueidx with the rawfilename "name"
		if (RM_RMstep == 2)
		{
			for (var i = 0 ; i <= RMrawfilesdata.length - 1; i++)
			{
				if (RMrawfilesdata[i][1] == name)
				{
					return RMrawfilesdata[i][valueidx];
				}
			}
		}
		else if (RM_RMstep == 3)
		{
			for (var i = 0 ; i <= RMtagsdata.length - 1; i++)
			{
				if (RMtagsdata[i][1] == name)
				{
					return RMtagsdata[i][valueidx];
				}
			}
		}
	}
	
	
	var RM_init_RMsteps_from_load_data = function() {
		
		// When RM is loaded from test data, the tables in RMsteps should be considered initialized
		stepRM2initialized = true;
		stepRM3initialized = true;
		
		// Refresh the check state of the check boxes in RM step 1 depending on RMbrepsRepInRawFiles, RMtrepsRepInRawFiles and RMconditionRepInRawFiles
		if (RMbrepsRepInRawFiles)
		{
			$('input[name=RMbreprepres][value=rawfiles]').prop('checked', 'checked');
		}
		else
		{
			$('input[name=RMbreprepres][value=tags]').prop('checked', 'checked');
		}
		
		if (RMtrepsRepInRawFiles)
		{
			$('input[name=RMtreprepres][value=rawfiles]').prop('checked', 'checked');
		}
		else
		{
			$('input[name=RMtreprepres][value=tags]').prop('checked', 'checked');
		}
		
		if (RMconditionsRepInRawFiles)
		{
			$('input[name=RMconditionrepres][value=rawfiles]').prop('checked', 'checked');
		}
		else
		{
			$('input[name=RMconditionrepres][value=tags]').prop('checked', 'checked');
		}
	}
	
	
	var RM_redraw_table = function() {
		
		// Visualize the RMrawfilesdata array back to the user
		// Get all the rows, selected and not
		var items = $('#RMrawfiles').find('.rawfiles_tbl_td_not_selected,.rawfiles_tbl_td_selected');
		
		// Unfortunately the column index of brep trep etc in a row is not stable and depends on options of step RM1
		var brepindex = 1;
		var trepindex = 2;
		var fracindex = 3;
		var condindex = 4;
		if (RM_RMstep == 3) condindex--; // if we are in RM# step the fractions column is not visible
		var lastelement;
		
		// RM_redraw_table is a modular function that works well in both RM step3 or 2
		// Here we will assume that we are in step RM2 but everything will work well if we are in RM3 as well
		// To compensate with the differences of these 2 steps simply toggle the representation variables
		// if we are in step 3 (this will be undone at the end of the function)
		
		if (RM_RMstep == 3)
		{
			RMbrepsRepInRawFiles = !RMbrepsRepInRawFiles;
			RMtrepsRepInRawFiles = !RMtrepsRepInRawFiles;
			RMconditionsRepInRawFiles = !RMconditionsRepInRawFiles;
		}
		
		// If breps are not represented as different rawfiles, treps, fracs and conds have one index less
		// in RMrawfiles table
		if (!RMbrepsRepInRawFiles)
		{
			trepindex--;
			fracindex--;
			condindex--;
		}
		
		// The same applies to treps
		if (!RMtrepsRepInRawFiles)
		{
			fracindex--;
			condindex--;
		}
		
		// Last element is the index of the column of the last attribute reporesented as different rawfiles
		lastelement = condindex;
		if (!RMconditionsRepInRawFiles)
		{
			lastelement--;
		}
		
		for (var i = 0; i < items.length; i++)
		{
			// For each RMrawfile row
			var items_tds = $(items[i]).find('td');
			
			// Get the name brep trep cond and frac of the row
			var items_name = items_tds[0];
			var items_biorep = items_tds[brepindex];
			var items_techrep = items_tds[trepindex];
			var items_frac = items_tds[fracindex]; //will not be used in step RM3
			var items_condition = items_tds[condindex];
			
			// Refresh the corresponding data in RMrawfile datatable using RM_get_row
			if (RMbrepsRepInRawFiles) $(items_biorep).text(RM_get_row($(items_name).text(), 2));
			if (RMtrepsRepInRawFiles) $(items_techrep).text(RM_get_row($(items_name).text(), 3));
			if (RM_RMstep == 2) $(items_frac).text(RM_get_row($(items_name).text(), 4));
			if (RMconditionsRepInRawFiles) $(items_condition).text(RM_get_row($(items_name).text(), 5));
			
			// If the corresponding row is not selected strike it out - otherwise display it as normal text
			if (!(RM_get_row($(items_name).text(), 6) == 'true'))//if not used
			{
				for (var j = 0; j <= lastelement; j++)
				{
					$(items_tds[j]).css("text-decoration", "line-through");
				}
			}
			else
			{
				for (var j = 0; j <= lastelement; j++)
				{
					$(items_tds[j]).css("text-decoration", "none");
				}
			}
		}
		
		// Undo the representation variables toggle in RM_RMstep 3
		if (RM_RMstep == 3)
		{
			RMbrepsRepInRawFiles = !RMbrepsRepInRawFiles;
			RMtrepsRepInRawFiles = !RMtrepsRepInRawFiles;
			RMconditionsRepInRawFiles = !RMconditionsRepInRawFiles;
		}
	}
	
	
	var RM_refresh_fractions = function() {
		
		// This function autocompletes tha fractions of RMrawfiles
		$.each(RMrawfilesdata, function (idx, my_row)
			{
				if (my_row[6] == false) // If the rawfile is not used
				{
					return;
				}
				var my_brep = my_row[2];
				var my_trep = my_row[3];
				var my_assigned_condition = my_row[5];
				var my_cur_fraction = 1;
				if (my_brep == "-" && RMbrepsRepInRawFiles) return;
				if (my_trep == "-" && RMtrepsRepInRawFiles) return;
				if (my_assigned_condition == "-" && RMconditionsRepInRawFiles) return;
				$.each(RMrawfilesdata, function (idxJ, my_rowJ)
					{
						// Set the same fraction to every rawfile that has the same attributes that are represented as different rawfiles
						if (my_rowJ[6] == false) return;
						if ((my_rowJ[2] == my_brep || !RMbrepsRepInRawFiles) && (my_rowJ[3] == my_trep || !RMtrepsRepInRawFiles) && (my_rowJ[5] == my_assigned_condition || !RMconditionsRepInRawFiles))
						{
							my_rowJ[4] = my_cur_fraction++;
						}
					});
			});
	}
	
	
	var RM_set_row = function(name, state, valueidx) {
		
		// Taking the name of the raw file as argument, this function finds the RMrawfilesdata record and changes the value with the index valueidx
		// Structure of RMrawfilesdata: id, name, brep, trep, frac, cond, used, selected
		// e.g. to change the state of "selected" to true of raw file "WTb12" write RM_set_row("WTb12", true, 7)
		// If state is undefined then we toogle the selected value
		// valueidx is default to 7 (selected)
		//
		// RM_set_row alters the data in RMtagsdata instead of RMrawfilesdata if we are on step RM_RMstep 3
		if (typeof(valueidx) === 'undefined') valueidx = 7;
		if (RM_RMstep == 2)
		{
			for (var i = 0 ; i <= RMrawfilesdata.length - 1; i++)
			{
				// Find the corresponding entry in RMrawfilesdata nand change the variable value
				if (RMrawfilesdata[i][1] == name)
				{
					if (typeof(state) === 'undefined')
					{
						RMrawfilesdata[i][valueidx] = !RMrawfilesdata[i][valueidx];
					}
					else
					{
						RMrawfilesdata[i][valueidx] = state;
					}
				}
			}
		}
		else if (RM_RMstep == 3)
		{
			for (var i = 0 ; i <= RMtagsdata.length - 1; i++)
			{
				// Find the corresponding entry in RMtagsdata nand change the variable value
				if (RMtagsdata[i][1] == name)
				{
					if (typeof(state) === 'undefined')
					{
						RMtagsdata[i][valueidx] = !RMtagsdata[i][valueidx];
					}
					else
					{
						RMtagsdata[i][valueidx] = state;
					}
				}
			}
		}
	}
	
	
	var rollbackStage = function (stageIndex) {
		// Executed when a "Back" type of button is clicked, rollbackStage is a togglefun for toggleNextClass
		// If the back button is pressed while in stage 1 PS resets itself like starting from the beggining
		// If stageindex = 3 (that is the analysis is currently running) PS resets the session i.e. the current
		// session files are copied to another work folder starting a new session waiting for new parameters from the user
		
		var ret = true;
		stageIndex = getItems("button.main", /s[0-9]+btnb/).length - stageIndex - 1; //hotfix - rollbackStage is called from toggleNextClass that counts stages from last to first
		// stageIndex here shows the target stage (see executeStage for more info)
		switch (stageIndex) {
			case 1:
			resetState();
			break;
			case 3:
			resetSession();
			break;
			default:
		}   
		return ret;
	}
	
	
	var SaveParams = function(output_file) {
		
		// Save the parameters of the current experiment
		// First create the contents of the text file and store them as a string
		// SaveParams takes output_file as an argument. If output_file is 0 show a messagebox
		// with a download link, otherwise simply save the parameters to the output folder on the server
		//
		// For more information on Parameters file see the comments on CreateParamsFile
		var mytext = CreateParamsFile();
		
		if (mytext !== "")
		{
			// Send the text to downloadparam file.php. downloadparam file is responsible to create a
			// test file on teh server and provide a download link
			
			var thedata = new FormData();
			thedata.append('expname', $('input[name="expid"]').val());
			thedata.append('texttodownload', mytext);
			thedata.append('session_id', sessionid);
			thedata.append('output', output_file);
			$.ajax({
				url: cgi_bin_path + 'download_param_file.php', //Server script to fire up data analysis
				type: 'POST',
				// Form data
				data: thedata,
				//Options to tell jQuery not to worry about content-type.
				processData: false,
				cache: false,
				contentType: false,
				}).done(function (data, textStatus, jqXHR) {
				if (output_file == 0)
				{
					msgbox('Download your parameters using this link: <a href="' + data.results_url + '" download target="_blank">' + $('input[name="expid"]').val() + ' parameters.txt</a>.')
				}
			});
		}
	}
	
	
	var select_labels_according_to_test_dataset = function() {
		
		// This function selects which conditions will be compared based on the test dataset's configuration
		// For each dataset, information about which conditions should be compared are stored in alocal
		// SQLite3 database server-side
		// While the test dataset is getting loaded, this info is stored im my_lbls_toselect:
		
		
		if (my_lbls_toselect.length > 0 )
		{
			$("#conditions_list option").each(function(idx, my_opt)
				{
					if($.inArray($(my_opt).val(), my_lbls_toselect) != -1)
					{
						$(my_opt).prop("selected", true);
					}
					else
					{
						$(my_opt).prop("selected", false);
					}
				});
		}
	}
	
	
	var setItemAttrColor = function (selector, attr, hexcolor) {
		
		// Set the color of an atribute (e.g. border) of an item given its id selector
		var cssstr = $(selector).css(attr);
		if (cssstr.length == 0) {
			cssstr = $(selector).css(attr + "-top-color");
			$(selector).css(attr + "-top-color", cssstr.replace(/rgb\([0-9]+, [0-9]+, [0-9]+\)$/g, hexcolor));
			cssstr = $(selector).css(attr + "-right-color");
			$(selector).css(attr + "-right-color", cssstr.replace(/rgb\([0-9]+, [0-9]+, [0-9]+\)$/g, hexcolor));
			cssstr = $(selector).css(attr + "-bottom-color");
			$(selector).css(attr + "-bottom-color", cssstr.replace(/rgb\([0-9]+, [0-9]+, [0-9]+\)$/g, hexcolor));
			cssstr = $(selector).css(attr + "-left-color");
			$(selector).css(attr + "-left-color", cssstr.replace(/rgb\([0-9]+, [0-9]+, [0-9]+\)$/g, hexcolor));
		}
		else
		{
			$(selector).css(attr, cssstr.replace(/rgb\([0-9]+, [0-9]+, [0-9]+\)$/g, hexcolor));
		}
	}
	
	
	function set_reps() {
		
		// set_reps defines the contents of rep_counts, an array that stores how many treps belong to an already defined
		// brep and how many fracs belong to an already defined trep
		// e.g. if the user is in the process of defining the experimental structure in Stage 2 and
		// has at the moment defined 6 rawfiles with the following simple structure:
		//
		// brep		trep	frac
		// 1		1		1
		// 1		2		1
		// 1		3		1
		// 1		4		1
		// 1		5		1
		// 1		6		1
		// 
		// rep_counts contain one Property (biorep) that has inside one element. This element contains an array called techrep
		// that has 6 elements (as many as the techreps). Each one of them contains an array of one element since there is one frac per techrep
		//
		// rep_counts are very important for auto-completion since in case the user wants now to assign
		// brep as 1 to two new rawfiles PS will automatically assign trep 7 and 8 to these files
		// rearranging the table in stage 2 to
		//
		// brep		trep	frac
		// 1		1		1
		// 1		2		1
		// 1		3		1
		// 1		4		1
		// 1		5		1
		// 1		6		1
		// 1		7		1
		// 1		8		1
		//
		// An example of rep_counts structure after defining the structure of an experiment with 2 breps, 2 treps each, 4 fractions each is this:
		//
		/*
			biorep
			1:
				techrep
				1:
					fraction
					1:
						1
					2:
						1
					3:
						1
					4:
						1
				2:
					fraction
					1:
						1
					2:
						1
					3:
						1
					4:
						1
			2:
				techrep
				1:
					fraction
					1:
						1
					2:
						1
					3:
						1
					4:
						1
				2:
					fraction
					1:
						1
					2:
						1
					3:
						1
					4:
						1
		*/
		//

		for (var i = 0; i < rawfiles_structure.length; i++)
		{
			var rep = "biorep";
			var max_level = 0;
			if (rawfiles_structure[i][rep] != 0)
			{
				if (!(rawfiles_structure[i][rep] in rep_counts[rep]))
				{
					rep_counts[rep][rawfiles_structure[i][rep]] = {techrep: []};
					bioreps++;
				}
			}
			rep = "techrep";
			if (rawfiles_structure[i][rep] != 0)
			{
				if (!(rawfiles_structure[i][rep] in rep_counts["biorep"][rawfiles_structure[i]["biorep"]][rep]))
				{
					rep_counts["biorep"][rawfiles_structure[i]["biorep"]][rep][rawfiles_structure[i][rep]] = {fraction: []};
					techreps++;
				}
				max_level = 1;
			}
			rep = "fraction";
			if (rawfiles_structure[i][rep] != 0)
			{
				if (!(rawfiles_structure[i][rep] in rep_counts["biorep"][rawfiles_structure[i]["biorep"]]["techrep"][rawfiles_structure[i]["techrep"]]["fraction"]))
				{
					rep_counts["biorep"][rawfiles_structure[i]["biorep"]]["techrep"][rawfiles_structure[i]["techrep"]]["fraction"][rawfiles_structure[i][rep]] = 1;
					fractions++;
				}
				max_level = 2;
			}
		}
		
		// Check if every element in Stage 2 is fully confighured so that we may advance to stage3:
		set_s22btnf();
		
		// Set reset to disabled if rawfiles_structure is empty
		$("#btnResetExpStructCoord").prop('disabled', rawfiles_structure.length == 0);
	}
	
	
	var set_RMisused = function(isused, showmessage) {
		
		// This function makes the appropriate changes in the interface in case Replication Multiplexing is decided to be used
		if (typeof(showmessage) === 'undefined') showmessage = true;
		if (!isIsobaricLabel) return;
		if (isused)
		{
			// If we decided to use Rep Mult hide the main data table and display an interactive text in its place
			$("#step2ToolsWrapper").css({"display": "none"});
			$("#step2RMuseinfo").css({"display": "inline-flex"});
			$("#step2InfoHeader").css({"display": "none"});
			$("#step2desc").css({"display": "none"});
			$("#s22btnf").prop('disabled', false);
			// Disable LS and Merge - Rename
			tempAllowLS = AllowLS
			tempAllowMergeLabels = AllowMergeLabels
			AllowLS = false;
			AllowMergeLabels = false;
			// Refresh the conditions
			RefreshConditionsList(RM_getRMconditions());
			$("#s3QuantitationFiltering").hide();
			$("#expquantfilt").prop("checked", false);
		}
		else
		{
			// If RM is not used show again the main datatable
			$("#step2ToolsWrapper").css({"display": "inline-flex"});
			$("#step2RMuseinfo").css({"display": "none"});
			$("#step2InfoHeader").css({"display": "block"});
			$("#step2desc").css({"display": "block"});
			$("#s3QuantitationFiltering").show();
			$("#expquantfilt").prop("checked", false);
			
			// Reset the variables that are related to RM
			RM_RMstep = 1;
			RMrawfilesdata = [];
			RMtagsdata = [];
			RMbrepsRepInRawFiles = true;
			RMtrepsRepInRawFiles = true;
			RMconditionsRepInRawFiles = true;
			RMtags = [];
			RMisused = false;
			stepRM2initialized = false;
			stepRM3initialized = false;
			rawfiles_tbl_allfiles_DT.columns.adjust().draw();
			if (showmessage) msgbox("Replication Multiplexing was discarded, please define the experimental coordinates of your raw files");
			set_s22btnf();
			AllowLS = tempAllowLS;
			AllowMergeLabels = tempAllowMergeLabels;
			
			// Refresh the conditions
			RefreshConditionsList(peptideLabelsNamesFromFile);
		}
		
		// Change the UI depending on whether LS or merging is allowed
		if(!AllowLS)
		{
			$("#s3showdivLabelSwapRow").hide();
		}
		else
		{
			$("#s3showdivLabelSwapRow").show();
		}
		$('#conds_list_container').contextmenu("option", "autoTrigger", AllowMergeLabels);
		
		// Set RMisused
		RMisused = isused;
	}
	
	
	function set_s22btnf() {
		
		// set_s22btnf enables the Next Button in step 2 in case every raw file is completely configured
		// and if some basic validations checks are passed
		
		if (rawfiles_structure.length == rawfiles.length)
		{	//if all files have been assigned to something
			
			// at least one rawfile should be used for the analysis
			var found_used_rec = false;
			$.each(rawfiles_structure, function (idx, my_raw_file)
				{
					if(my_raw_file.used == true)
					{
						found_used_rec = true;
						return;
					}
				});
				
				// Populate the variable below with all bioreps that are used in the analysis
				var bioreps_used = [];
				$.each(rawfiles_structure, function (idx, my_raw_file)
					{
						if(my_raw_file.used == true)
						{
							if (bioreps_used.indexOf(my_raw_file.biorep) == -1)
							{
								bioreps_used.push(my_raw_file.biorep);
							}
						}
					});
					
					
					if(isLabelFree == false)
					{
						// If the experiment is not LabelFree the only prerequisites are that at least one raw file should be selected
						// all rawfiles are configured and that at least 2 bioreps are used
						$("#s22btnf").prop('disabled', !(rawfiles.length > 0 && rawfiles_structure.length == rawfiles.length && found_used_rec && bioreps_used.length > 1));
					}
					else
					{
						//In case of label free data The next button should be enabled only if al used raw files have been assigned to a condition
						var all_files_have_label_assigned = true;
						// Also measure how many Label Free labels we have (they must be at least 1)
						var all_LFQ_labels_found = [];
						$.each(rawfiles_structure, function (idx, my_raw_file)
							{
								
								if(my_raw_file.used == true)
								{
									//foreach used file in rawfiles structure see if there is a corresponding one in
									//Rawfileconditions then all rawfiles are assigned to a condition
									var found_corresponding_file = false;
									$.each(RawFileConditions, function (idx, my_cond)
										{
											if (my_raw_file.rawfile == my_cond.name)
											{
												all_LFQ_labels_found.push(my_cond.condition);
												found_corresponding_file = true;
												return;
											}
										});
										
										if (found_corresponding_file == false)
										{
											all_files_have_label_assigned = false;
											return;
										}
								}
							});
							all_LFQ_labels_found = ArrNoDupe(all_LFQ_labels_found);
							//Here if all used files have a corresponding label (condition) assigned all_files_have_label_assigned is set to true otherwise to false
							$("#s22btnf").prop('disabled', !(rawfiles.length > 0 && rawfiles_structure.length == rawfiles.length && all_files_have_label_assigned == true && found_used_rec && bioreps_used.length > 1 && all_LFQ_labels_found.length > 1));
					}
					
					
		}
		else
		{
			// If there are rawfiles taht are not configured disable Next button in Stage 2
			$("#s22btnf").prop('disabled', true);
		}
	}
	
	
	function showHideAdvancedParameters() {
		
		//Toggle visibility of HTML division with advanced parameters
		
		var txt = $(this).text();
		if ($("#s3advparams").hasClass("hidden"))
		{
			$("#s3advparams").removeClass("hidden");
			$(this).text(txt.replace("Show", "Hide"));
			var parentdiv = $(this).closest("div");
			parentdiv.scrollTop(parentdiv.prop("scrollHeight"));
		}
		else
		{
			$("#s3advparams").addClass("hidden");
			$(this).text(txt.replace("Hide", "Show"));
		}
	}
	
	
	function showLabelSwapDialogBox() {
		
		// Show the LS Dialog Box and initialize it
		
		if(!AllowLS) return;
		$(".expparamsDlg").css({"left" : ($("body").width()/2) - ($("#divLabelSwap").width()/2)});
		$('body').append('<div id="mask"></div>');
		$("#divLabelSwap").fadeIn(300);
		$('#mask').fadeIn(300);
		InitializeLS();
	}
	
	
	var showRepMultDialog = function() {
		
		// Show the replication multiplexing dialog, set RM step to 1 and run showstepRMDialog to display
		// the elements that correspond to the first step of RM
		// for more information concerning RM please refer to our documentation
		
		RM_RMstep = 1;
		showstepRMDialog();
		
		// Show the box:
		$(".RepMultDlg").css({"left": ($("body").width() / 2) - ($("#ReplicationMultiplexingDialog").width() / 2)});
		$(".RepMultDlg").css({"top": ($("body").height() / 2) - ($("#ReplicationMultiplexingDialog").height() / 2)});
		$('body').append('<div id="maskRepMult"></div>');
		$("#ReplicationMultiplexingDialog").fadeIn(300);
		$('#maskRepMult').fadeIn(300);
	}
	
	
	var showstepRMDialog = function() {
		
		// For more information on RM please refer to our documentation and the comments on
		// onRMDialogNext_click
		
		// showstepRMDialog show the RM step designated by RM_RMstep in the rm dialog box
		// Notice that although RM steps are 3 there are simply 2 divs (#RM1 and #RM2)
		// see the comment on onRMDialogNext_click to understand why this happens
		
		switch(RM_RMstep)
		{
			case 1:
			// Show the first div and disable the back button
			$("#RM1").removeClass("hidden");
			$("#RM2").addClass("hidden");
			$("#RMDialogBack").attr("disabled", true);
			$("#RMDialogNext").attr("disabled", false);
			break;
			case 2:
			// When advancing from step RM1 to RM2 check that at least one experiment attribute is represented as tags
			// if not prompt the user and exit the function
			if ($('input[name=RMbreprepres]:checked').val() == "rawfiles" && $('input[name=RMtreprepres]:checked').val() == "rawfiles" && $('input[name=RMconditionrepres]:checked').val() == "rawfiles")
			{
				msgbox("Something must be represented in your experiment's tags");
				RM_RMstep = 1;
				return;
			}
			
			// Disable next button in RM2 and empty all input text boxes
			$("#RMDialogBack").attr("disabled", false);
			$("#RMDialogNext").attr("disabled", true);
			$("#stepRM2txtbrep").val("");
			$("#stepRM2txttrep").val("");
			$("#stepRM2txtcondition").val("");
			
			// The following three variables show how breps, treps and conditoions are represented in RM
			RMbrepsRepInRawFiles = $('input[name=RMbreprepres]:checked').val() == "rawfiles";
			RMtrepsRepInRawFiles = $('input[name=RMtreprepres]:checked').val() == "rawfiles";
			RMconditionsRepInRawFiles = $('input[name=RMconditionrepres]:checked').val() == "rawfiles";
			
			// Make the subtitle of RM2 more informative
			var infotodisplay = "";
			if (!RMbrepsRepInRawFiles && !RMtrepsRepInRawFiles && !RMconditionsRepInRawFiles)
			{
				// If everything is represented as tags (something quite unusual...) simply ask the user to advance to the next level
				infotodisplay = "<strong>ProteoSign</strong> considered all your raw files as fractions of the same MS/MS run, please click the Next button.";
			}
			else
			{
				// Give the necessary information to the user
				infotodisplay = "Please select one or more files, define ";
				if(RMbrepsRepInRawFiles) infotodisplay += "their <i>biological replicates</i>, ";
				if(RMtrepsRepInRawFiles) infotodisplay += 'their <i>technical replicates</i> (if you do not have technical replication assign "1" to all), ';
				if(RMconditionsRepInRawFiles) infotodisplay += "their <i>conditions</i>, ";
				infotodisplay += "using the text boxes below and click <strong>Assign</strong>. The <i>fractions</i> will be assigned automatically.<br><strong>Right click</strong> on the table for extended functionality.";
			}
			// Hide the tools on the bottom of the screen if nothing is represented as rawfiles
			if (!RMbrepsRepInRawFiles && !RMtrepsRepInRawFiles && !RMconditionsRepInRawFiles)
			{
				$("#stepRM2AssignmentTools").css({"display": "none"});
			}
			else
			{
				$("#stepRM2AssignmentTools").css({"display": "block"});
			}
			$("#stepRM2Info").empty();
			$("#stepRM2Info").append(infotodisplay);
			
			// Conditions mught be represented as different raw files, in this case the condition column should be shown
			if (!RMconditionsRepInRawFiles)
			{
				RMrawfilesDT.column("4").visible(false);
				$("#stepRM2txtcondition").css({"display": "none"});
			}
			else
			{
				RMrawfilesDT.column("4").visible(true);
				$("#stepRM2txtcondition").css({"display": "block"});
			}
			
			// Show the trep column only if treps are represented as different rawfiles
			if (!RMtrepsRepInRawFiles)
			{
				RMrawfilesDT.column("2").visible(false);
				$("#stepRM2txttrep").css({"display": "none"});
			}
			else
			{
				RMrawfilesDT.column("2").visible(true);
				$("#stepRM2txttrep").css({"display": "block"});
			}
			
			// Show the brep column only if breps are represented as different rawfiles
			if (!RMbrepsRepInRawFiles)
			{
				RMrawfilesDT.column("1").visible(false);
				$("#stepRM2txtbrep").css({"display": "none"});
			}
			else
			{
				RMrawfilesDT.column("1").visible(true);
				$("#stepRM2txtbrep").css({"display": "block"});
			}
			
			// Fractions column should always be visible
			RMrawfilesDT.column("3").visible(true);
			RMrawfilesDT.clear();
			$($("#RMrawfiles").DataTable().column(0).header()).text("Raw File");
			
			// RMrawfilesdata is a two dimension array that has 8 columns:
			// counter, rawfilename, brep, trep, frac, cond, used, selected
			// An example is:
			//
			// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
			// | counter | rawfilename                                               | brep | trep | frac | cond | used | selected |
			// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
			// | 1       | M456-E04-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
			// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
			// | 2       | M456-D10-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
			// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
			// | 3       | M456-E08-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
			// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
			// | 4       | M456-C09-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
			// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
			// | 5       | M456-D05-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
			// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
			//
			// the columns of the attributes represented as different rawfiles will be populated using the datatable
			// in RM2
			
			// If step2 is not initialized, initialize it and populate the datable's rows:
			if (!stepRM2initialized) RMrawfilesdata = [];
			var mycounter = 1;
			$.each(rawfiles, function (idx, filename_i) {
				RMrawfilesDT.row.add(
					{
						'fname': filename_i,
						'brep': '-',
						'trep': '-',
						'frac': '-',
						'cond': '-',
						'DT_RowClass': "rawfiles_tbl_td_not_selected", //Notice: all rows are given the rawfiles_tbl_td_not_selected as a class whether they are selected or not, the row is selected if it has the rawfiles_tbl_td_selected class and not if it does not have it
						'DT_RowId': 'RMtr_' + filename_i //DT_RowID is the id of the corresponding row and it is the most appendable element of the row, so in case we want to refer to the row e.g. the user has selected we will use this ID
					}
				);
				// Structure of RMrawfilesdata: id, name, brep, trep, frac, cond, used, selected
				if (!stepRM2initialized) RMrawfilesdata.push([mycounter++, filename_i, '-', '-', '-', '-', 'true', false]);
			});
			
			// Show RM2 and draw the table
			$("#RM1").addClass("hidden");
			$("#RM2").removeClass("hidden");
			RMrawfilesDT.columns.adjust().draw();
			RMrawfilesDT.draw();
			
			// Handle RMrawfiles (datatable) clicks
			// If a rawfile row was selected and now another one was clicked, toogle_selected all rows up to the second clicked
			//
			// RM_set_row takes valueidx == 7 as default (selected) and toogles boolean value if state is undefined
			// RM_set_row is smart enough to alter the values of RMrawfilesdata if we are on step RM2
			// and RMtagsdata if we are on RM3
			// see its comments for more information
			lastclicked_RMrawfiles_tbl_tr = null;
			$('#RMrawfiles tbody tr').click(function (event) {
				if (event.shiftKey)
				{
					if (lastclicked_RMrawfiles_tbl_tr !== null)
					{
						var i1 = $('#RMrawfiles tbody tr').index(lastclicked_RMrawfiles_tbl_tr);
						var i2 = $('#RMrawfiles tbody tr').index(this);
						var trs = $('#RMrawfiles tbody tr');
						if (i2 > i1)
						{
							for (var i = (i1 + 1); i <= i2; i++)
							{
								$(trs[i]).toggleClass('rawfiles_tbl_td_selected');
								// Toogle the selected value for the corresponding rows in RMrawfilesdata
								RM_set_row($(trs[i]).attr("id").substr(5));
							}
						}
						else
						{
							for (var i = (i1 - 1); i >= i2; i--)
							{
								$(trs[i]).toggleClass('rawfiles_tbl_td_selected');
								// Toogle the selected value for the corresponding rows in RMrawfilesdata
								RM_set_row($(trs[i]).attr("id").substr(5));
							}
						}
					}
				}
				else
				{
					$(this).toggleClass('rawfiles_tbl_td_selected');
					// Toogle the selected value for the corresponding rows in RMrawfilesdata
					RM_set_row($(this).attr("id").substr(5));
				}
				lastclicked_RMrawfiles_tbl_tr = this;
			}); // END on rawfiles click
			// Refresh fractions redraw the table and check if the next button should be enabled
			RM_refresh_fractions();
			RM_redraw_table();
			RM_check_next_enable();
			stepRM2initialized = true;
			break;
			case 3:
			
			// RM3 has the same controls as RM2 but has different rows in RMrawfiles table
			// Notice that even though the table still has the name RMrawfiles, it now contains the tags of the experiment to assign the attributes that remain to
			
			// First disable the next button and empty the assign text boxes
			$("#RMDialogNext").attr("disabled", true);
			$("#stepRM2txtbrep").val("");
			$("#stepRM2txttrep").val("");
			$("#stepRM2txtcondition").val("");
			
			// Make the subtitle of RM3 more informative
			var infotodisplay = "";
			infotodisplay = "Please select one or more tags, define ";
			if(!RMbrepsRepInRawFiles) infotodisplay += "their <i>biological replicates</i>, ";
			if(!RMtrepsRepInRawFiles) infotodisplay += 'their <i>technical replicates</i>, ';
			if(!RMconditionsRepInRawFiles) infotodisplay += "their <i>conditions</i>, ";
			infotodisplay += "using the text boxes below and click <strong>Assign</strong>.<br><strong>Right click</strong> on the table for extended functionality. When ready, click <strong>Next</strong> to submit your experimental structure.";
			$("#stepRM2Info").empty();
			$("#stepRM2Info").append(infotodisplay);
			
			// Make sure the assign tools in the bottom of the screen are visible
			$("#stepRM2AssignmentTools").css({"display": "block"});
			// Show the user the assign text boxes that correspond to the attributes that are represented as different tags
			if (RMconditionsRepInRawFiles)
			{
				RMrawfilesDT.column("4").visible(false); // conditions column
				$("#stepRM2txtcondition").css({"display": "none"});
			}
			else
			{
				RMrawfilesDT.column("4").visible(true); // conditions column
				$("#stepRM2txtcondition").css({"display": "block"});
			}
			if (RMtrepsRepInRawFiles)
			{
				RMrawfilesDT.column("2").visible(false); // treps column
				$("#stepRM2txttrep").css({"display": "none"});
			}
			else
			{
				RMrawfilesDT.column("2").visible(true); // treps column
				$("#stepRM2txttrep").css({"display": "block"});
			}
			if (RMbrepsRepInRawFiles)
			{
				RMrawfilesDT.column("1").visible(false); // breps column
				$("#stepRM2txtbrep").css({"display": "none"});
			}
			else
			{
				RMrawfilesDT.column("1").visible(true); // breps column
				$("#stepRM2txtbrep").css({"display": "block"});
			}
			
			// The fractions are not used here so hide their column
			RMrawfilesDT.column("3").visible(false); // fractions column
			// Change the first column name
			$($("#RMrawfiles").DataTable().column(0).header()).text("Tag");
			
			// Clear and repopulate the table
			RMrawfilesDT.clear();
			
			// Initialize the RMtagsdata if RM3 is not initialized
			// RMtagsdata has teh same structure as RMrawfilesdata (see the comments above)
			
			if (!stepRM3initialized) RMtagsdata = [];
			RMtags = peptideLabelsNamesFromFile; // get a copy of the user's dataset's tags
			var mycounter = 1;
			$.each(RMtags, function (idx, tagname_i) {
				RMrawfilesDT.row.add(
					{
						'fname': tagname_i,
						'brep': '-',
						'trep': '-',
						'frac': '-',
						'cond': '-',
						'DT_RowClass': "rawfiles_tbl_td_not_selected", //Notice: all rows are given the rawfiles_tbl_td_not_selected as a class whether they are selected or not, the row is selected if it has the rawfiles_tbl_td_selected class and not if it does not have it
						'DT_RowId': 'RMtr_' + tagname_i //DT_RowID is the id of the corresponding row and it is the most appendable element of the row, so in case we want to refer to the row e.g. the user has selected we will use this ID
					}
				);
				// Structure of RMtagsdata: id, name, brep, trep, frac, cond, used, selected
				// Initialize RMtagsdata if necessary
				if (!stepRM3initialized) RMtagsdata.push([mycounter++, tagname_i, '-', '-', '-', '-', 'true', false]);
			});
			
			// Redraw the table
			RMrawfilesDT.columns.adjust().draw();
			RMrawfilesDT.draw();
			
			// Redefine the click handler:
			lastclicked_RMrawfiles_tbl_tr = null;
			$('#RMrawfiles tbody tr').click(function (event) {
				if (event.shiftKey) {
					if (lastclicked_RMrawfiles_tbl_tr !== null)
					{
						var i1 = $('#RMrawfiles tbody tr').index(lastclicked_RMrawfiles_tbl_tr);
						var i2 = $('#RMrawfiles tbody tr').index(this);
						var trs = $('#RMrawfiles tbody tr');
						if (i2 > i1)
						{
							for (var i = (i1 + 1); i <= i2; i++)
							{
								$(trs[i]).toggleClass('rawfiles_tbl_td_selected');
								// Toogle the selected value for the corresponding rows in RMtagsdata
								RM_set_row($(trs[i]).attr("id").substr(5));
							}
						}
						else
						{
							for (var i = (i1 - 1); i >= i2; i--)
							{
								$(trs[i]).toggleClass('rawfiles_tbl_td_selected');
								// Toogle the selected value for the corresponding rows in RMtagsdata
								RM_set_row($(trs[i]).attr("id").substr(5));
							}
						}
					}
				}
				else
				{
					$(this).toggleClass('rawfiles_tbl_td_selected');
					// Toogle the selected value for the corresponding rows in RMtagsdata
					RM_set_row($(this).attr("id").substr(5));
				}
				lastclicked_RMrawfiles_tbl_tr = this;
			});
			// Notice that now RMrawfilesDT contains the tags of the dataset
			// Redraw the table and check if the next button should be enabled
		    RM_redraw_table();
		    RM_check_next_enable();
			stepRM3initialized = true;
			break;
			case 4:
			// If the next button set RM_RMstep to 4 this means that we just finished step 3 so RM is completed
			// Close the dialog and inform PS that RM mode is on
			dlgFadeoutRepMult();
			set_RMisused(true);
		}
		// Always resize the dialog box when changing RMsteps
		$(".RepMultDlg").css({"top": ($("body").height() / 2) - ($("#ReplicationMultiplexingDialog").height() / 2)});
	}
	
	
	function shuffle(o) {
		
		// Used to shuffle elements of an array
		//http://dzone.com/snippets/array-shuffle-javascript
		
		for (var j, x, i = o.length; i; j = Math.floor(Math.random() * i), x = o[--i], o[i] = o[j], o[j] = x)
		;
		return o;
		
	}
	
	
	function starHover(starNumber) {
		
		//bound to rating stars mousenter event
		star_hovered = starNumber;
		display_star_score();
	}
	
	
	var stepRM2Assign_click = function() {
		
		// Bound to the assign button in step RM2 
		// it sets the correct values in RMrawfilesdata (or RMtagsdata) and refreshes the table
		//
		// RMrawfilesdata (or RMtagsdata) is a two dimension array that has 8 columns:
		// counter, rawfilename, brep, trep, frac, cond, used, selected
		//
		// An example is:
		//
		// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
		// | counter | rawfilename                                               | brep | trep | frac | cond | used | selected |
		// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
		// | 1       | M456-E04-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
		// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
		// | 2       | M456-D10-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
		// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
		// | 3       | M456-E08-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
		// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
		// | 4       | M456-C09-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
		// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
		// | 5       | M456-D05-O207-T87S1-T87S2-T87S3-T88S1-T88S2-T88S3-P5329-1 | -    | -    | -    | -    | true | false    |
		// +---------+-----------------------------------------------------------+------+------+------+------+------+----------+
		//
		var items = $('#RMrawfiles').find('.rawfiles_tbl_td_selected'); //get all the selected rows
		for (var i = 0; i < items.length; i++)
		{
			var items_tds = $(items[i]).find('td');
			var items_name = items_tds[0]; // get the name of the selected row
			// Set brep, trep, condition their new values - RM_set_row is smart enough to set RMrawfilesdata or RMtagsdata
			// according to the current RMStep
			if ($("#stepRM2txtbrep").val() != "") RM_set_row($(items_name).text(), $("#stepRM2txtbrep").val(), 2);
			if ($("#stepRM2txttrep").val() != "") RM_set_row($(items_name).text(), $("#stepRM2txttrep").val(), 3);
			if ($("#stepRM2txtcondition").val() != "") RM_set_row($(items_name).text(), $("#stepRM2txtcondition").val(), 5);
			RM_set_row($(items_name).text(), 'true', 6); // Set used to true
		}
		if (RM_RMstep == 2) RM_refresh_fractions(); // Refresh fractions if necessary
		
		// Redraw the table
		RM_redraw_table();
		
		// Deselect all the rows and set selected to false in RMrawfilesdata or RMtagsdata
		$('#RMrawfiles tbody tr').removeClass('rawfiles_tbl_td_selected');
		if (RM_RMstep == 2)
		{
			for (var i = 0 ; i <= RMrawfilesdata.length - 1; i++)
			{
				RMrawfilesdata[i][7] = false;
			}
		}
		else if (RM_RMstep == 3)
		{
			for (var i = 0 ; i <= RMtagsdata.length - 1; i++)
			{
				RMtagsdata[i][7] = false;
			}
		}
		
		// Clear the input text boxes
		$("#stepRM2txtbrep").val("");
		$("#stepRM2txttrep").val("");
		$("#stepRM2txtcondition").val("");
		RM_check_next_enable();
	}
	
	
	function StringisNumber(n) {
		
		// from http://stackoverflow.com/questions/1303646/check-whether-variable-is-number-or-string-in-javascript
		// Returns true if a string is number
		return /^-?[\d.]+(?:e-?\d+)?$/.test(n);
	} 
	
	
	var string_to_array = function(my_string, curdimension) {
		
		// string_to_array folds a one-line string to an array using the assunptions described in array_to_string
		//
		// see the comments of array_to_string for more information
		// Notice that the string should not contain the character ~ and should contain the character |
		// only to denote dimension breaks
		
		var curdelimiter = "~";
		var olddelimiter = "|";
		// First reverse the dimensionality of my_string to fold the array easily:
		//
		// e.g. transform "OT2_Terhune_2012-10-09_DMC-HLNF01_01||1||1||1|OT2_Terhune_2012-10-09_DMC-HLNF01_02||1||2||1"
		// to "OT2_Terhune_2012-10-09_DMC-HLNF01_01~1~1~1~~OT2_Terhune_2012-10-09_DMC-HLNF01_02~1~2~1"
		//
		// Notice that string_to_array should be called from other functions without defining the curdimension var
		
		// The following block will be executed during the first recursion cycle only
		if (typeof(curdimension) === 'undefined')
		{
			curdimension = 1;
			var mycounter = 1;
			// Start creating repetitions of ||...| until you find the lowest dimension
			for (var i = 15; i > 0; i--)
			{
				var repeated_delimiter = olddelimiter.repeat(i);
				if (my_string.includes(repeated_delimiter))
				{
					// Replace the lowest dimension delimiter (e.g. "||") to the new delimiter (e.g. "~")
					var temp_regex = new RegExp(("\\" + olddelimiter).repeat(i), "g");
					my_string = my_string.replace(temp_regex, curdelimiter.repeat(mycounter++));
				}
			}
			mycounter--;
			curdimension = mycounter;
		}
		
		
		var retarray = [];
		// Split the string to a temporary array
		var temp_array = my_string.split(curdelimiter.repeat(curdimension));
		
		// For each elelent of temp_array recursively call string_to_array to fold even more the element
		for (var i = 0 ; i <= temp_array.length - 1; i++)
		{
			if (curdimension != 1 && temp_array[i].includes(curdelimiter.repeat(curdimension - 1)))
			{
				// If the element contains the delimiter of a higher dimension fold the element to an array
				retarray.push(string_to_array(temp_array[i], curdimension - 1));
			}
			else
			{
				// Otherwise simply add a single element value to the array
				retarray.push(temp_array[i]);
			}
		}
		return retarray;
	}
	
	
	var SwitchLFQList = function(e) {
		
		// Bound to New condition text box (Label Free dialog box - Stage 2) when a key is pressed (onkeyup)
		// it disables the list below if the text box is not empty
		var srcelem = e.target || e.srcElement;
		$("#s2LFQConditionsList").attr("disabled", srcelem.value != "");
	}
	
	
	var toggleCurrentSectionSpinner = function() {
		
		// toggleCurrentSectionSpinner toggles a turning gear in heading 2 of the current stage (section)
		
		if(sectionSpinnerOn)
		{
			// if the gear is displayed replace heading 2 with the html code the heading had before
			// displaying the gear (stored temporarily in global var sectionSpinnerText)
			$('.main_div .main_section:not(.hidden) h2').html(sectionSpinnerText);
		}
		else
		{
			// if the gear is not currently shown, temporarily save heading 2 text in sectionSpinnerText
			// and add to it the classes fa fa-cog and fa-spin - an action enough to display the turning gear
			// (see PS's CSS for more information)
			
			sectionSpinnerText = $('.main_div .main_section:not(.hidden) h2').text();
			$('.main_div .main_section:not(.hidden) h2').html(sectionSpinnerText + " <i class='fa fa-cog fa-spin'></i>");
		}
		sectionSpinnerOn = !sectionSpinnerOn;
	}
	
	
	var toggleNextClass = function (itms, className, mustnothaveClass, toggleFun) {
		
		// Getting an array of item ids (itms), toggleNextClass searches through the items for the
		// first item that does not have the classname [className]. The function, then, adds the classname to 
		// the item and removes it from the following one in the array.
		
		// If mustnothaveClass is set to false it does exactly the oppposite. it searches through the items for the
		// first to have the classname, removes it and adds it to the following one
		
		// when the one item that does not (or does) have the classname [className] is found and if toggleFun
		// is specified, toggleNextClass calls toggleFun(index) where index is a zero based index of the
		// element in itms
		
		// this allows us to perform validity checks before safely toggling next class and execute index-specific code
		// toggleNextClass will stop execution if toggleFun returns false
		
		// toggleNextClass is an excellent function for implementing step by step tab-like interfaces
		// PS typically uses togglenextclass to advance to further stages of analysis and hide current stage HTML divisions
		// while showing the next one.
		// It does that by toggling the "hiddenn" class in successive HTML divisions. toggleFun is "executeStage" in this case, and performs a validity check to make sure that
		// there is not any evident user input error before advancing to the next step
		
		// Comments in this function take for granted that mustnothaveClass is true from now on!
		
		
		var nToggle = 0;
		var toggleClassFunction = (mustnothaveClass ? ["removeClass", "addClass"] : ["addClass", "removeClass"]);
		
		
		itms.some(function (itm, index) {
			
			var cond = (mustnothaveClass ? $(itm).hasClass(className) : !$(itm).hasClass(className));
			
			if (cond) { // if the current item has the class [className]
				if (nToggle > 0) { // and if we have already found the one item that did not have the class [className] and sdded it
					$(itm)[toggleClassFunction[0]](className); // remove the class from this item
					nToggle++;
				}
				} else { // if the current item is the one that does not have the class [className]
				if (toggleFun != null)
				{ // if toggleFun is set then make sure it returns true for the specific element before toggling its class
					if (toggleFun(index))
					{
						$(itm)[toggleClassFunction[1]](className); // add the class from this item
						nToggle++;
					}
					} else { //if toggle fun is null simply add the class [className] to the element
					
					$(itm)[toggleClassFunction[1]](className); // add the class from this item
					nToggle++;
				}
			}
			return (nToggle === 2);
		});
		// Toggle 1st item again if toggled all in succession (valid when we have more than one items in the list)
		if (nToggle == 1 && itms.length > 1) {
			var itm = itms[0];
			var cond = (mustnothaveClass ? $(itm).hasClass(className) : !$(itm).hasClass(className));
			if (cond) {
				$(itm)[toggleClassFunction[0]](className);
				} else {
				$(itm)[toggleClassFunction[1]](className);
			}
		}
	}
	
	
	var uploadParameters = function() {
		
		// Bound to upload parameters hidden upload button
		if ($("#__upldparams").val() == "") {
			return;
		}
		var fnames = "";
		if (this.files.length > 1) return;
		if (this.files.length > 0) 
		{
			uploadingFile = this.files[0];
			// Check that the browser supports all necessary APIs
			if (window.File && window.FileReader && window.FileList && window.Blob) {
				// Determine all oversized files
				oversized_files_idxs = [];
				$.each(uploadingFiles, function (idx, file_i)
					{
						if(file_i.size > 2147483648)
						{
							oversized_files_idxs.push(idx);
						}
					});				
			}
			if(oversized_files_idxs.length == 0)
			{
				// Use upload_param_file to uplaod the parameter file to the server
				var thedata = new FormData();
				thedata.append('thefile', uploadingFile);
				thedata.append('session_id', sessionid);
				$.ajax({
					url: cgi_bin_path + 'upload_param_file.php', //Server script to receive file
					type: 'POST',
					// Form data
					data: thedata,
					//Options to tell jQuery not to process data or worry about content-type.
					cache: false,
					contentType: false,
					processData: false,
					}).done(function (data, textStatus, jqXHR) {
					// Load the parameters of the uploaded parameter file
					LoadParams(data.restext);
					$("#__upldparams").val("");
				});
			}
		}
	}
	
	
//}
//--------------------------------------
// HELPER ROUTINES __ END



//--------------------------------------
//$(document).ready is executed when the web page is loaded and initializes ProteoSign
//--------------------------------------

$(document).ready(ProteoSignInit);
