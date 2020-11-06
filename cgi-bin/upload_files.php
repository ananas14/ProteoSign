<?php
	error_reporting(E_ALL ^ E_WARNING);
	require 'get_labels.php';
	require 'get_rawfiles_names.php';
	
	
	
	// upload_files.php does not only upload the files to the server but at the same time finds the rawfiles contained in the uploaded file and
	// the labels (if any) that are contained in the file using get_rawfiles_names() and get_labels() from the corresponding php files that are loaded above.
	// upload_files also guesses the experiment type - that is if the experiment was Metabolically labelled (L), Isobarically labelled (IL) or Label Free (LF)
	
	// upload_files can be run using command line, it accepts the following arguments:
	// the first arg should always be "CL" standing for command line
	// the second contains the path of the file to read
	
	$CL_run = false;
	if (isset($argv[1]) && $argv[1] == "CL")
	{
		print("Upload files is run using Command Line\n");
		$_POST["server_side"] = "true";
		$_POST['thefile'] = $argv[2];
		$_POST["session_id"] = 1;
		$CL_run = true;
		$file_copied_successfully = true;
	}
	
	
	// PS assumes that some patterns in headers lead to a specific conclusion. For example if the headers contain "Spectrum File" it is assumed that this
	// file is a PD file since only PD files contain this kind of  (the uploaded file could be: MQ for MQ_evidence files, MQP for MQ_proteinGroups files and PD for PD psm files)
	//
	// This table shows which main patterns lead to which conlusion:
	
	//
	// +------------------------------+------------------------------------+
	// |            Pattern           |             Conclusion             |
	// +------------------------------+------------------------------------+
	// | Spectrum File                | PD file                            |
	// +------------------------------+------------------------------------+
	// | Raw file                     | MQ file                            |
	// +------------------------------+------------------------------------+
	// | Peptide IDs                  | MQP file                           |
	// +------------------------------+------------------------------------+
	// | [label]/[label]              | Label header in PD file            |
	// +------------------------------+------------------------------------+
	// | Ratio [label]/[label]        | Label header in MQ file            |
	// +------------------------------+------------------------------------+
	// | Reporter intensity [label]   | Label header in MQ file            |
	// +------------------------------+------------------------------------+
	// | Any header containing "file" | Raw files header in MQ or PD files |
	// +------------------------------+------------------------------------+
	// | Quan Channel                 | L experiment type in PD files      |
	// +------------------------------+------------------------------------+
	// | [digits]/[digits]            | IL experiment type in PD files     |
	// +------------------------------+------------------------------------+
	// | Labeling State               | L experiment type in MQ files      |
	// +------------------------------+------------------------------------+
	// | Reporter intensity           | IL experiment type in MQ files     |
	// +------------------------------+------------------------------------+
	//
	
	// _get_labels_aguments_sets contain unique regular expressions to get the labels from the header of the files. The resulting labels will be stored in "peptideLabelsNamesFromFile" in javascript
	// It is an array of arrays built so that it can provide alternative regexes for many filetypes providing the necessary labels. Notice that here the notation "L" stands for Isobaric labelling and labelling as well
	// e.g. an "L" file from PD might contain the labels in the form "{LABEL}\{LABEL}" or in the form "Abundance: {LABEL}". This is why if we expand _get_labels_aguments_sets we will get: TOP -> PD -> L -> {an array of two arrays the first element of which have regexes for the two
	// patterns discussed above} notice that the terminal nodes are always three siblings the first of which is the pattern for label recognition and the other 2 are obsolete
	
	// For more information about _get_labels_aguments_sets and unique_patterns see the comments below
	
	$_get_labels_aguments_sets = array(
	'PD' => array(
		'L' => array(
			array('/^([^\s]+)\/([^\s]+)$/i', '/Modifications/i', '/\((?:[^:]+?:)?(.+?)\)(;|$)/'), //Proteome Discoverer
			array('/Abundance: (\d+)/i', '/Modifications/i', '/\((?:[^:]+?:)?(.+?)\)(;|$)/'), //Proteome Discoverer
		)
	),
	'MQ' => array(
		'L' => array(
				array('/^Ratio ([^\s]+)\/([^\s]+)/i', null, null), //MaxQuant
				array('/^Reporter intensity ([0-9]+)/i', null, null), //MaxQuant, MS/MS-multiplexing (reporter ions), e.g. iTRAQ
			)
		)
	);
	
	$unique_patterns = array(
			'PD' => array(
				'L' => array('/Quan Channel/i'),
				'LF' => array(null), // there is no unique header for PD-LF
				'IL' => array('/\d+\/\d+/', 'Abundance: \d+')
		),
			'MQ' => array(
				'L' => array('/Labeling State/i'),
				'LF' => array(null), // there is no unique header for MQ-LF
				'IL' => array('/Reporter intensity/i')
		)
	);
	
	$server_side_file = ($_POST["server_side"] === "true");
	
	if ($server_side_file) {
		$name = $_POST['thefile'];
		$tmp_name = $_POST['thefile'];
		} else {
		$name = $_FILES['thefile']['name'];
		$tmp_name = $_FILES['thefile']['tmp_name'];
		
	}
	$server_response = [];
	//Overall success, will be checked by client.
	$server_response['success'] = false;
	$server_response['msg'] = '';
	$server_response['mkdir_msg'] = '';
	$server_response['peptide_labels'] = [];
	$server_response['peptide_labels_names'] = [];
	$server_response['skipped_labels'] = [];
	$server_response['raw_filesnames'] = [];
	$server_response['file_type'] = '';
	$server_response['exp_type'] = '';
	$server_response['ret_session'] = $_POST["session_id"];
	$document_root = dirname(__DIR__);
	
	if (isset($name)) {
		if (!empty($name)) {
			
			if (!$CL_run)
			{
				// Create a new session directory
				$location = $document_root . '/uploads/' . $_POST["session_id"];
				if (!file_exists($location) && !is_dir($location)) {
					if (!mkdir($location, 0777, true)) {
						$server_response['mkdir_msg'] = "The directory $location already exists, probably created from another upload_files instance.";
						//Sometimes upload files tries to create an existing directory even if file_exists($location) returns true, check again if the file exists and display the error if so
						if (!file_exists($location) && !is_dir($location)) {
							$server_response['msg'] = "The directory $location could not be created ('mkdir' returned FALSE).";
							goto end;
						}
					}
				}
				
				// Copy the file to the session folder: if the file is server sided
				// simply copy it from test data folder. If ot was manually uploaded
				// from the user, apache server uses a temp_name to store it and move_uploaded_file is the
				// dedicated php function to copy it to another folder - here the session folder
				$file_copied_successfully = false;
				if ($server_side_file)
				{
					$file_copied_successfully = copy($document_root . "/test data/" . $tmp_name, $location . '/' . $name);
				}
				else
				{
					if (!file_exists($location . "/" . $name))
					{
						$file_copied_successfully = move_uploaded_file($tmp_name, $location . "/" . $name);
					}
					else
					{
						$server_response['success'] = false;
						$server_response['msg'] = "A file with the same name ($name) has already been uploaded, try again later";
						goto end;
					}
				}
			}
			
			if (!$CL_run)
			{
				$full_path = $location . '/' . $name;
			}
			else
			{
				$full_path = $name;
			}
			
			$handle = null;
			if ($file_copied_successfully && ($handle = fopen($full_path, "r")))
			{
				
				
				// First decide if the processing program is PD or MQ and thye filetype (MQ, MQP or PD)
				// Get the headers of the file
				$first_line = fgetcsv($handle, 0, "\t");
				fclose($handle);
				// The file could be: MQ for MQ_evidence files, MQP for MQ_proteinGroups files and PD for PD psm files
				// if the headers contain "Spectrum File" it is from proteome discoverer. If they contain "Raw file" it is an MQ
				// file and if it contains "Peptide IDs" it is an MQP file
				
				$tmp = preg_grep("/Spectrum File/i", $first_line);
				if (count($tmp) > 0) {
					$dtype = 'PD';
					} else {
					$tmp = preg_grep("/Raw file/i", $first_line);
					if (count($tmp) > 0){
						$dtype = 'MQ';
					}
					else{
						$tmp = preg_grep("/Peptide IDs/i", $first_line);
						if (count($tmp) > 0){
							$dtype = 'MQP';
						}
						else{
							$dtype = 'unknown';
						}
					}
				}
				
				// Now according to the filetype, search the headers against unique patterns to find the experiment type (L for labelled - IL for isobarically labelled and LF for label free)
				// the unique_patterns variable is an array of three dimensions. The first is the processing program, the second is a filetype and the third are the patterns.
				// for example unique_patterns['PD']['L'] = array("/^Quan Channel$/") that means that this pattern can be found only in Metabolically labelled datasets from all datasets derived from PD
				// Notice that there are no unique patterns for LFQ datasets in both MQ and PD so if no unique patterns for Metabolically labelled and Isobarically labelled datasets are found
				// the dataset will be considered to be an LFQ dataset
				
				$filetype = "";
				if ($dtype == 'MQ' || $dtype == 'PD')
				{
					foreach ($unique_patterns[$dtype]['L'] as &$Lpattern)
					{
						$tmp = preg_grep($Lpattern, $first_line);
						if (count($tmp) > 0)
						{
							$filetype = "L";
						}
					}
					if ($filetype == "")
					{
						foreach ($unique_patterns[$dtype]['IL'] as &$ILpattern)
						{
							$tmp = preg_grep($ILpattern, $first_line);
							if (count($tmp) > 0)
							{
								$filetype = "IL";
							}
						}
						if ($filetype == "")
						{
							$filetype = "LF";
						}
					}				
				}
				$server_response['exp_type'] = $filetype;
				
				// If it is not MQP file try to find the labels contained (if it is a labelled dataset)
				if ($dtype == 'MQ' || $dtype == 'PD')
				{
					$server_response['file_type'] = $dtype;
					
					// get_labels_aguments_sets is an array three dimensions the first dimension corresponds to different processing programs
					// the second one to different filetypes (at the moment only L - Labelled filetypes are supported but in the future more specific triads might be introduced responding to e.g. LFQ or isobarically labelled data)
					// and the third one is a pattern triad - processing program and filetype specific:
					//
					// Each triad of patterns corresponds to: (1) a pattern of headers containing labels
					// (2) a pattern for headers containing the peptide modifications and (3) a pattern to discover all modifications applied to the peptides in this experiment
					//
					// notice that in variables names the modifications are also called "label definitions"
					//
					// For example for PD: get_labels_aguments_sets['PD'] might get the following values:
					/*
						Array
						(
						[0] => Array
						(
						[0] => /^([^\s]+)\/([^\s]+)$/i               (label header pattern)
						[1] => /Modifications/i                      (label definition header pattern)
						[2] => /\((?:[^:]+?:)?(.+?)\)(;|$)/          (label definition pattern)
						)
						)
					*/
					
					// for example the array [0] in the example above would match:
					// [0][0]: the header "Medium/Light"
					// [0][1]: the header "Modifications"
					// [0][2]: the value  "N-Term(Dimethyl)" under the column Modifications which is a valid "label definition"
					// this way if this triad was fed to a function to test all patterns above it would be able to find all labels, the label-definition column and all valid label definitions in this column
					// this is exactly what get_labels(...) does
					//
					// ProteoSign does not currently use label definitions, however since this might be a very useful future feature we keep a get_labels definition that is able to find all
					// modifications but keep the finding operation commented out
					// In case you want to add a new pattern triad to fetch the labels from the file you can add something like ('/[pattern]/', null, null) to avoid label definition searching
					
					$get_labels_aguments_sets = array_values($_get_labels_aguments_sets[$dtype]['L']);
					
					// Keep calling "get_labels" until you get some labels
					for ($i = 0; $i < count($get_labels_aguments_sets); $i++) {
						$argset = $get_labels_aguments_sets[$i];
						// For example an argset that targets a labelled dataset may be
						/*
							Array
							(
							[0] => /^([^\s]+)\/([^\s]+)$/i               (label header pattern)
							[1] => /Modifications/i                      (label definition header pattern)
							[2] => /\((?:[^:]+?:)?(.+?)\)(;|$)/          (label definition pattern)
							)
						*/
						
						$tmp = get_labels($full_path, $argset[0], $argset[1], $argset[2]);
						
						// tmp is an array of two arrays.
						// tmp[0] contains all labels detected in the file e.g.
						/*
							Array
							(
							[0] => Medium
							[1] => Light
							)
						*/
						// and tmp[1] is null
						
						// The following block checks for the validity of label definitions, comment in in case you want to use label definitions in ProteoSign:
						/*
							if (count($tmp[0]) > 0)
							{
							// If get_labels returned at least one label definition
							if (count($tmp[1]) > 0)
							{
							$okdefs = 0;
							// For each label definition:
							foreach ($tmp[1] as $lbldef)
							{
							// Make sure that the label definition does not contain space characters,
							// commas, semicolons and colons
							if (preg_match('/[\s,;\:]/i', $lbldef) == 0)
							{
							$okdefs++;
							}
							else
							{
							array_push($server_response['skipped_labels'], $lbldef);
							}
							}
							// If all definitions are OK stop calling get_labels
							if ($okdefs == count($tmp[1]))
							{
							break;
							}
							}
							else
							{
							// If no definitions were returned stop calling get_labels
							break;
							}
							}
						*/
						// If we found labels just break testing patterns
						if (count($tmp[0]) > 0)
						{
							break;
						}
					}
					
					// peptide_labels_names is the variable that will be sent back to JS containing all labels. peptide_labels is another variable left for
					// returning modifications in case we uncomment the respective feature in get_labels
					
					$server_response['peptide_labels_names'] = $tmp[0];
					$server_response['peptide_labels'] = $tmp[1];
					
					// In a similar manner, get_rawfiles_names gets all raw_files contained in the file: its definition is almost identical to this of get_labels
					
					$server_response['raw_filesnames'] = get_rawfiles_names($full_path, '/file/i');
					
					
					if (count($server_response['raw_filesnames']) == 0 && !$CL_run)
					{
						error_log("[client: " . $_SERVER['REMOTE_ADDR'] . "] Could not retrieve replicate information (raw files names) from data file " . $name);
					}
					// Rename the files to msdiffexp_peptide.txt fro PD and MQ files and to msdiffexp_protein.txt for MQP files
					if (!$CL_run)
					{
						rename($location . '/' . $name, $location . '/msdiffexp_peptide.txt');
					}
				}
				elseif ($dtype == 'MQP')
				{
					$server_response['file_type'] = $dtype;
					if (!$CL_run)
					{
						rename($location . '/' . $name, $location . '/msdiffexp_protein.txt');
					}
				}
				else
				{
					$server_response['file_type'] = "unknown";
					if (!$CL_run)
					{
						unlink($location . '/' . $name);
					}
					$server_response['success'] = true;
					$server_response['msg'] = "The file $name is not valid";
					goto end;
				}
				$server_response['success'] = true;
			}
			else
			{
				$server_response['msg'] = "The file $name could not be moved ('move_uploaded_file' returned FALSE).";
			}
		}
		else
		{
			$server_response['msg'] = "The variable 'name' was empty (empty() returned TRUE).";
		}
	}
	else
	{
		$server_response['msg'] = "The variable $name was not set ('isset' returned FALSE).";
	}
	
	end:
	if (!$CL_run)
	{
		error_log("[client: " . $_SERVER['REMOTE_ADDR'] . "] upload_files.php [" . $_POST["session_id"] . " " . $name . " ]> Success: " . ($server_response['success'] ? 'Yes' : 'No') . " | Message: " . $server_response['msg']);
		//Send info back to the client
		header('Content-type: application/json');
		echo json_encode($server_response);
	}
	else
	{
		// In case the program was run by command line print the necessary results back to the user
		print("Filetype:\n");
		print($server_response['file_type']);
		print("\nExperiment type:\n");
		print($server_response['exp_type']);
		print("\n\nPeptide label names from file:\n");
		print_r($server_response['peptide_labels_names']);
		print("\n\nRaw files detected:\n");
		print_r($server_response['raw_filesnames']);
	}
?>