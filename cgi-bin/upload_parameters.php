 <?php
	$server_response = [];
	$server_response['success'] = false;
	$server_response['msg'] = "";
	define('DS', DIRECTORY_SEPARATOR);
	$session_folder = $_POST["session_id"];
	$document_root = dirname(__DIR__);
	$upload_dir = $document_root . DS . "uploads" . DS . $session_folder;
	$upload_parameter_file = $upload_dir . DS . "MSdiffexp_definitions.R";
	$upload_experimental_structure_file = $upload_dir . DS . "exp_struct.txt";
	$upload_LFQ_data_file = $upload_dir . DS . "LFQ_conditions.txt";
	$upload_Rename_Array_file = $upload_dir . DS . "Rename_array.txt";
	$upload_LS_file = $upload_dir . DS . "LS_array.txt";
	$upload_RMrawfiles_file = $upload_dir . DS . "RMrawfiles.txt";
	$upload_RMtags_file = $upload_dir . DS . "RMtags.txt";
	$parameters_template = "parameters_template.R";
	$the_parameters["REPLACE1"] = $_POST["exppddata"];
	$the_parameters["REPLACE2"] = "\"" . $_POST["expid"] . "\"";
	$the_parameters["REPLACE3"] = "\"" . $_POST["exptpoint"] . "\"";
	$the_parameters["REPLACE4"] = $_POST["expprotquant"];
	$the_parameters["REPLACE5"] = $_POST["expquantfilt"];
	$the_parameters["REPLACE6"] = $_POST["expquantfiltprot"];
	$the_parameters["REPLACE7"] = "\"" . $_POST["expquantfiltlbl"] . "\"";
	$the_parameters["REPLACE8"] = $_POST["labelfree"];
	$the_parameters["REPLACE9"] = "\"" . $upload_dir . '\\msdiffexp_wd' . "\"";
	$the_parameters["REPLACE10"] = $_POST["IsIsobaricLabel"];
	$the_parameters["REPLACE11"] = $_POST["All_MQ_Labels"];
	$the_parameters["REPLACE12"] = $_POST["AllowMergeLabels"];
	$the_parameters["REPLACE13"] = $_POST["AllowLS"];
	$the_parameters["REPLACE14"] = $_POST["LeastBreps"];
	$the_parameters["REPLACE15"] = $_POST["LeastPeps"];
	$the_parameters["REPLACE16"] = $_POST["PThreshold"];
	$the_parameters["REPLACE17"] = $_POST["RMisused"];
	$the_parameters["REPLACE18"] = $_POST["RMbrepsRepInRawFiles"];
	$the_parameters["REPLACE19"] = $_POST["RMtrepsRepInRawFiles"];
	$the_parameters["REPLACE20"] = $_POST["RMconditionsRepInRawFiles"];
	$the_parameters["REPLACE21"] = $_POST["GOenrichment"];
	// Especially for the goorganism parameter we need to maotch the organism description to the corresponding organism ID contained
	// in the valid_GO_Organisms.txt file. Skim through the text file and find the corresponding ID
	$the_parameters["REPLACE22"] = "";
	if($ff = fopen("valid_GO_Organisms.txt", 'r')){
		while(!feof($ff)) {
			$line = fgets($ff);
			$line = rtrim($line);
			$myarray = explode("\t", $line);
			if ($myarray[0] == $_POST["GOorganism"])
			{
				$the_parameters["REPLACE22"] = '"' . $myarray[1] . '"';
				break;
			}
		}
		if ($the_parameters["REPLACE22"] == "")
		{
			$server_response['msg'] = "No valid organism found!";
			$server_response['success'] = false;
			goto end;
		}
	}
	else
	{
		$server_response['msg'] = "Could not open the valid GO organisms list file!";
		$server_response['success'] = false;
		goto end;
	}
	
	if($_POST["explbl00"] == "T" && isset($_POST["explbl0"]) && strlen($_POST["explbl0"]) > 0){
		$the_parameters["APPEND0"] = "addLabel(\"" . $_POST["explbl0"] . "\",c(\"\"))";
	}
	$lbl_i = 1;
	while(isset($_POST["explbl" . $lbl_i . "name"])){
		$the_parameters["APPEND" . $lbl_i] = "addLabel(\"" . $_POST["explbl" . $lbl_i . "name"] . "\",c(\"\"))";
		$lbl_i++;
	}
	$parameters_file_contents = file_get_contents("$parameters_template") ;
	if($parameters_file_contents){
		while (list($key, $param_value) = each($the_parameters)) {
			if(preg_match('/^REPLACE/',$key)){
				$parameters_file_contents = preg_replace('/' . $key . "([[:space:]])" .  '/', $param_value . "$1", $parameters_file_contents);
			}else if(preg_match('/^APPEND/',$key)){
				$parameters_file_contents = $parameters_file_contents . $param_value . PHP_EOL;
			}
		}
	}else{
		$server_response['msg'] = "The file $parameters_template could not be read ('file_get_contents' returned FALSE)";
	}

	if($ff = fopen($upload_parameter_file, 'w')){
		$canwrite = fwrite($ff, $parameters_file_contents);
		if(!$canwrite){
			$server_response['msg'] = "The file $upload_parameter_file could not be written ('fwrite' returned FALSE)";
		}
		if($canwrite && !fclose($ff)){
			$server_response['msg'] = "The file $upload_parameter_file could not be closed ('fclose' returned FALSE)";
		}else{
			$server_response['success'] = true;
		}
	}else{
		$server_response['msg'] = "The file $upload_parameter_file could not be opened ('fopen' returned FALSE)";
	}
	
	if($server_response['success']){
		$server_response['success'] = false;
		if($ff = fopen($upload_experimental_structure_file, 'w')){
			$canwrite = fwrite($ff, $_POST["exp_struct"]);
			if(!$canwrite){
				$server_response['msg'] = "The file $upload_experimental_structure_file could not be written ('fwrite' returned FALSE)";
			}
			if($canwrite && !fclose($ff)){
				$server_response['msg'] = "The file $upload_experimental_structure_file could not be closed ('fclose' returned FALSE)";
			}else{
				$server_response['success'] = true;
			}
		}else{
			$server_response['msg'] = "The file $upload_experimental_structure_file could not be opened ('fopen' returned FALSE)";
		}
		if($ff = fopen($upload_LFQ_data_file, 'w')){
			$canwrite = fwrite($ff, $_POST["LFQ_conds"]);
			if(!$canwrite){
				$server_response['msg'] = "The file $upload_LFQ_data_file could not be written ('fwrite' returned FALSE)";
			}
			if($canwrite && !fclose($ff)){
				$server_response['msg'] = "The file $upload_LFQ_data_file could not be closed ('fclose' returned FALSE)";
			}else{
				$server_response['success'] = true;
			}
		}else{
			$server_response['msg'] = "The file $upload_LFQ_data_file could not be opened ('fopen' returned FALSE)";
		}	
		if ($_POST["AllowMergeLabels"] == "T")
		{
			if($ff = fopen($upload_Rename_Array_file, 'w')){
				$canwrite = fwrite($ff, $_POST["Rename_Array"]);
				if(!$canwrite){
					$server_response['msg'] = "The file $upload_Rename_Array_file could not be written ('fwrite' returned FALSE)";
				}
				if($canwrite && !fclose($ff)){
					$server_response['msg'] = "The file $upload_Rename_Array_file could not be closed ('fclose' returned FALSE)";
				}else{
					$server_response['success'] = true;
				}
			}else{
				$server_response['msg'] = "The file $upload_Rename_Array_file could not be opened ('fopen' returned FALSE)";
			}	
		}
		if ($_POST["AllowLS"] == "T")
		{
			if($ff = fopen($upload_LS_file, 'w')){
				$canwrite = fwrite($ff, $_POST["LabelSwapArray"]);
				if(!$canwrite){
					$server_response['msg'] = "The file $upload_LS_file could not be written ('fwrite' returned FALSE)";
				}
				if($canwrite && !fclose($ff)){
					$server_response['msg'] = "The file $upload_LS_file could not be closed ('fclose' returned FALSE)";
				}else{
					$server_response['success'] = true;
				}
			}else{
				$server_response['msg'] = "The file $upload_LS_file could not be opened ('fopen' returned FALSE)";
			}	
		}
		if ($_POST["RMisused"] == 'T')
		{
			if($ff = fopen($upload_RMrawfiles_file, 'w')){
				$canwrite = fwrite($ff, $_POST["RMRawFiles"]);
				if(!$canwrite){
					$server_response['msg'] = "The file $upload_RMrawfiles_file could not be written ('fwrite' returned FALSE)";
				}
				if($canwrite && !fclose($ff)){
					$server_response['msg'] = "The file $upload_RMrawfiles_file could not be closed ('fclose' returned FALSE)";
				}else{
					$server_response['success'] = true;
				}
			}else{
				$server_response['msg'] = "The file $upload_RMrawfiles_file could not be opened ('fopen' returned FALSE)";
			}
			if($ff = fopen($upload_RMtags_file, 'w')){
				$canwrite = fwrite($ff, $_POST["RMTags"]);
				if(!$canwrite){
					$server_response['msg'] = "The file $upload_RMtags_file could not be written ('fwrite' returned FALSE)";
				}
				if($canwrite && !fclose($ff)){
					$server_response['msg'] = "The file $upload_RMtags_file could not be closed ('fclose' returned FALSE)";
				}else{
					$server_response['success'] = true;
				}
			}else{
				$server_response['msg'] = "The file $upload_RMtags_file could not be opened ('fopen' returned FALSE)";
			}
		}
	}
end:
	error_log("[client: " . $_SERVER['REMOTE_ADDR'] . "] upload_parameters.php [" . $_POST["session_id"] . "]> Success: " . ($server_response['success'] ? 'Yes' : 'No') . " | Message: " . $server_response['msg']);
	//Send info back to the client
	header('Content-type: application/json');
	echo json_encode($server_response);	
?>
        
       