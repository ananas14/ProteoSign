<?php
	//
	// OBSOLETE FILE
	//
	
	$document_root = dirname(__DIR__);
	$location = $document_root . '/uploads/' . $_POST["session_id"];
	// //DEBUG:
	// $location = "C:/xampp/htdocs/ProteoSign/uploads/1477755182769";
	$name = 'msdiffexp_peptide.txt';
	$data_file = $location . '/' . $name;
	$server_response['labels'] = [];
	$server_response['isIsobaricLabel'] = false;
	$server_response['success'] = false;
	$server_response['msg'] = "";
	$server_response['special_IL_labels'] = [];
	$patternfound = false;
	if (!file_exists($location) && !is_dir($location)) {
		$server_response['msg'] = "The directory $location could not be found! at get_labels_to_js.";
		goto end;
	}
	// $temp = "";
	if (($handle = fopen($data_file, "r")) !== FALSE) {
		$titles = explode("\t", fgets($handle));
		$pattern = '/^Intensity .*/i';
		$pattern2 = '/^Reporter intensity/i';
		$pattern3 = '/^\d*\/\d*$/';
		$pattern4 = '/^\d+[a-zA-Z]?\s?$/';
		$pattern5 = '/^Abundance: \d*/';
		$counter = 0;
		$pattern4counter = 0;
		$pattern5counter = 0;
		foreach ($titles as &$value) {
			// $temp .= $value;
			$value = preg_replace("~^\"~", "", $value);
			$value = preg_replace("~\n$~", "", $value);
			$value = preg_replace("~\r$~", "", $value);
			$value = preg_replace("~\"$~", "", $value);
			preg_match($pattern, $value, $matches);
			if (count($matches) > 0)
			{
				$label_name = substr($matches[0], 10);
				if(preg_match('/\+/', $label_name, $my_matches) == 0)
				{
					$server_response['labels'][$counter] = $label_name;
					$counter++;
				}
			}
			preg_match($pattern2, $value, $matches);
			if (count($matches) > 0)
			{
				$server_response['isIsobaricLabel'] = true;
			}
			preg_match($pattern3, $value, $matches);
			if (count($matches) > 0)
			{
				$server_response['isIsobaricLabel'] = true;
				$patternfound = true;
			}
			preg_match($pattern4, $value, $matches);
			if (count($matches) > 0)
			{
				if ($patternfound == false)
				{
					//here a special case of Isobaric Label (PD) data is present, the headers do not contain the (tag/tag e.g 125/126) header but only contain the (tag e.g. 125 etc.) header. In this case return isIsobaricLabel true and return a special variable (special_IL_labels) with the tags:
					//notice that if some conditions have their own header and then a ratio is present as header (e.g. 125 126 125/126) the special_IL_labels will only contain the first two instances (125, 126). Nevertheless, a ratio will have been found and the special_IL_labels variable will not be used
					$server_response['isIsobaricLabel'] = true;
					$matches[0] = preg_replace("/\s$/", "", $matches[0]);
					$server_response['special_IL_labels'][$pattern4counter] = $matches[0];
					$pattern4counter++;
				}
			}
			preg_match($pattern5, $value, $matches);
			if (count($matches) > 0)
			{
				$server_response['isIsobaricLabel'] = true;
				if ($patternfound == false)
				{
					$server_response['special_IL_labels'] = [];
				}
				$patternfound = true;
				$server_response['isIsobaricLabel'] = true;
				$matches[0] = preg_replace("/^Abundance: /", "", $matches[0]);
				$server_response['special_IL_labels'][$pattern5counter++] = $matches[0];
			}
		}
		fclose($handle);
		$server_response['success'] = true;
	}
	end:
	// $server_response['msg'] = $temp;
	error_log("[client: " . $_SERVER['REMOTE_ADDR'] . "] get_labels_to_js.php [" . $_POST["session_id"] . " " . $name . " ]> Success: " . ($server_response['success'] ? 'Yes' : 'No') . " | Message: " . $server_response['msg']);
	//Send info back to the client
	header('Content-type: application/json');
	echo json_encode($server_response);
?>