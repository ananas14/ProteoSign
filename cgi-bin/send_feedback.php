<?php
	/*
		This php file opens the Feedbacks tab-seperated file in cgi-bin and appends the sessionid, the current time and the feedback of the user
	*/
	
	function wp_normalize_path( $path ) {
		$path = str_replace( '\\', '/', $path );
		$path = preg_replace( '|(?<=.)/+|', '/', $path );
		if ( ':' === substr( $path, 1, 1 ) ) {
			$path = ucfirst( $path );
		}
		return $path;
	}
	
	
	// Make sure that the file will be saved in a subdirectory of the uploads folder
	
	$document_root = dirname(__DIR__);
	
	$document_root .= "/uploads";
	$file_location = $document_root . "/" . $_POST['feedback_file'];
	$document_root = realpath($document_root);
	$document_root = wp_normalize_path($document_root);
	$realfilepath = realpath("$file_location");
	$realfilepath = wp_normalize_path($realfilepath);


	if ($document_root != substr($realfilepath, 0, strlen($document_root)))
	{
		$server_response['success'] = false;
		$server_response['msg'] = "The file location was not valid!";
		goto end;
	}
	
	
	$file_location = $realfilepath;
	$server_response['success'] = true;
	$server_response['msg'] = '';
	$user_feedback = htmlentities($_POST["texttoappend"]);
	$user_feedback = preg_replace("/\t/", "   ", $user_feedback);
	$user_feedback = preg_replace("/\n/", " <br> ", $user_feedback);
	$user_feedback = preg_replace("/\r/", "", $user_feedback);
	$texttowrite = date("l d/m/Y H:i:s") . "\t" . $_POST["session_id"] . "\t" . $_POST["stars_score"] . "\t" . $user_feedback . "\tNo\n";
	if ($ff = fopen($file_location, 'a'))
	{
		$canwrite = fwrite($ff, $texttowrite);
		if(!$canwrite){
			$server_response['msg'] = "The file $file_location could not be written ('fwrite' returned FALSE)";
		}
		if($canwrite && !fclose($ff)){
			$server_response['msg'] = "The file $file_location could not be closed ('fclose' returned FALSE)";
			}else{
			$server_response['success'] = true;
		}
	}
	end:
	error_log("[client: " . $_SERVER['REMOTE_ADDR'] . "] send_feedback.php [" . $_POST["session_id"] . " ]> Success: " . ($server_response['success'] ? 'Yes' : 'No') . " | Message: " . $server_response['msg']);
	//Send info back to the client
	header('Content-type: application/json');
	echo json_encode($server_response);
?>