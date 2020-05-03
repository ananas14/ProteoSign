 <?php
	/*
		This php file appends data to a log file in the session folder
	*/
	$document_root = dirname(__DIR__);
	$location = $document_root . '/uploads/' . $_POST["session_id"];
	$file_location = $location . "/php_log.txt";
	$server_response['success'] = true;
	$server_response['msg'] = '';
	$texttoappend = $_POST["texttoappend"];
	if ($ff = fopen($file_location, 'a'))
	{
		$canwrite = fwrite($ff, $texttoappend);
		if(!$canwrite){
			$server_response['msg'] = "The file $file_location could not be written ('fwrite' returned FALSE)";
		}
		if($canwrite && !fclose($ff)){
			$server_response['msg'] = "The file $file_location could not be closed ('fclose' returned FALSE)";
		}else{
			$server_response['success'] = true;
		}
	}
	error_log("[client: " . $_SERVER['REMOTE_ADDR'] . "] log_php.php [" . $_POST["session_id"] . " ]> Success: " . ($server_response['success'] ? 'Yes' : 'No') . " | Message: " . $server_response['msg']);
	//Send info back to the client
	header('Content-type: application/json');
	echo json_encode($server_response);
 ?>