 <?php
	/*
		This php file opens the Feedbacks tab-seperated file in cgi-bin and appends the sessionid, the current time and the feedback of the user
	*/
	$document_root = dirname(__DIR__);
	//WARNING! the line below is hardcoded for security reasons
	$file_location = $_POST['feedback_file'];
	$server_response['success'] = true;
	$server_response['msg'] = '';
	$user_feedback = $_POST["texttoappend"];
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
	error_log("[client: " . $_SERVER['REMOTE_ADDR'] . "] send_feedback.php [" . $_POST["session_id"] . " ]> Success: " . ($server_response['success'] ? 'Yes' : 'No') . " | Message: " . $server_response['msg']);
	//Send info back to the client
	header('Content-type: application/json');
	echo json_encode($server_response);
 ?>