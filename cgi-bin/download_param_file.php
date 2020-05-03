<?php
	$session_folder = $_POST["session_id"];
	$createOutputParam_file = $_POST["output"];//1 for T 0 for F
	$document_root = dirname(__DIR__);
	$upload_dir = $document_root . "/uploads/" . $session_folder;
	$temp_path = "uploads/" . $session_folder;
	$experiment_name = $_POST["expname"];
	$my_text = $_POST["texttodownload"];
	if ($createOutputParam_file == 0)
	{
		$myfile = fopen($upload_dir . "/" .  $experiment_name . "_parameters.txt", "w");
	}
	else
	{
		$myfile = fopen($upload_dir . "/" .  $experiment_name . "_parameters_from_session_" . $session_folder . ".txt", "w");
	}
	fwrite($myfile, $my_text);
	fclose($myfile);
	$server_response['success'] = true;
	$server_response['results_url'] = $temp_path . "/" .  $experiment_name . "_parameters.txt";
	header('Content-type: application/json');
echo json_encode($server_response);
?>