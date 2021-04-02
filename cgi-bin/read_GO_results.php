 <?php
 

 $document_root = dirname(__DIR__);
	$session_folder = $_POST["session_id"];
	$upload_dir = $document_root . "/uploads/" . $session_folder;
	$path = $upload_dir . "/msdiffexp_wd/msdiffexp_out/" . $_POST["GOfile"];

	$server_response = [];
	$server_response['GOdata'] = [];
	
	# Open the file and get rid of the header
	$myfile = fopen($path, "r");
	fgets($myfile);
	$counter = 1;
	
	
	
	# For each line break it in tabs and set a value for each cell
	while(!feof($myfile)) {
	  $line = fgets($myfile);
	  $tabs = explode("\t", "$line\t\t\t\t\t\t\t\t");
	  $temp = [];
	  $temp[] = ['id'=>$counter++, 'source'=>$tabs[0], 'function'=>$tabs[1], 'term_id'=>$tabs[2], 'p_value'=>$tabs[3], 'term_size'=>$tabs[4], 'query_size'=>$tabs[5], 'intersection_size'=>$tabs[6], 'intersection'=>$tabs[7]];
  
	  
	  array_push($server_response['GOdata'], $temp);
	}
	fclose($myfile);
	$server_response['msg'] = "";
	$server_response['success'] = true;
	
	error_log("[client: " . $_SERVER['REMOTE_ADDR'] . "] read_GO_results.php [" . $_POST["session_id"] . " ]> Success: " . ($server_response['success'] ? 'Yes' : 'No') . " | Message: " . $server_response['msg']);
	//Send info back to the client
	header('Content-type: application/json');
	echo json_encode($server_response);
 ?>