 <?php
	/*
	// DEBUG
	$test_data_dir = "C:/Users/G Efstathiou/Documents/GitHub/ProteoSign/test data/";
	$descriptions_requested = false;
	$dataset_requested = "dimethyl 2-plex (PD)";
	*/
	$zipname = preg_replace("/[\s\?\\\$\:\(\):<>\|\*\/\"\']+/","_",$_GET["dataset_info_requested"]) . '.zip';
	$test_data_dir = dirname(__DIR__) . "/test_data";
	
	$server_response = [];
	$server_response['success'] = false;
	$server_response['msg'] = "";
	try{
		$db = new SQLite3("$test_data_dir/testdatadb");
	}
	catch (Exception $exception) {
		$server_response['msg'] = "Failed to connect to the database: " . $exception->getMessage();
		goto end;
	}

	$qres = $db->query('select file from files inner join dataset_files on dataset_files.file_id=files.id inner join dataset on dataset.id=dataset_files.dataset_id where dataset.desc="' .  $_GET["dataset_info_requested"] . '"');
	
	$requested_files = [];
    while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
		foreach($qrow as $col => $val){
			$requested_files[] = $val;
		}
	}
	$qres->finalize();
	$db->close();
		
	$filename = "$test_data_dir/" . $zipname;
	clearstatcache();
	if(! file_exists($filename)){
		$zip = new ZipArchive();
		if ($zip->open($filename, ZipArchive::CREATE) !== TRUE) {
			$server_response['msg'] = 'Failed to create archive "' . $filename . '".';
			goto end;
		}
		foreach($requested_files as $file_i) {
			$zip->addFile("$test_data_dir/" . $file_i, $file_i);
		}
		$zip->close();
	}
	$server_response['success'] = true;
 end:
	header('Content-Description: File Transfer');
	header('Content-Type: application/octet-stream');
	header('Content-Disposition: attachment; filename=' . $zipname);
	header('Expires: 0');
	header('Cache-Control: must-revalidate');
	header('Pragma: public');
	header('Content-Length: ' . filesize($filename));
	ob_clean();
    flush();	
	readfile($filename);
	exit;
?>
