<?php
	// This file gets all valid GO enrichment organisms for gprofiler from a server file
	$session_folder = $_POST["session_id"];
	$server_response['validOrganisms'] = [];
	$document_root = dirname(__DIR__);
	$myfile = fopen($document_root . "/cgi-bin/valid_GO_Organisms.txt", "r");
	while(!feof($myfile)) {
		$row = fgets($myfile);
		$cells = explode("\t", $row);
		// the table valid_GO_Organisms has two columns the first has all valid organisms and the second matches each organism to an id for gprofiler
		// we obviously need all the values from the first column
		if ($cells[0] != "") array_push($server_response['validOrganisms'], $cells[0]);
	}
	fclose($myfile);
	$server_response['success'] = true;
	header('Content-type: application/json');
	echo json_encode($server_response);
?>