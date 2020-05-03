<?php
	error_reporting(E_ALL ^ E_WARNING);
	$document_root = dirname(__DIR__);
	$server_response['text'] = '';
	$location = $document_root . '/uploads/' . $_POST["session_id"];
	$name = $_FILES['thefile']['name'];
	$tmp_name = $_FILES['thefile']['tmp_name'];
	$file_copied_successfully = move_uploaded_file($tmp_name, $location . "/" . $name);
	$server_response['success'] = $file_copied_successfully;
	$server_response['restext'] = file_get_contents ($location . "/" . $name);
	header('Content-type: application/json');
	echo json_encode($server_response);
?>