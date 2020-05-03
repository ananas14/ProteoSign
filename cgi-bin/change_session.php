 <?php
	 /*
		This php file makes a copy of a session to a new one so that the same files can be run 
	 */
	 $document_root = dirname(__DIR__);
	 $old_location = $document_root . '/uploads/' . $_POST["old_session_id"];
	 $new_location = $document_root . '/uploads/' . $_POST["session_id"];
	 $server_response['success'] = true;
	 $server_response['msg'] = 'old session: ' . $_POST["old_session_id"] . ' new session: ' . $_POST["session_id"] . ' ';
	 if (!file_exists($new_location))
	 {
		 mkdir($new_location, 0777, true);
	 }
	 if (file_exists($old_location . '/msdiffexp_wd/msdiffexp_peptide.txt'))
	 {
		 copy($old_location . '/msdiffexp_wd/msdiffexp_peptide.txt', $new_location . '/msdiffexp_peptide.txt');
	 }
	 else
	 {
		 //if we did not find the file in its appropriate location it might have been left behind, one folder up
		 if (file_exists($old_location . '/msdiffexp_peptide.txt'))
		 {
			 copy($old_location . '/msdiffexp_peptide.txt', $new_location . '/msdiffexp_peptide.txt');
		 }
		 else
		 {
			$server_response['success'] = false;
			$server_response['msg'] .= 'msdiffexp_peptide not found ';
		 }

	 }
	 if (file_exists($old_location . '/msdiffexp_wd/msdiffexp_protein.txt'))
	 {
		 copy($old_location . '/msdiffexp_wd/msdiffexp_protein.txt', $new_location . '/msdiffexp_protein.txt');
	 }
	 else{
		if (file_exists($old_location . '/msdiffexp_protein.txt'))
		{
			copy($old_location . '/msdiffexp_protein.txt', $new_location . '/msdiffexp_protein.txt');
		}
		else
		{
			$server_response['msg'] .= 'msdiffexp_protein not found, probably working on a PD dataset ';
		}
	 }
	 if (file_exists($old_location . '/php_log.txt'))
	 {
		 copy($old_location . '/php_log.txt', $new_location . '/php_log.txt');
	 }
	 else{
		 $server_response['msg'] .= 'php_log not found ';
	 }
	error_log("[client: " . $_SERVER['REMOTE_ADDR'] . "] change_session.php [" . $_POST["old_session_id"] . "->" . $_POST["session_id"] . " ]> Success: " . ($server_response['success'] ? 'Yes' : 'No') . " | Message: " . $server_response['msg']);
	//Send info back to the client
	header('Content-type: application/json');
	echo json_encode($server_response);
 ?>