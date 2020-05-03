 <?php
	// number of days the website source (html, php, js and css files) was changed
	$softchanges = 153;
	$ver = [1,0,0,0];
	$server_response = [];
	/*
	// Calculation of "version string" (x.x.x.x) based on the number of file edits (not creations) since the beginning of the project
	// Examples: 1234 translates to version 2.2.3.4, 4109 translates to version 5.1.0.9, 12078 translates to version 13.0.7.8
	$server_response['version'] = '';
	$i=0;
	while($i<count($ver) && $softchanges > 0){
		if(pow(10,(count($ver)-$i-1)) <= $softchanges){
			$intlog = (int)(log10($softchanges));
			if($intlog > count($ver)-1){
				$intlog = count($ver)-1;
			}
			$digit = (int)($softchanges/(pow(10,$intlog)));
			$ver[$i] = $ver[$i] + $digit;
			$softchanges -= $digit*pow(10,$intlog);
		}
		$server_response['version'] = $server_response['version'] . '.' . $ver[$i];
		$i++;
	}
	
	$server_response['version'] = substr($server_response['version'],1);
	*/
	$server_response['version'] = '1.0';
	$server_response['hostaddr'] = $_SERVER['REMOTE_ADDR'];
	$server_response['hostname'] = gethostbyaddr($server_response['hostaddr']);
	
	header('Content-type: application/json');
	echo json_encode($server_response);		
 ?>