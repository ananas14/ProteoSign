 <?php
	$server_response = [];
	$opts = array('http'=>array('header' => "User-Agent:MyAgent/1.0\r\n"));
	$context = stream_context_create($opts);
	$html = file_get_contents($_POST['rssurl'],false,$context);
	preg_match('/<div id="first">([\s\S]+?)<\/div>/i', $html, $matches);
	$myresult = $matches[1];
	$server_response['html_code'] = $myresult;
	header('Content-type: application/json');
	echo json_encode($server_response);
 ?>