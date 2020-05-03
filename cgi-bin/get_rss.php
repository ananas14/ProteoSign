 <?php
	$server_response = [];
	$url = $_POST["rssurl"];
	$xml = simplexml_load_file($_POST["rssurl"]);
	$i = 0;
	if(!empty($xml))
	{
		foreach ($xml->channel->item as $item)
		{
			$server_response[trim((string)$item->title)] = trim((string)$item->link);
			$i++;
			if ($i > 40) break;
		}
	}
	header('Content-type: application/json');
	echo json_encode($server_response);
	
 ?>