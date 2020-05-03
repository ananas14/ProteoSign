<?php
	$path = "../testdatadb";
	try{
		$db = new SQLite3($path);
	}
		catch (Exception $exception) {
		echo "failed to open database";
		end;
	}
	echo "enter your query:\n";
	$my_q = trim(fgets(STDIN));
	beginning:
	echo "execute the following query? (y/n)";
	echo $my_q;
	echo "\n";
	$my_answer = trim(fgets(STDIN));
	if ($my_answer == "n")
	{
		die;
	}
	elseif ($my_answer == "y")
	{
		$qres = $db->query($my_q);
	}
	else
	{
		echo "Please try again";
		goto beginning;
	}
?>