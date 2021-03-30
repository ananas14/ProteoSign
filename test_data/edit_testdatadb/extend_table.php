<?php
	$path = "../testdatadb";
	try{
		$db = new SQLite3($path);
	}
	catch (Exception $exception) {
		echo "failed to open database";
		end;
	}
	echo "Table to extend: ";
	$table = trim(fgets(STDIN));
	echo "\nFile to get rows from: ";
	$FileName = trim(fgets(STDIN));
	if ($FileName == "") $FileName = "dat.myxls";
	echo "\n";
	$myfile = fopen($FileName, "r") or die("Unable to open file!");
	if ($table == "dataset")
	{
		$columnnames = "id, desc";
	}
	elseif($table == "param_value")
	{
		$columnnames = "id, param_id, value, dataset_id";
	}	
	elseif($table == "processed_files")
	{
		$columnnames = "id, name";
	}	
	elseif($table == "experimental_structure")
	{
		$columnnames = "id, processed_file_id, brep, trep, frac, dataset_id, used, condition";
	}	
	elseif($table == "dataset_files")
	{
		$columnnames = "id, dataset_id, file_id";
	}	
	elseif($table == "files")
	{
		$columnnames = "id, file";
	}
	elseif($table == "extra_options")
	{
		$columnnames = "id, dataset_id, option, opt_value";
	}
	else
	{
		echo "there is no such table";
		die;
	}
	
	while(!feof($myfile)) {
		$myline = fgets($myfile);
		if ($myline == "")
		{
			echo "Skipped empty line!";
			continue;
		}
		$my_array = explode("\t", $myline);
		$my_sql_vals = "";
		foreach($my_array as &$my_val)
		{
			$my_val = trim($my_val);
			$my_sql_vals = $my_sql_vals . "'" . $my_val . "',";
		}
		$my_sql_vals = rtrim($my_sql_vals, ",");
		$thequery = 'insert into ' . $table . " (" . $columnnames . ") values (" . $my_sql_vals . ");";
		$db->query('insert into ' . $table . " (" . $columnnames . ") values (" . $my_sql_vals . ");");
		unset($my_val);
		
	}
	fclose($myfile);
?>