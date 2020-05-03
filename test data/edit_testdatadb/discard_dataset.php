<?php
	// Used to discard a dataset and clean the db from everyu record pertaining to it
	
	$tables_with_dataset_id = array("dataset_files", "param_value", "experimental_structure", "extra_options");
	$tables_to_refresh_autoinc = array("dataset", "dataset_files", "param_value", "experimental_structure", "extra_options", "processed_files", "files");
	$path = "../testdatadb";
	$DatasetId = 0;
	try{
		$db = new SQLite3($path);
	}
	catch (Exception $exception) {
		echo "failed to open database";
		end;
	}
	
	
	function get_max_id($table)
	{
		global $db;
		$qres = $db->query('select * from ' . $table . ' where id = (select max(id) from ' . $table . ');');
		$res = $qres->fetchArray(SQLITE3_ASSOC);
		if ($res)
		{
			foreach($res as $col => $val){
				if ($col == "id")
				{
					return $val;
				}
			}
		}
		else
		{
			return 0;
		}
	}
	
	function del_where_dataset_id($table)
	{
		global $db;
		global $DatasetId;
		$db->query('delete from ' . $table . ' where dataset_id = ' . $DatasetId . ';');
	}
	
	function refresh_autoinc($table)
	{
		global $db;
		$db->query('update sqlite_sequence set seq = ' . get_max_id($table) . ' where name = "' . $table . '";');
	}
	
	echo "Please select the id of the dataset to erase\nBelow you will find all ids and their respective descriptions:\n";
	$qres = $db->query('select * from dataset');
	while($res = $qres->fetchArray(SQLITE3_ASSOC)){
		foreach($res as $col => $val){
			if ($col == "id")
			{
				echo $val;
			}
			else
			{
				echo "\t" . $val . "\n";
			}
		}
	}
	echo "Selected id: ";
	// Get the dataset description
	$DatasetId = trim(fgets(STDIN));
	$qres = $db->query('select * from dataset where id = ' . $DatasetId);
	$res = $qres->fetchArray(SQLITE3_ASSOC);
	if ($res)
	{
		foreach($res as $col => $val){
			if ($col == "desc")
			{
				echo "\nDelete " . $val . "? (y/n)";
				$my_answer = trim(fgets(STDIN));
				if ($my_answer == "n")
				{
					die;
				}
				elseif ($my_answer == "y")
				{
					// Delete the dataset with the id $DatasetId and erase all relative entries in all tables, then refresh the autoincrement counters for all tables
					// First query the dataset_files table to find all associated with the current dataset files and delete them from the files table
					$qres = $db->query('select * from dataset_files where dataset_id = ' . $DatasetId . ';');
					while($res = $qres->fetchArray(SQLITE3_ASSOC)){
						foreach($res as $col => $val){
							if ($col = "id")
							{
								$db->query('delete from files where id = ' . $val . ";");
							}
						}
					}
					// Do the same thing with the processed_files table - get the association of the dataset from experimental_structure and erase the entries from processed_files
					$qres = $db->query('select processed_file_id from experimental_structure where dataset_id = ' . $DatasetId . ';');
					while($res = $qres->fetchArray(SQLITE3_ASSOC)){
						foreach($res as $col => $val){
							if ($col = "processed_file_id")
							{
								$db->query('delete from processed_files where id = ' . $val . ";");
							}
						}
					}
					// Erase all entries with the current dataset_id in all tables where this is possible
					foreach($tables_with_dataset_id as $my_table)
					{
						del_where_dataset_id($my_table);
					}
					// Delete the entry from the dataset table
					$db->query('delete from dataset where id = ' . $DatasetId . ';');
					// Also refresh the autoincrement counter for all the tables
					foreach($tables_to_refresh_autoinc as $my_table)
					{
						refresh_autoinc($my_table);
					}
				}
				else
				{
					echo "Invalid answer";
					die;
				}
			}
		}
	}
	else
	{
		echo "Invalid ID";
		die;
	}
	echo "\nDataset " . $DatasetId . " was discarded successfully!";
	
?>	