<?php
	
	$path = "../testdatadb";
	try{
		$db = new SQLite3($path);
	}
	catch (Exception $exception) {
		echo "failed to open database";
		end;
	}
	$option = "s";
	if ($option == "s")
	{
		echo "Datasets (dataset):\r\n";
		$qres = $db->query('select * from dataset');
		if (!$qres)
		{
			echo "(Empty table)";
		}
		else
		{
			$qrow = $qres->fetchArray(SQLITE3_ASSOC);
			foreach($qrow as $col => $val){
				echo $col . "\t";
			}
			echo "\r\n";
			foreach($qrow as $col => $val){
				echo $val . "\t";
			}
			echo "\r\n";
			while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
				foreach($qrow as $col => $val){
					echo $val . "\t";
				}
				echo "\r\n";
			}	
		}
		
		
		
		echo "\nDataset files (dataset_files):\r\n";
		$qres = $db->query('select * from dataset_files');
		if (!$qres)
		{
			echo "(Empty table)";
		}
		else
		{
			$qrow = $qres->fetchArray(SQLITE3_ASSOC);
			foreach($qrow as $col => $val){
				echo $col . "\t";
			}
			echo "\r\n";
			foreach($qrow as $col => $val){
				echo $val . "\t";
			}
			echo "\r\n";
			while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
				foreach($qrow as $col => $val){
					echo $val . "\t";
				}
				echo "\r\n";
			}		
		}
		
		
		
		echo "\nFiles (files):\r\n";
		$qres = $db->query('select * from files');
		if (!$qres)
		{
			echo "(Empty table)";
		}
		else
		{
			$qrow = $qres->fetchArray(SQLITE3_ASSOC);
			foreach($qrow as $col => $val){
				echo $col . "\t";
			}
			echo "\r\n";
			foreach($qrow as $col => $val){
				echo $val . "\t";
			}
			echo "\r\n";
			while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
				foreach($qrow as $col => $val){
					echo $val . "\t";
				}
				echo "\r\n";
			}	
		}
		
		
		
		echo "\nParameters (param_value):\r\n";
		$qres = $db->query('select * from param_value');
		if (!$qres)
		{
			echo "(Empty table)";
		}
		else
		{
			$qrow = $qres->fetchArray(SQLITE3_ASSOC);
			foreach($qrow as $col => $val){
				echo $col . "\t";
			}
			echo "\r\n";
			foreach($qrow as $col => $val){
				echo $val . "\t";
			}
			echo "\r\n";
			while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
				foreach($qrow as $col => $val){
					echo $val . "\t";
				}
				echo "\r\n";
			}		
		}
		
		
		
		echo "\nProcessed files (processed_files):\r\n";
		$qres = $db->query('select * from processed_files');
		if (!$qres)
		{
			echo "(Empty table)";
		}
		else
		{
			$qrow = $qres->fetchArray(SQLITE3_ASSOC);
			foreach($qrow as $col => $val){
				echo $col . "\t";
			}
			echo "\r\n";
			foreach($qrow as $col => $val){
				echo $val . "\t";
			}
			echo "\r\n";
			while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
				foreach($qrow as $col => $val){
					echo $val . "\t";
				}
				echo "\r\n";
			}	
		}
		
		
		echo "\nExperimental structure (experimental_structure):\r\n";
		$qres = $db->query('select * from experimental_structure');
		if (!$qres)
		{
			echo "(Empty table)";
		}
		else
		{
			$qrow = $qres->fetchArray(SQLITE3_ASSOC);
			foreach($qrow as $col => $val){
				echo $col . "\t";
			}
			echo "\r\n";
			foreach($qrow as $col => $val){
				echo $val . "\t";
			}
			echo "\r\n";
			while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
				foreach($qrow as $col => $val){
					echo $val . "\t";
				}
				echo "\r\n";
			}	
		}
		
		
		echo "\nParameters (param):\r\n";
		$qres = $db->query('select * from param');
		if (!$qres)
		{
			echo "(Empty table)";
		}
		else
		{
			$qrow = $qres->fetchArray(SQLITE3_ASSOC);
			foreach($qrow as $col => $val){
				echo $col . "\t";
			}
			echo "\r\n";
			foreach($qrow as $col => $val){
				echo $val . "\t";
			}
			echo "\r\n";
			while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
				foreach($qrow as $col => $val){
					echo $val . "\t";
				}
				echo "\r\n";
			}	
		}
		
		echo "\nExtra options (extra_options):\r\n";
		$qres = $db->query('select * from extra_options');
		if (!$qres)
		{
			echo "(Empty table)";
		}
		else
		{
			$qrow = $qres->fetchArray(SQLITE3_ASSOC);
			foreach($qrow as $col => $val){
				echo $col . "\t";
			}
			echo "\r\n";
			foreach($qrow as $col => $val){
				echo $val . "\t";
			}
			echo "\r\n";
			while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
				foreach($qrow as $col => $val){
					echo $val . "\t";
				}
				echo "\r\n";
			}	
		}
		
	}
	// else
	// {
	// echo "Datasets (dataset):\r\n";
	// $qres = $db->query('select * from dataset');
	// while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
	// foreach($qrow as $col => $val){
	// echo $col . ": " . $val . "\r\n";
	// }
	// echo "\r\n";
	// }			
	
	
	
	// echo "Dataset files (dataset_files):\r\n";
	// $qres = $db->query('select * from dataset_files');
	// while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
	// foreach($qrow as $col => $val){
	// echo $col . ": " . $val . "\r\n";
	// }
	// echo "\r\n";
	// }			
	
	
	// echo "Files (files):\r\n";
	// $qres = $db->query('select * from files');
	// while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
	// foreach($qrow as $col => $val){
	// echo $col . ": " . $val . "\r\n";
	// }
	// echo "\r\n";
	// }	
	
	
	// echo "Parameters (param_value)\r\n";
	// $qres = $db->query('select * from param_value');
	// while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
	// foreach($qrow as $col => $val){
	// echo $col . ": " . $val . "\r\n";
	// }
	// echo "\r\n";
	// }	
	
	
	// echo "Processed files (processed_files):\r\n";
	// $qres = $db->query('select * from processed_files');
	// while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
	// foreach($qrow as $col => $val){
	// echo $col . ": " . $val . "\r\n";
	// }
	// echo "\r\n";
	// }		
	
	
	// echo "Experimental structure (experimental_structure):\r\n";
	// $qres = $db->query('select * from experimental_structure');
	// while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
	// foreach($qrow as $col => $val){
	// echo $col . ": " . $val . "\r\n";
	// }
	// echo "\r\n";
	// }
	
	// echo "Parameters (param):\r\n";
	// $qres = $db->query('select * from param');
	// while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
	// foreach($qrow as $col => $val){
	// echo $col . ": " . $val . "\r\n";
	// }
	// echo "\r\n";
	// }			
	
	// echo "Extra Options (extra_options):\r\n";
	// $qres = $db->query('select * from extra_options');
	// while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
	// foreach($qrow as $col => $val){
	// echo $col . ": " . $val . "\r\n";
	// }
	// echo "\r\n";
	// }
	
	// }
?>