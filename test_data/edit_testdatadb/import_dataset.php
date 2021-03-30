<?php
	
	$path = "../testdatadb";
	try{
		$db = new SQLite3($path);
	}
	catch (Exception $exception) {
		echo "failed to open database";
		end;
	}
	
	$dataset_id = 1;
	$first_processed_file = 1;
	$isLabelFree = false;
	
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
	
	function extend_table($table, $tab_string)
	{
		global $db, $dataset_id;
		// This function extends a table with the values of a tabular string
		if ($table == "dataset")
		{
			$columnnames = "desc";
		}
		elseif($table == "param_value")
		{
			$columnnames = "param_id, value, dataset_id";
		}	
		elseif($table == "processed_files")
		{
			$columnnames = "name";
		}	
		elseif($table == "experimental_structure")
		{
			$columnnames = "processed_file_id, brep, trep, frac, dataset_id, used, condition";
		}	
		elseif($table == "dataset_files")
		{
			$columnnames = "dataset_id, file_id";
		}	
		elseif($table == "files")
		{
			$columnnames = "file";
		}
		elseif($table == "extra_options")
		{
			$columnnames = "dataset_id, option, opt_value";
		}
		else
		{
			return false;
		}
		$my_lines = explode("\n", $tab_string);
		foreach($my_lines as &$my_line)
		{
			if ($my_line == "")
			{
				continue;
			}
			$my_values = explode("\t", $my_line);
			$my_sql_vals = "";
			foreach($my_values as &$value)
			{
				$value = trim($value);
				$my_sql_vals = $my_sql_vals . "'" . $value . "',";
			}
			$my_sql_vals = rtrim($my_sql_vals, ",");
			$db->query('insert into ' . $table . " (" . $columnnames . ") values (" . $my_sql_vals . ");");
		}
	}
	// Get the dataset description
	echo "\nName of the new test dataset: ";
	$DatasetName = trim(fgets(STDIN));
	
	
	// Open the parameters file
	echo "\nParameters file: ";
	$FileName = trim(fgets(STDIN));
	if ($FileName == "") $FileName = "dat.myxls";
	echo "\n";
	$myparamfile = fopen($FileName, "r") or die("Unable to open file!");
	
	// Set the dataset name and retrieve its id
	
	$db->query('insert into dataset (desc) values ("' . $DatasetName . '");');
	$dataset_id = get_max_id("dataset");
	
	// Ask for filenames and associate them with the new dataset
	echo "\nAdd a new data file (type nothing and hit enter to stop adding files): ";
	
	$DataFileName = trim(fgets(STDIN));
	while($DataFileName != "")
	{
		extend_table("files", $DataFileName);
		$temp_id = get_max_id("files");
		$temp_tab_string = $dataset_id . "\t" . $temp_id;
		extend_table("dataset_files", $temp_tab_string);
		echo "\nAdd a new data file (type nothing and hit enter to stop adding files): ";
		$DataFileName = trim(fgets(STDIN));
	}
	
	$myLFQconds = array();
	// First skim through the parameter file and search for LFQconditions
	while(!feof($myparamfile))
	{
		$myline = fgets($myparamfile);
		$myline = trim($myline);
		if ($myline == "LFQconditions")
		{
			$myval = fgets($myparamfile);
			$myval = trim($myval);
			if ($myval != "")
			{
				$my_tab_vals = explode("\t\t", $myval);
				foreach($my_tab_vals as &$tab_val)
				{
					$tab_val = trim($tab_val);
					if ($tab_val != "")
					{
						$mypieces = explode("\t", $tab_val);
						$myLFQconds[trim($mypieces[0])] = trim($mypieces[1]);
					}
				}
			}
		}
	}
	fclose($myparamfile);
	// Re-open the file
	$myparamfile = fopen($FileName, "r") or die("Unable to open file!");
	// Get all lines oif the parameters file and use the necessary information
	while(!feof($myparamfile))
	{
		$myline = fgets($myparamfile);
		$myline = trim($myline);
		
		
		
		if ($myline == "isLabelFree")
		{
			$myval = fgets($myparamfile);
			$myval = trim($myval);
			if ($myval == "false")
			{
				$isLabelFree = false;
			}
			else
			{
				$isLabelFree = true;
			}
		}
		if ($myline == "rawfiles_structure")
		{
			$first_processed_file = get_max_id("processed_files") + 1;
			$mytabstring = fgets($myparamfile);
			$mytabstring = str_replace("\t\t", "\n", $mytabstring);
			$mytabstring_structure = $mytabstring;
			$mytabstring = preg_replace('/\t.*?\n/', "\n", $mytabstring);
			extend_table("processed_files", $mytabstring);
			$mytabstring_structure = preg_replace('/^.*?\t/', "", $mytabstring_structure, 1);
			$mytabstring_structure = preg_replace('/(\n).*?\t/', "$1", $mytabstring_structure);
			
			// $mytabstring_structure contains basic info for the structure. If its not LFQ the structure per line is "brep trep frac used"
			$my_tab_lines = explode("\n", $mytabstring_structure);
			$my_final_tabstring = "";
			$i = 0;
			foreach($my_tab_lines as &$tab_line)
			{
				$tab_line = trim($tab_line);
				if ($tab_line == "")
				{
					continue;
				}
				$my_tab_tabs = explode("\t", $tab_line);
				$my_final_tabstring .= ($first_processed_file + $i) . "\t";
				if (!$isLabelFree)
				{
					$j = 0;
					foreach($my_tab_tabs as &$tab_tab)
					{
						$j++;
						if ($j < 4)
						{
							// add brep trep and frac
							$tab_tab = trim($tab_tab);
							$my_final_tabstring .= $tab_tab . "\t";
						}
						else
						{
							// add dataset id, used and condition
							$tab_tab = trim($tab_tab);
							$my_final_tabstring .= $dataset_id . "\t";
							$my_final_tabstring .= $tab_tab . "\t";
							$my_final_tabstring .= '-' . "\n";
						}
					}
					$i++;
				}
				else
				{
					$j = 0;
					foreach($my_tab_tabs as &$tab_tab)
					{
						$j++;
						if ($j < 4)
						{
							// add brep trep and frac
							$tab_tab = trim($tab_tab);
							$my_final_tabstring .= $tab_tab . "\t";
						}
						elseif ($j = 4)
						{
							// add dataset id, used
							$tab_tab = trim($tab_tab);
							$my_final_tabstring .= $dataset_id . "\t";
							$my_final_tabstring .= $tab_tab . "\t";
							// Finally add the condition - get the current processed file and get the respective condition from the associative array built above
							$qres = $db->query('select name from processed_files where id = ' . ($first_processed_file + $i));
							if(!$qres){
								echo "failed to query the db...\n";	
							}
							$curcondition = "";
							while($qrow = $qres->fetchArray(SQLITE3_ASSOC)){
								foreach($qrow as $col => $val){
									// The name of the processed_file we are looking for is under $val
									if (array_key_exists($val, $myLFQconds))
									{
										$curcondition = $myLFQconds[$val];
									}
									else
									{
										$curcondition = '-';
									}
								}
							}
							$my_final_tabstring .= $curcondition . "\n";
						}
					}
					$i++;
				}
				
			}
			extend_table("experimental_structure", $my_final_tabstring);
		}
		elseif ($myline == "procprogram")
		{
			$myproc = trim(fgets($myparamfile));
			if ($myproc == "PD")
			{
				$mytabstring = "1" . "\t" . "1" . "\t" . $dataset_id;
			}
			else
			{
				$mytabstring = "1" . "\t" . "0" . "\t" . $dataset_id;
			}
			extend_table("param_value", $mytabstring);
		}
		elseif ($myline == "exptpoint")
		{
			$myexptpoint = trim(fgets($myparamfile));
			$mytabstring = "2" . "\t" . $myexptpoint . "\t" . $dataset_id;
			extend_table("param_value", $mytabstring);
		}		
		elseif ($myline == "quantitation_filtering")
		{
			$myopt = trim(fgets($myparamfile));
			if ($myopt == "F")
			{
				$mytabstring = "9" . "\t" . "0" . "\t" . $dataset_id;
			}
			else
			{
				$mytabstring = "9" . "\t" . "1" . "\t" . $dataset_id;
			}
			extend_table("param_value", $mytabstring);
		}
		elseif ($myline == "filtering_label")
		{
			$myopt = trim(fgets($myparamfile));
			$mytabstring = "10" . "\t" . $myopt . "\t" . $dataset_id;
			extend_table("param_value", $mytabstring);
		}
		elseif ($myline == "peptide_level_filtering")
		{
			$myopt = trim(fgets($myparamfile));
			if ($myopt == "F")
			{
				$mytabstring = "11" . "\t" . "0" . "\t" . $dataset_id;
			}
			else
			{
				$mytabstring = "11" . "\t" . "1" . "\t" . $dataset_id;
			}
			extend_table("param_value", $mytabstring);
		}
		elseif ($myline == "!Rename")
		{
			$myopt = trim(fgets($myparamfile));
			$myopt = rtrim($myopt, "|");
			$mytabstring = $dataset_id . "\t" . "Rename" . "\t" . $myopt;
			extend_table("extra_options", $mytabstring);
		}
		elseif ($myline == "conditions_to_compare")
		{
			$myopt = trim(fgets($myparamfile));
			$myopt = str_replace("\t", "|", $myopt);
			$myopt = rtrim($myopt, "|");
			$mytabstring = $dataset_id . "\t" . "Select_Labels" . "\t" . $myopt;
			extend_table("extra_options", $mytabstring);
		}
		elseif ($myline == "!LS_Array")
		{
			$myopt = trim(fgets($myparamfile));
			$myopt = rtrim($myopt, "|");
			if ($myopt == "")
			{
				continue;
			}
			$mytabstring = $dataset_id . "\t" . "LS_Array" . "\t" . $myopt;
			extend_table("extra_options", $mytabstring);
		}
		elseif ($myline == "!LS_c_p_Add")
		{
			$myopt = trim(fgets($myparamfile));
			$myopt = rtrim($myopt, "|");
			if ($myopt == "")
			{
				continue;
			}
			$my_parts = explode("||", $myopt);
			foreach($my_parts as &$part)
			{
				if (preg_match("/\d\|/", $part))
				{
					$my_parts2 = explode("|", $part);
					$part = $my_parts2[0] . "|" . $my_parts2[count($my_parts2) - 1];
				}
			}
			$myopt = implode("||", $my_parts);
			$myopt = str_replace("||", "|", $myopt);
			$mytabstring = $dataset_id . "\t" . "LS_c_p_Add" . "\t" . $myopt;
			extend_table("extra_options", $mytabstring);
		}
	}
	fclose($myparamfile);
	echo $DatasetName . " was imported successfully!\n";
?>		