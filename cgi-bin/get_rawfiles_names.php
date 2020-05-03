<?php
	function get_rawfiles_names($data_file, $labeldefcol_re){
		$rawfiles_names = [];
		$ret = [];
		$labeldefcol = -1;
		if (($handle = fopen($data_file, "r")) !== FALSE) {
			if (($data = fgetcsv($handle, 0, "\t")) !== FALSE) {
				$ncols = count($data);
				// find the column which contains the raw data file names
				for ($c=0; $c < $ncols; $c++) {
					$matches = [];
					if(preg_match($labeldefcol_re, $data[$c])){
						$labeldefcol = $c;
						break;
					}
				}
				if($labeldefcol > -1){
					while (($data = fgetcsv($handle, 0, "\t")) !== FALSE) {
						if (count($data) >= $labeldefcol)
						{
							if(! isset($rawfiles_names[$data[$labeldefcol]])){
								$rawfiles_names[$data[$labeldefcol]] = 1;
							}
						}
					}
					$ret = array_keys($rawfiles_names);
				}
			}
			fclose($handle);
		}
		return $ret;
	}
?>