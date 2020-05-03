<?php
	// This php file contains get_labels function and is loaded at the beginning of upload_files. The function gets 4 arguments: (1 - data_file) the location of the 
	// tabular file, (2 - labelmatch_re) a pattern of headers containing labels (3 - labelmatch_re) a pattern for headers containing the peptide modifications and
	// (4 - labeldef_re) a pattern to discover all modifications applied to the peptides in this experiment
	
	function get_labels($data_file, $labelmatch_re, $labeldefcol_re, $labeldef_re) {
		$labels = [];
		$labels_defs = [];
		$labeldefcol = -1;
		$row = 0;
		if (($handle = fopen($data_file, "r")) !== FALSE)
		{
			if (($data = fgetcsv($handle, 0, "\t")) !== FALSE)
			{
				// data here stores the first line of the tabular file - the headers
				$row++;
				$ncols = count($data);
				// In label free data usually labelmatch_re and labeldef_re are set to null
				//
				// If we have labelled data (i.e. isobarically or metabolically labelled data)
				if (isset($labelmatch_re))
				{
					for ($c = 0; $c < $ncols; $c++)
					{
						// Loop all headers
						$matches = [];
						// For each header check if the labelmatch_re pattern matches the header:
						if (preg_match_all($labelmatch_re, $data[$c], $matches) > 0)
						{
							// Each header might contain up to 2 different labels (in form "[label1]/[label2]") so check for two matches and add the labels to the labels array
							// if they are not already there
							foreach ($matches[1] as $match)
							{
								if (!array_key_exists($match, $labels) && preg_match('/\+/', $match, $my_matches) == 0)
								{
									$labels[$match] = 1;
								}
							}
							if(isset($matches[2]))
							{
								foreach ($matches[2] as $match)
								{
									if (!array_key_exists($match, $labels) && preg_match('/\+/', $match, $my_matches) == 0) {
										$labels[$match] = 1;
									}
								}
							}
						}
						// If the header matches the label definition header regex assign labeldefcol to the current header index
						if (isset($labeldefcol_re) && $labeldefcol < 0 && preg_match($labeldefcol_re, $data[$c]))
						{
							$labeldefcol = $c;
						}
					}
				}
				else
				{
					// if we are looking only for the label definition header execute only the corresponding code block
					/*
						for ($c = 0; $c < $ncols; $c++)
						{
							$matches = [];
							if (preg_match($labeldefcol_re, $data[$c]))
							{
								$labeldefcol = $c;
								break;
							}
						}
					*/
				}
				
				// If we have labelled data (i.e. metabolically or isobarically labelled)
				if (isset($labelmatch_re))
				{
					// Search for all modifications in the file:
					// This feature is currently commented out in ProteoSign:
					/*
						while ($labeldefcol > -1 && ($data = fgetcsv($handle, 0, "\t")) !== FALSE)
						{
							$matches = [];
							$row++;
							if (count($data) >= $labeldefcol)
							{
								if (preg_match_all($labeldef_re, $data[$labeldefcol], $matches) > 0)
								{
									foreach ($matches[1] as $match)
									{
										# Label definitions must be more than one letter long (Ad hoc)
										if (strlen($match) > 1 && !array_key_exists($match, $labels_defs))
										{
											$labels_defs[$match] = 1;
										}
									}
								}
							}
						}
					*/
					//
					// An example of this block's return would be $labels_def ==
					/*
						Array
							(
								[Dimethyl] => 1
								[2H(4)] => 1
								[Oxidation] => 1
								[Carbamidomethyl] => 1
							)
					*/
				}
			}
			fclose($handle);
		}
		return [array_keys($labels), array_keys($labels_defs)];
	}
	
?>