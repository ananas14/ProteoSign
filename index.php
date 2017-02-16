<!DOCTYPE html>

<?php
	//To prevent cache from storing the website: (from http://docstore.mik.ua/orelly/webprog/php/ch07_05.htm)
	header("Cache-Control: no-cache, no-store, must-revalidate"); // HTTP 1.1.
	header("Pragma: no-cache"); // HTTP 1.0.
	header("Expires: 0"); // Proxies.
?>

<html>
	<head>
		<title>ProteoSign</title>
		<link rel="stylesheet" type="text/css" href="css/style.css">
		<link rel="stylesheet" type="text/css" href="css/jquery-ui.min.css">
		<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.5.0/css/font-awesome.min.css">
		<meta name="keywords" content="ProteoSign, mass spectrometry, proteomics, data, statistics, statistical analysis, bioinformatics, limma, comparative proteomics, silac, labeling, differential expression">
		<meta charset="UTF-8">
		<link rel="shortcut icon" href="favicon.ico" type="image/x-icon">
	</head>
	<body>
		<div class="main_div">
			<div class="main_div main_header">
				<span class="logo">ProteoSign</span><br>
				<span class="logo" style="font-size: 100%; left: 12%">comparative proteomics</span>
				<a href="help.html" target="_blank"><i class="fa fa-question-circle helplink" aria-hidden="true"></i></a><a class="helplink" href="help.php" target="_blank"> HELP</a>
			</div>
			<div class="main_div main_nav">
				<p id="scrollingtext" class="scrolling_p"></p>
			</div>
			<div class="main_div main_section" id="s1div">
				<p><strong>ProteoSign</strong> is an online service for protein <strong>differential expression (or abundance) analysis</strong> designed with the end-proteomics user in mind.</p>
				<p>Borrowing from the experience of the microarray field, ProteoSign software utilizes the well-established <a href="http://bioinf.wehi.edu.au/limma/" target="_blank"><strong>Linear Models For Microarray Data (LIMMA)</strong></a> methodology for assessing statistical significance in protein abundance differences between two or more proteome states.</p>
				<p>ProteoSign aims to fully automate the process of statistically evaluating differential expression in mass spectrometry-based bottom-up (shotgun) quantitative proteomics data by requiring <strong>minimal user interaction time</strong> and generating <strong>publication-quality data plots.</strong></p>
				<br><hr class="style-one"><br>
				<p style="color: #000000">Click <em>Start</em> to begin your analysis.</p>
				<br>
				<div class="buttonbar">
					<button type="button" class="main" id="s1btnf">Start</button>
				</div>
			</div>
			<div class="main_div main_section hidden" id="s2div">
				<h2>Step 1</h2>
				<p style="color: #000000">Please read the text below. When you are ready, upload your dataset and click <em>Next</em>. Alternatively, you can <a href="#" onclick="onChooseFromTestDatasets();"><u>choose</u></a> from available test datasets.</p>
				<p>
				<dl style="font-family: 'Courier New'; font-size: 75%; padding-left: 5%;">
					<dt>Input file support:</dt>
					<dd>&#8212 <a href="http://www.google.com/search?q=Thermo%20Proteome%20Discoverer&btnI" target="_blank">Proteome Discoverer&#8482</a> (PD) v1.4+ PSM file.<span class="tooltip" style="position: relative">?<div style="left: -25em; min-width: 50em;"><img class="callout" src="images/callout_black.gif" style="margin-left: 24.3em"/>The <em>peptide spectrum match</em> (PSM) file, comprising the necessary information required for protein quantitation, has to be produced manually through the PD desktop application (see the <em>Export&#x2192;To Text</em> option under the <em>File</em> menu).<p>The PSM file must be created from within a PD report which combines the information from all relevant MS analyses (i.e. replicates, conditions etc). Such a PD report is sometimes required to be created based on separate reports and is referred to as a multi-consensus report.</p><p>Before exporting the information, please <strong>disable</strong> peptide grouping and edit the report's <em>Quantification&nbsp;Method</em> by configuring the <em>Ratio&nbsp;Calculation</em> parameters as shown below. The following screen-shots show where the aforementioned changes can be made from within the PD application.</p><figure style="float: left; width: 49%; margin: 0 auto;"><img src="images/pdpepgroup.png" style="width: 100%;"/><figcaption style="text-align: center;">Disabling peptide grouping.</figcaption></figure><figure style="float: right; width: 49%; margin: 0 auto;"><img src="images/pdquantmeth.png" style="width: 100%;"/><figcaption style="text-align: center;">Editing ratio calculation parameters.</figcaption></figure></div></span></dd>
					<dd>&#8212 <a href="http://www.google.com/search?q=Max%Planck%MaxQuant&btnI" target="_blank">MaxQuant</a> (MQ) v1.3.0.5+ 'proteinGroups' and 'evidence' files.</dd>
					<dt><br>Experiment support:</dt>
					<dd>&#8212 Label-free, <a href="http://www.ncbi.nlm.nih.gov/pubmed/?term=15385600" target="_blank">iTRAQ</a>, <a href="http://www.ncbi.nlm.nih.gov/pubmed/?term=12713048" target="_blank">TMT</a>, <a href="http://www.ncbi.nlm.nih.gov/pubmed/?term=12118079" target="_blank">SILAC</a>, <a href="http://www.ncbi.nlm.nih.gov/pubmed/?term=19053139" target="_blank">pulsed SILAC</a>, <a href="http://www.ncbi.nlm.nih.gov/pubmed/?term=19300442" target="_blank">dimethyl</a>.</dd>
				</dl>
				</p>
				 <p style="color: #000000">The input file(s) should be made out of merging multiple LC-MS/MS runs (fractionation, replication) of a <em>single</em>, possibly multiplexed, experiment and should not be the result of merging different or parallel experiments. Thus, <em>multi-experiment analysis is not supported.</em> In addition, only single factor (one-way) analysis is currently supported.</p>
				<!-- Actual file chooser button, hidden (it can't be avoided) -->
				<input type="file" id="__s2btnupld" style="display: none" multiple/>
				<p>
					<button class="main upload" type="file" id="s2btnupld"/>Upload file(s)</button>
				</p>
				<!-- File uploaders will be rendered below when user selects file(s) -->
				<div class="sectionForm" id="s2uluploaders" style="height: 20%">
					<table style="border-spacing: 15px 5px;">
					</table>
				</div>
				<div class="buttonbar">
					<button type="button" class="main" id="s2btnb">Back</button>
					<button type="button" class="main" id="s2btnf" disabled>Next</button>
				</div>
			</div>
			
			<div class="main_div main_section hidden" id="s22div">
				<h2>Step 2</h2>
				<p style="color: #000000">Please define the experimental design by assigning the appropriate structure coordinate to each file. 
					<p id = "step2desc" style="font-size: 14px;">A coordinate specifies the biological replicate, technical replicate and fraction a raw file belongs to. Select one or more files from the table on the left, type the biological and technical replicate numbers they belong to in the boxes on the right, and click the <strong>Assign</strong> button. If multiple files are fractions of the same replicate, simply assign the same biological and technical replicate number to them, the order of the fractions does not matter. Use <strong>Shift</strong> to select a range of items. <strong>Right click</strong> on the table for extended functionality.</p>
					<p style="font-size: 14px;" id="ExtraInfoForLFQ">Please use the <strong>"Assign Condition"</strong> option in the right click menu to tell ProteoSign which Raw Files correspond to which condition.</p>
				</p>
				<hr class="style-one">
            <div style="text-align: center; height: 45%; width: 100%; display: -webkit-flex; display: inline-flex;">
			<div class="expstructcoord_tbl" style="margin: 10px; height: 100%; overflow: auto; max-height: 100%; width: 70%;">
                  <!-- <div style="margin-top: 1em;"> -->
						<div style="overflow:auto; height: 200px; margin-top: 1em;">
                     <table class="compact row-border rawfiles_tbl" cellspacing="0" width="100%" id="rawfiles_tbl_allfiles">
					 <!--<thead>
									<tr>
										<th><span style="position:absolute; left: 0;top: 0; border-bottom: 1px solid #aaa;">Raw files</span></th>
									</tr>
								</thead>-->
								<tbody>
								</tbody>
							</table>
						</div>
					<!-- </div> -->
					</div>
					<div id="expstructcoord" style="margin: 10px; margin-top: 30px; width: 30%; position: relative;" class="expstructcoord_tbl">
                  <div style="margin-top: 1em">Experimental structure coordinate</div>
                  <div style="position: absolute; top: 40%; left: 0; right: 0;">
                     <input id="expstructcoord_biorep" type="text" placeholder="Bio. replicate #" onkeypress="inputChCheck(event, '^[0-9]+$', 2)">
                     <input id="expstructcoord_techrep" type="text" placeholder="Tech. replicate #" onkeypress="inputChCheck(event, '^[0-9]+$', 2)">
                  </div>
                  <div style="position: absolute; bottom: 1em; left: 0; right: 0;">
                     <button type="button" class="main inline" style="font-size: 80%; padding: 5px 15px 6px 15px;" id="btnAssignExpStructCoord">Assign</button>
                     <button type="button" class="main inline" style="font-size: 80%; padding: 5px 15px 6px 15px;" id="btnResetExpStructCoord" disabled="">Reset</button>
                  </div>
               </div>
				</div>
				<div class="buttonbar">
					<button type="button" class="main" id="s22btnb">Back</button>
					<button type="button" class="main" id="s22btnf" disabled onclick="ons22btnfclick();">Next</button>
				</div>
			</div>
			
			<div class="main_div main_section hidden" id="s3div">
				<h2>Step 3</h2>
				<p style="color: #000000">Define the experimental parameters and click <em>Submit</em>.</p>
				<div class="sectionForm" id="s3expparams">
					<table name = "tbl_main_params" id = "tbl_main_params" style="border-spacing: 5px 5px;">
						<tr>
							<td id="quantitation_prog_lbl">&#8212 Raw data were quantified with Proteome Discoverer&#8482</td>
							<td><input name="exppddata" type="checkbox"  onclick="return false" style="visibility: hidden"/></td>
							<td></td>
							<td></td>
							<td></td>
						</tr>
						<tr>
							<td>&#8212 Experiment ID</td>
							<td></td>
							<td><input data-required="true" name="expid" type="text" onkeypress="inputChCheck(event,'^(?!_)[a-zA-Z0-9_]+$',25)" placeholder="Character rules apply"/></td>
							<td></td>
							<td></td>
							<td></td>
						</tr>
						<tr>
							<td>&#8212 Experiment description </td>
							<td></td>
							<td><input data-required="true" name="exptpoint" type="text" onkeypress="inputChCheck(event,'^(?!_)[a-zA-Z0-9_]+$',20)" placeholder="Character rules apply"/></td>
							<td></td>
							<td></td>
							<td></td>
						</tr>
						<tr>
							<td>&#8212 Conditions to compare </td>
							<td></td>
							<td id="conds_list_container">
								<select name = "conditions_list" class="conditions_list" onclick="onlistclick()" id = "conditions_list" style="width: 100%" multiple oncontextmenu="return false;" onmouseup="oncondslistclick()">
								
									
								</ul>
							</td>
						</tr>
					</table>
					<div id="adv_parameters_div"> <!-- WARNING!! in case we add another advanced parameter other than Quantitation filtering change script.js function addformlabel so that if label-free data are detected, only Quantitation filtering will be hidden and not advanced parameters as a whole -->
					<u class="astyle" id="s3showhideadvparams">Show advanced parameters</u>
					<table id="s3advparams" class="hidden" style="border-spacing: 5px 5px;">
						<tr class="hidden">
							<td>&#8212 Protein (as opposed to peptide) quantitation?</td>
							<td></td>
							<td><input name="expprotquant" type="checkbox" checked/></td>
						</tr>
 						<tr class='hidden'>
							<td>&#8212 Non-labelled "background" species present?</td>
							<td></td>
							<td><input name="explbl00" type="checkbox"/></td>
						</tr>
						<tr class='hidden'>
							<td>&#8212 "Background" species</td>
							<td></td>
							<td><select id="explbl0name_" name="explbl0"></select></td>
						</tr>
						<tr id="s3QuantitationFiltering">
							<td>&#8212 Quantitation filtering?</td>
							<td><span class="tooltip">?<span><img class="callout" src="images/callout_black.gif" /><strong>Warning!</strong><br><u>For labelled experiments only</u>. If enabled, two options will become available: <strong>a)</strong> singlets (peptides with no signal at the MS level in all but one labelled version) can be excluded from protein quantitation (this is referred to as <u>peptide-level filtering</u>).<br><strong>b)</strong> Proteins identified just by peptides with a certain label can be excluded from the analysis (this is referred to as <u>protein-level filtering</u>).</span></span></td>
							<td><input name="expquantfilt" id="expquantfilt" type="checkbox"/></td>
						</tr>
						<tr class='hidden'>
							<td>&#8212 Filtering based on which label?</td>
							<td><span class="tooltip2">?<span><img class="callout2" src="images/callout_black.gif" />The singlet label.</span></span></td>
							<td><select name="expquantfiltlbl" id = "expquantfiltlbl" onclick="onexpquantfiltlblclick()"></select></td>
						</tr>
						<tr class='hidden'>
							<td>&#8212 Type of filtering: Peptide-level (as opposed to protein-level)?</td>
							<td><input name="expquantfiltprot" type="checkbox"/></td>
						</tr>
						<!--
						<tr>
							<td>&#8212 <u class = "astyle" id="s3showcondsadvanced">Advanced conditions options...</u>
						</tr>
						-->
					</table>
					</div>
				</div>
				<div class="buttonbar">
					<button type="button" class="main" id="s3btnb">Back</button>
					<button type="button" class="main" id="s3btnf">Submit</button>
				</div>
			</div>
			<div class="main_div main_section hidden" id="s4div">
				<h2>Step 4</h2>
				<p id="results_p" style="color: #000000"></p>
				<br>
				<hr class="style-one">
				<br>
				<div id="server_feedback" class="image_array_div">
				</div>
				<div class="buttonbar">
					<button type="button" class="main" id="s4btnb">Back</button>
					<button type="button" class="main" id="s4btnf" disabled>Next</button>
				</div>
			</div>						
			<div class="main_div main_section hidden" id="s5div">
				<h2>End</h2>
				<p style="color: #000000">Please download your results. Click <em>Reset</em> to start again. Thank you for using ProteoSign.</p>
				<div id="s5dndres">
					<br>
					<hr class="style-one">
					<br>
					<div id="dndres" style="text-align: center"></div>
				</div>				
				<div class="buttonbar">
					<button type="button" class="main" id="s5btnb">Back</button>
					<button type="button" class="main" id="s5btnf">Reset</button>
				</div>
			</div>						
			<div class="main_div main_footer">
				<p>
					<em>Reference:</em><br>
					<a href="#" target="_blank"><strong>Efstathiou G., Antonakis A. N., Trudgian D. C., Aivaliotis M., Thomas B., Pavlopoulos G., Papanikolaou N., Acuto O. and Iliopoulos I.<br>ProteoSign: An end-user online differential proteomics statistical analysis service.</strong><br></a>
				<p>
				<span id="proteosignversion" style='font-family: "Courier New"; display: block; float: left'></span>
				<span id="githubproject" style='font-family: "Courier New"; margin: 0 18%;'><a href="https://github.com/yorgodillo/ProteoSign" target="_blank">Project hosted at GitHub <img src="images/github_icon.png" height="30" style="position: relative; vertical-align: top; margin: -0.7em 0;"></img></a></span>
				<span style='font-family: "Courier New"; display: block; float: right'><a id='contact' onclick='$("#contact").attr("href","mailto:msdiffexp@gmail.com?Subject=Session%20"+sessionid); return true;' target='_blank'><u>Contact</u></a></span>
			</div>
		<div>

		<div class="expparamsDlg" id="s3expparamsDlgLabels">
			<div style="text-align: center">
				<span style="font-size: 80%">Select one or more (by holding the <em>Ctrl</em> key).</span>
				<select id="s3expparamsDlgLabelsSelection" style="min-width:90%" multiple>
				</select>
			</div>
			<div class="buttonbar" style="text-align: center; position: relative; width: auto; margin-top: 1em">
				<button type="button" class="main inline" style="width: 4.2em" id="dlglabelsBtnOK">OK</button>
				<button type="button" class="main inline" style="width: 4.2em; margin-left: 1em" id="dlglabelsBtnCancel">Cancel</button>
			</div>			
		</div>
		
		<div class="expparamsDlg" id="s1TestDatasets">
			<div style="text-align: center;">
				<span style="font-size: 85%;">Select a dataset and press <strong>OK</strong>.<br>You can also <strong>download</strong> it to your machine.<br></span>
				<select id="s1TestDatasetsSelection" style="min-width:90%">
				</select>
			</div>
			<div class="buttonbar" style="text-align: center; position: relative; width: auto; margin-top: 1em">
				<button type="button" class="main inline" style="font-size: 80%;" id="dlgTestDatasetsBtnOK">OK</button>
				<button type="button" class="main inline" style="font-size: 80%; margin-left: 1em" id="dlgTestDatasetsBtnCancel">Cancel</button>
				<button type="button" class="main inline unrelated" style="font-size: 80%; margin-left: 1em; " id="dlgTestDatasetsBtnDownload">Download</button>
			</div>			
		</div>		

		<div class="expparamsDlg" id="s2LFQConditions" style="width: 340px">
			<div style="text-align: center;">
				<span style="font-size: 85%;">Select a condition to assign to the selected raw files<br> or type a <strong>New condition</strong> in the box below.</span>
				<input id="s2LFQConditionsNew" placeholder="New Condition" style="width:70%" onkeypress="inputChCheck(event,'^(?!_)[a-zA-Z0-9_]+$',10);" onkeyup="SwitchLFQList(event);"></input>
				<select id="s2LFQConditionsList" style="width:70%">
				</select><br>
			</div>
			<div class="buttonbar" style="text-align: center; position: relative; width: auto; margin-top: 1em">
				<button type="button" class="main inline" style="font-size: 80%;" id="s2LFQConditionsOK" onclick="ons2LFQConditionsOK_click();">OK</button>
				<button type="button" class="main inline" style="font-size: 80%; margin-left: 1em" id="s2LFQConditionsCancel">Cancel</button>
			</div>			
		</div>	
		
		<div class="expparamsDlg" id="s4UserInfo" style="width: 600px">
			<div style="text-align: center;">
				<span style="font-size: 100%;"><strong>ProteoSign</strong> says:</span>
				<div style="text-align: left; font-size: 85%;" id= "s4UserInfoText"></div>
				</select><br>
			</div>
			<div class="buttonbar" style="text-align: center; position: relative; width: auto; margin-top: 1em">
				<button type="button" class="main inline" style="font-size: 80%;" id="s4UserInfoOK">Discard</button>
			</div>			
		</div>
		
		<div class="expparamsDlg" id="condsadvancedDlg" style="width: 600px">
			<div style="text-align: center;">
				<span style="font-size: 85%;">Please type the Name of the <strong>New condition:</strong></span>
				<input id="s3AdvNewCondition" placeholder="New Condition" style="width:70%" onkeypress="inputChCheck(event,'^(?!_)[a-zA-Z0-9_]+$',10);"></input>
			</div>
			<div class="buttonbar" style="text-align: center; position: relative; width: auto; margin-top: 1em">
			<button type="button" class="main inline" style="font-size: 80%;" id="s3AdvancedOK">OK</button>
			<button type="button" class="main inline" style="font-size: 80%; margin-left: 1em" id="s3AdvancedCancel">Cancel</button>
			</div>	
		</div>	
		<script src="http://code.jquery.com/jquery-2.1.0.js"></script>
		<script src="js/jquery.dataTables.min.js"></script>
 		<script src="js/jquery-ui.min.js"></script>
		<script src="js/jquery.ui-contextmenu.min.js"></script>
		<script src="js/dataTables.select.min.js"></script>
		<script src="js/script.js"></script>
</body>
</html>
