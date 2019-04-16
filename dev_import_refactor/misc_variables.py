




#Definitions for different run mode types

#Auto-annotations
run_mode_pecaan_dict = {\
	'use_basename':'no',\
	'custom_gene_id':'no',\
	'ignore_gene_id_typo':'no',\
	'ignore_description_field_check':'no',\
	'ignore_replace_warning':'no',\
	'ignore_trna_check':'yes',\
	'ignore_locus_tag_import':'yes',\
	'ignore_phage_name_typos':'yes',\
	'ignore_host_typos':'yes',\
	'ignore_generic_author':'yes',\
	'ignore_description_check':'yes'\
	}

#Manual annotations
run_mode_phagesdb_dict = {\
	'use_basename':'no',\
	'custom_gene_id':'no',\
	'ignore_gene_id_typo':'no',\
	'ignore_description_field_check':'no',\
	'ignore_replace_warning':'no',\
	'ignore_trna_check':'no',\
	'ignore_locus_tag_import':'yes',\
	'ignore_phage_name_typos':'no',\
	'ignore_host_typos':'no',\
	'ignore_generic_author':'no',\
	'ignore_description_check':'no'\
	}


#SEA-PHAGES NCBI records
run_mode_ncbi_auto_dict = {\
	'use_basename':'no',\
	'custom_gene_id':'no',\
	'ignore_gene_id_typo':'yes',\
	'ignore_description_field_check':'yes',\
	'ignore_replace_warning':'yes',\
	'ignore_trna_check':'yes',\
	'ignore_locus_tag_import':'no',\
	'ignore_phage_name_typos':'yes',\
	'ignore_host_typos':'no',\
	'ignore_generic_author':'yes',\
	'ignore_description_check':'yes'\
	}


#Misc NCBI records
run_mode_ncbi_misc_dict = {
	'use_basename':'yes',\
	'custom_gene_id':'yes',\
	'ignore_gene_id_typo':'yes',\
	'ignore_description_field_check':'no',\
	'ignore_replace_warning':'yes',\
	'ignore_trna_check':'yes',\
	'ignore_locus_tag_import':'no',\
	'ignore_phage_name_typos':'yes',\
	'ignore_host_typos':'yes',\
	'ignore_generic_author':'yes',\
	'ignore_description_check':'yes'\
	}




#Custom QC settings. User can select the settings, so it is initialized as
#an empty dictionary that only gets filled if there is a ticket indicating
#a custom set of parameters is needed.
run_mode_custom_dict = {}


#A dictionary that holds all the other dictionaries.
#Import tables will use the keys to retrieve the right combination of parameters.
#If new options needed to be created, they need to be added to this dictionary.
#'none': reserved for import tickets that do not need a run mode specified (such as UPDATE tickets)
#'other': reserved for when users manually create import tickets and do not know
#which is the best option for their needs. Currently, this defaults to the 'phagesdb'
#run mode, since that is the most stringest criteria.
#'custom': reserved for when the user wants to specify a unique combination of options
#that are not reflected in the other run modes.
run_mode_options_dict = {\
	'none':'',\
	'pecaan':run_mode_pecaan_dict,\
	'phagesdb':run_mode_phagesdb_dict,\
	'ncbi_auto':run_mode_ncbi_auto_dict,\
	'ncbi_misc':run_mode_ncbi_misc_dict,\
	'other':run_mode_phagesdb_dict,\
	'custom':run_mode_custom_dict}







if run_mode_custom_total > 0:

	print "\n\nOne or more import tickets indicate custom QC parameters."
	print "The following options will be applied to all tickets requiring a custom QC:"

	#1. Use the file's basename as the PhageID instead of the phage name in the file
	print "\n\n\n\nNormally, the PhageID is determined from the Organism field in the record."
	run_mode_custom_dict['use_basename'] = select_option(\
		"\nInstead, do you want to use the file name as the PhageID? (yes or no) ", \
		set(['yes','y','no','n']))

	#2. Create GeneIDs from locus tags or PhageID_GeneNumber concatenation
	#Many times non-Hatfull authored genomes do not have consistent or specific locus tags, which complicates
	#assigning GeneIDs. This option provides a locus tag override and will assign GeneIDs
	#by joining PhageID and CDS number.
	print "\n\n\n\nNormally, the GeneIDs are assigned by using the locus tags."
	run_mode_custom_dict['custom_gene_id'] = select_option(\
		"\nInstead, do you want to create GeneIDs by combining the PhageID and gene number? (yes or no) ", \
		set(['yes','y','no','n']))

	#3. Ensure GeneIDs have phage name spelled correctly?
	#New SEA-PHAGES annotated genomes should be check for spelling, but maybe not
	#for other types of genomes.
	print "\n\n\n\nNormally, the GeneIDs are required to contain the phage name without typos."
	run_mode_custom_dict['ignore_gene_id_typo'] = select_option(\
		"\nInstead, do you want to allow missing or mispelled phage names in the GeneID? (yes or no) ", \
		set(['yes','y','no','n']))

	#4. Gene descriptions are set to import table field regardless
	#This should be run for new SEA-PHAGES annotated genomes, but may want to
	#be skipped if importing many genomes from NCBI
	print "\n\n\n\nNormally, the gene descriptions are verified to be present in the import table qualifier."
	run_mode_custom_dict['ignore_description_field_check'] = select_option(\
		"\nInstead, do you want to use the import table qualifier without verification? (yes or no) ", \
		set(['yes','y','no','n']))

	#5. Replacing final with final warning
	#Once a genome gets in NCBI, it is expected that a Final status genome is
	#replaced with another Final status genome, so it could get annoying to
	#keep getting the warning, so it can be turned off.
	print "\n\n\n\nNormally, a warning is indicated if a Final status genome is being replaced."
	run_mode_custom_dict['ignore_replace_warning'] = select_option(\
		"\nInstead, do you want to silence the status warnings? (yes or no) ", \
		set(['yes','y','no','n']))

	#6. tRNA QC
	#Many genomes from NCBI, including SEA-PHAGES, may not have consistently
	#annotated tRNAs. So the tRNA QC can be skipped.
	print "\n\n\n\nNormally, tRNAs are checked only in new manually annotated genomes."
	run_mode_custom_dict['ignore_trna_check'] = select_option(\
		"\nInstead, do you want to ignore the tRNA quality checks? (yes or no) ", \
		set(['yes','y','no','n']))

	#7. Retain locus tags
	#The locus tag field in the database should only reflect 'official' locus tags from
	#bona fide Genbank records, and not simply Genbank-formatted files.
	#Locus tags from Pecaan auto-annotated and SMART team manually annotated
	#genomes should not be retained.
	print "\n\n\n\nNormally, CDS locus tags are retained only for bona fide Genbank records."
	run_mode_custom_dict['ignore_locus_tag_import'] = select_option(\
		"\nInstead, do you want to ignore all locus tags? (yes or no) ", \
		set(['yes','y','no','n']))

	#8. Ignore phage name typos in the header
	#The phage name can be found in several header fields. This should be
	#checked for new manual annotations. Since NCBI doesn't like to change these
	#types of typos, parsing NCBI records should skip this step.
	print "\n\n\n\nNormally, the phage name is verified only in new manual annotations."
	run_mode_custom_dict['ignore_phage_name_typos'] = select_option(\
		"\nInstead, do you want to ignore any phage name typos in the header? (yes or no) ", \
		set(['yes','y','no','n']))

	#9. Ignore host name typos in the header
	#The host name can be found in several header fields. This should be
	#checked for new manual annotations and in SEA-PHAGES NCBI records.
	print "\n\n\n\nNormally, the host name is verified in new manual annotations or SEA-PHAGES NCBI records."
	run_mode_custom_dict['ignore_host_typos'] = select_option(\
		"\nInstead, do you want to ignore any host name typos in the header? (yes or no) ", \
		set(['yes','y','no','n']))

	#10. Ignore generic author in the author list
	#Sometimes the generic author 'Lastname' or 'Firstname' gets added in DNA Master
	#This should be checked for new manual annotations only.
	print "\n\n\n\nNormally, the author list is checked to ensure the generic author Lastname,Firstname is absent in new manual annotations."
	run_mode_custom_dict['ignore_generic_author'] = select_option(\
		"\nDo you want to ignore any generic authors in the author list? (yes or no) ", \
		set(['yes','y','no','n']))

	#11. Ignore gene description checks
	#The gene description may contain errors.
	#This should be checked for new manual annotations only.
	print "\n\n\n\nNormally, gene descriptions are checked in new manual annotations."
	run_mode_custom_dict['ignore_description_check'] = select_option(\
		"\nDo you want to ignore gene description checks? (yes or no) ", \
		set(['yes','y','no','n']))





	#Now add the customized dictionary of options to the dictionary of all run modes
	run_mode_options_dict['custom'] = run_mode_custom_dict

















#Set up dna and protein alphabets to verify sequence integrity
dna_alphabet_set = set(IUPAC.IUPACUnambiguousDNA.letters)
protein_alphabet_set = set(IUPAC.ExtendedIUPACProtein.letters)


#Create set of all types of actions allowed using this script
#Add = add a new genome without removing another.
#Remove = delete a genome without adding another.
#Replace = delete a genome and replace it with another. Genome names can be different, but the DNA sequence cannot be different.
#Update = make changes to HostStrain, Cluster, Subcluster, or status fields of phages already in the database.
action_set = set(["add","remove","replace","update"])


#Create set of most common gene description genbank qualifiers
description_set = set(["product","note","function"])


#Create list of potential Host Names in the Genbank file to ignore.
#This is primarily for databases that contain phages of all host phyla and not just Actinobacteria
host_ignore = ['enterobacteria','phage','bacteriophage','cyanophage']




#Dictionary for storing authorship info
author_dictionary = {'0':'gbk','1':'hatfull'}


#Phagesdb API to retrieve genome information
api_prefix = "https://phagesdb.org/api/phages/"
api_suffix = "/?format=json"











#phageName typo correction Dictionary
#Key = Phage Name as it is spelled in the Genbank record
#Value = Phage Name as it is spelled in phagesdb and/or Phamerator, and thus
#how it is spelled in the import ticket.
#The phageName parsed from the Genbank record gets reassigned this corrected phageName
#Reasons for the exceptions:

phage_name_typo_dict = {

	#phagesdb is unable to handle underscores
	'ATCC29399B_C':'ATCC29399BC',\
	'ATCC29399B_T':'ATCC29399BT',\

	#ELB20 was reported as 'ELB20' in the original publication, but spelled
	#'phiELB20' in the Genbank record.
	'phiELB20':'ELB20',\


	#Names are spelled differently in Genbank record compared to original publication
	'P100_1':'P100.1',\
	'P100_A':'P100A',\

	#'LeBron' was changed to 'Bron' by ICTV. They won't change it back.
	'Bron':'LeBron',\

	#Inadvertent typos that Genbank won't correct since ICTV uses the typos
	'BBPiebs31':'BPBiebs31',\
	'CaptnMurica':'CapnMurica',\
	'Fionnbarth':'Fionnbharth',\

	#Weird capitalizations in the real phage names that ICTV and now GenBank have "undone"
	'Baka':'BAKA',\
	'CJW1':'Cjw1',\
	'Dlane':'DLane',\
	'Kssjeb':'KSSJEB',\
	'Littlee':'LittleE',\
	'Billknuckles':'BillKnuckles',\
	'Packman':'PackMan',\
	'Mrgordo':'MrGordo',\
	'Ericb':'EricB',\
	'lockley':'Lockley',\
	'Heldan':'HelDan',\
	'Ta17a':'TA17a'\
	}
