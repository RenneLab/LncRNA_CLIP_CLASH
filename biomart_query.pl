#!/usr/bin/env perl
# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5
# use strict;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

my $confFile = "registry.xml";
#my $confFile = "PATH TO YOUR REGISTRY FILE UNDER biomart-perl/conf/. For Biomart Central Registry navigate to
#						http://www.biomart.org/biomart/martservice?type=registry";
#
# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#

@chr = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

foreach (@chr)

{
#OM: Only use 'clean' for re-downloading the registry.
#my $action='clean';
my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

		
	$query->setDataset("hsapiens_gene_ensembl");
	$query->addFilter("chromosome_name", [$_]);
	$query->addFilter("transcript_biotype", ["3prime_overlapping_ncRNA","antisense","bidirectional_promoter_lncRNA","IG_C_pseudogene","IG_J_pseudogene","IG_pseudogene","IG_V_pseudogene","lincRNA","macro_lncRNA","misc_RNA","nonsense_mediated_decay","non_coding","non_stop_decay","polymorphic_pseudogene","processed_pseudogene","processed_transcript","pseudogene","retained_intron","sense_intronic","sense_overlapping","transcribed_processed_pseudogene","transcribed_unitary_pseudogene","transcribed_unprocessed_pseudogene","TR_J_pseudogene","TR_V_pseudogene","unitary_pseudogene","unprocessed_pseudogene"]);
	$query->addFilter("biotype", ["3prime_overlapping_ncRNA","antisense","bidirectional_promoter_lncRNA","IG_C_pseudogene","IG_J_pseudogene","IG_pseudogene","IG_V_pseudogene","lincRNA","macro_lncRNA","misc_RNA","non_coding","polymorphic_pseudogene","processed_pseudogene","processed_transcript","pseudogene","sense_intronic","sense_overlapping","transcribed_processed_pseudogene","transcribed_unitary_pseudogene","transcribed_unprocessed_pseudogene","TR_J_pseudogene","TR_V_pseudogene","unitary_pseudogene","unprocessed_pseudogene"]);
	$query->addAttribute("ensembl_gene_id");
	$query->addAttribute("ensembl_transcript_id");
	#$query->addAttribute("cdna");
	$query->addAttribute("external_gene_name");
	$query->addAttribute("description");
	$query->addAttribute("start_position");
	$query->addAttribute("end_position");
	$query->addAttribute("gene_biotype");
	$query->addAttribute("chromosome_name");
	$query->addAttribute("strand");
	$query->addAttribute("transcript_biotype");
	$query->addAttribute("transcript_start");
	$query->addAttribute("transcript_end");
	$query->addAttribute("transcript_length");
	$query->addAttribute("transcription_start_site");

$query->formatter("FASTA");

my $query_runner = BioMart::QueryRunner->new();
############################## GET COUNT ############################
# $query->count(1);
# $query_runner->execute($query);
# print $query_runner->getCount();
#####################################################################


############################## GET RESULTS ##########################
# to obtain unique rows only
# $query_runner->uniqueRowsOnly(1);

$query_runner->execute($query);
$query_runner->printHeader();
$query_runner->printResults();
$query_runner->printFooter();
#####################################################################

}
