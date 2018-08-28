#!/usr/bin/env python3
# -coregenedir /srv/project_wgmlst/seqSphere_listeria_cgMLST_test/targets/ -inputdir /srv/project_wgmlst/samples_listeria_test  -outputdir /srv/results
import argparse
import sys
import io
import os
import re
import logging
from logging.handlers import RotatingFileHandler

from datetime import datetime
import glob
import pickle
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import Seq

from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from Bio.Blast import NCBIXML




import subprocess
#from subprocess import check_output
import shutil 

def open_log(log_name):
    working_dir = os.getcwd()
    log_name=os.path.join(working_dir, log_name)
    #def create_log ():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    #create the file handler
    handler = logging.handlers.RotatingFileHandler(log_name, maxBytes=40000, backupCount=5)
    handler.setLevel(logging.DEBUG)

    #create a Logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    #add the handlers to the logger
    logger.addHandler(handler)

    return logger

def check_program_is_exec_version (program, version, logger):
    # The function will check if the program is installed in your system and if the version
    # installed matched with the pre-requisites
    if shutil.which(program) is not None :
        # check version
        version_str= str(subprocess.check_output([program, '-version']))
        if not re.search(version, version_str):
            logger.info('%s program does not have the right version ', program)
            print ('Exiting script \n, Version of ' , program, 'does not fulfill the requirements')
            return False
        return True
    else:
        logger.info('Cannot find %s installed on your system', program)
        return False
        

def check_prerequisites (logger):
    pre_requisite_list = [['blastp', '2.6'], ['makeblastdb' , '2.6']]
    # check if blast is installed and has the minimum version 
    for program, version in pre_requisite_list :
        if not check_program_is_exec_version (program , version, logger):
            return False
    return True

def check_arg(args=None):
    
    parser = argparse.ArgumentParser(prog = 'get_subset', description="This program will make the Allele Calling using a predefined core Schema.")
    
    parser.add_argument('-coregenedir', help = 'Directory where the core gene files are located ')
    parser.add_argument('-inputdir', help ='Directory where are located the sample fasta files')
    parser.add_argument('-outputdir', help = 'Directory where the result files will be stored')
    parser.add_argument('-cpus', required= False, help = 'Number of CPUS to be used in the program', default = 3)
    parser.add_argument('-updateschema' , required=False, help = 'Create a new schema with the new locus found', default = True)
    return parser.parse_args()

def is_fasta_file (file_name):
    with open (file_name, 'r') as fh:
        fasta = SeqIO.parse(fh, 'fasta')
        return any(fasta)

def write_first_allele_seq(file_sequence, store_dir, logger):
    #with open (file_name, 'r' ) as fh :
    #seq_record = SeqIO.parse(open(file_name), "genbank").next()
    first_allele_directory = 'first_alleles'
    # split file_sequence into directory and filename
    f_name = os.path.basename(file_sequence)
    full_path_first_allele = os.path.join(store_dir, first_allele_directory)
    if not os.path.exists(full_path_first_allele):
        try:
            os.makedirs(full_path_first_allele)
            logger.info('Directory %s has been created', full_path_first_allele)
        except:
            print ('Cannot create the directory ', full_path_first_allele)
            logger.info('Directory %s cannot be created', full_path_first_allele)
            exit (0)
    first_record = SeqIO.parse(file_sequence, "fasta").__next__()
    # build the fasta file name to store under first_allele_firectory
    fasta_file = os.path.join(full_path_first_allele, f_name)   
    SeqIO.write(first_record, fasta_file, "fasta")
    return True

def create_blastdb (file_name, db_name,db_type, logger ):
    f_name = os.path.basename(file_name).split('.')
    db_dir = os.path.join(db_name,f_name[0])
    output_blast_dir = os.path.join(db_dir, f_name[0])
    if not os.path.exists(db_dir):
        try:
            os.makedirs(db_dir)
            logger.debug(' Created local blast directory for Core Gene %s', f_name[0])
        except:
            logger.info('Cannot create directory for local blast database on Core Gene file %s' , f_name[0])
            print ('Error when creating the directory %s for blastdb. ', db_dir)
            exit(0)
        
        blast_command = ['makeblastdb' , '-in' , file_name , '-parse_seqids', '-dbtype',  db_type, '-out' , output_blast_dir]
        blast_result = subprocess.run(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if blast_result.stderr:
            logger.error('cannot create blast db for %s ', f_name[0])
            logger.error('makeblastdb returning error code %s', blast_result.stderr)
            return False
        
    else:
        logger.info('Skeeping the blastdb creation for %s, as it is already exists', f_name[0])
    return True

def check_blast (reference_allele, sample_files, db_name, logger) :
    for s_file in sample_files:
        f_name = os.path.basename(s_file).split('.')
        dir_name = os.path.dirname(s_file)
        blast_dir = os.path.join(dir_name, db_name,f_name[0])
        blast_db = os.path.join(blast_dir,f_name[0])
        if not os.path.exists(blast_dir) :
            logger.error('Blast db folder for sample %s does not exist', f_name)
            return False
        cline = NcbiblastnCommandline(db=blast_db, evalue=0.001, outfmt=5, max_target_seqs=10, max_hsps=10,num_threads=1, query=reference_allele)
        out, err = cline()
        
        psiblast_xml = StringIO(out)
        blast_records = NCBIXML.parse(psiblast_xml)
        
        for blast_record in blast_records:
            locationcontigs = []
            for alignment in blast_record.alignments:
                # select the best match
                for match in alignment.hsps:
                    alleleMatchid = int((blast_record.query_id.split("_"))[-1])
    return True

def get_fasta_file_list (check_directory,  logger):
    if not os.path.isdir(check_directory):
        logger.info('directory %s does not exists', check_directory)
        return False
    filter_files = os.path.join(check_directory, '*.fasta')
    list_filtered_files =  glob.glob(filter_files)
    list_filtered_files.sort()
    if len (list_filtered_files) == 0 :
        logger.info('directory %s does not have any fasta file ', check_directory)
        return False
    valid_files = []
    for file_name in list_filtered_files:
        if is_fasta_file( file_name):
            valid_files.append(file_name)
        else:
            logger.info('Ignoring file  %s .Does not have a fasta format', file_name)
    if len(valid_files) == 0:
        logger.info('There are not valid fasta files in the directory %s', check_directory)
        logger.debug('Files in the directory are:  $s', list_filtered_files)
        return False
    else:
        return valid_files

def parsing_fasta_file_to_dict (fasta_file, logger):
    fasta_dict = {}
    for contig in SeqIO.parse(fasta_file, "fasta", generic_dna):
            fasta_dict[contig.id] = str(contig.seq.upper())
    logger.debug('file %s parsed to dictionary', fasta_file)
    return fasta_dict

def prepare_core_gene(core_gene_file_list, store_dir, logger):
    #check if directory exists and files have fasta files
    #valid_core_gene_files = get_fasta_file_list(core_gene_dir, logger)
    #if not valid_core_gene_files :
    #    return False
    #logger.debug('Schema files to be processed are : %s', valid_core_gene_files)
    #processing the files in the schema
    file_list = []
    blast_dir = os.path.join(store_dir,'blastdb')
    logger.info('start preparation  of core genes files')
    for fasta_file in core_gene_file_list:
        # parsing fasta file and get in the dictionary the id and the sequence
        fasta_file_parsed_dict = parsing_fasta_file_to_dict(fasta_file, logger)
        f_name = os.path.basename(fasta_file).split('.')
        file_list.append(os.path.join(store_dir, f_name[0]))
        # dump fasta file into pickle file
        with open (file_list[-1],'wb') as f:
            pickle.dump(fasta_file_parsed_dict, f)
        # create the first allele for each core gene file    
        write_first_allele_seq(fasta_file, store_dir, logger)
        # create local blast db for each core gene fasta file

        if not create_blastdb(fasta_file, blast_dir, 'nucl' ,logger):
            print('Error when creating the blastdb for core gene files. Check log file for more information. \n ')
            return False

    return file_list
    
def prepare_samples( sample_file_list, store_dir, logger):
    file_list = []
    blast_dir = os.path.join(store_dir,'blastdb')
    
    for fasta_file in sample_file_list:
        # parsing fasta file and get in the dictionary the id and the sequence
        fasta_file_parsed_dict = parsing_fasta_file_to_dict(fasta_file, logger)
        f_name = os.path.basename(fasta_file).split('.')
        file_list.append(os.path.join(store_dir, f_name[0]))
        # dump fasta file into pickle file
        with open (file_list[-1],'wb') as f:
            pickle.dump(fasta_file_parsed_dict, f)
        
        # create local blast db for each core gene fasta file
        if not create_blastdb(fasta_file, blast_dir, 'nucl' ,logger):
            print('Error when creating the blastdb for core gene files. Check log file for more information. \n ')
            return False

    return file_list
    


def save_files_into_tmp(valid_schema_files, logger):
    
    return True

def allele_call_nucleotides ( core_gene_dict_files, reference_query_directory,  sample_dict_files, blast_db_directory, outputdir, cpus , logger ):
    full_gene_list = []
    samples_matrix_dict = {}
    for core_file in core_gene_dict_files:
        print ( 'core file is : ', core_file)
        full_gene_list.append(os.path.basename(core_file))
        logger.info('Processing core gene file %s ', core_file)
        reference_query = os.path.join(reference_query_directory, str(os.path.basename(core_file) + '.fasta'))
        with open (core_file, 'rb') as core_f:
            core_dict = pickle.load(core_f)
        logger.debug('load in memory the core file %s ', core_file)
        #create new_allele_dict to infer
        for sample_file in sample_dict_files:
            print('sample file is: ', sample_file)
            with open (sample_file,'rb') as sample_f :
                sample_dict = pickle.load(sample_f)
            logger.debug('loaded in memory the sample file %s' , sample_file)
            
            sample_value = os.path.basename(sample_file)
            #intersection = set(core_dict.values()).intersection(gene_dict.values())
            blast_db_name = os.path.join(blast_db_directory, os.path.basename(sample_file),os.path.basename(sample_file))
            #blast_db_name = '/srv/project_wgmlst/samples_listeria/RA-L2073/blastdb'
            #reference_query = '/srv/project_wgmlst/lmo_test.fasta'
            cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 90, outfmt=5, max_target_seqs=10, max_hsps=10,num_threads=1, query=reference_query)
            #cline = NcbiblastnCommandline(db=Gene_Blast_DB_name, evalue=0.001, outfmt=5, max_target_seqs=10, max_hsps=10,num_threads=1, query='/srv/project_wgmlst/seqSphere_listeria_cgMLST_test/targets/lmo0001.fasta')
            out, err = cline()
            psiblast_xml = StringIO(out)
            blast_records = NCBIXML.parse(psiblast_xml)
            score_list = []
            #counter_blast = 0
            for blast_record in blast_records:
                locationcontigs = []
                #import pdb; pdb.set_trace()
                for alignment in blast_record.alignments:
                    # select the best match
                    #print('lenght of found alignments = ', len(alignment.hsps))
                    for match in alignment.hsps:
                        #counter_blast+=1
                        #print ('counter blast is : ', counter_blast)
                        #print ('aligment is : ', alignment, 'match is : ',match,'\n')
                        #print ('sequence is : ', match.sbjct)
                        #print ('sequence start at  : ', match.sbjct_start, ' sequence_end at: ', match.sbjct_end)
                        #print ('score is : ', match.score)
                        #import pdb; pdb.set_trace()
                        score_list.append(match.score)
                        # replace  the "-" character that blast insert
                        sequence_found = match.sbjct.replace('-','')
                        
            intersection = set([sequence_found]).intersection(core_dict.values())
            if len (intersection) > 1:
                print( 'Intersection lenght is : ', len(intersection))
                if sample_value in samples_matrix_dict:
                    samples_matrix_dict[sample_value].append('NIPHEM')
                else:
                    samples_matrix_dict[sample_value] = ['NIPHEM'] 
                logger.debug('NIPHEM was set for core gene %s and sample %s',  os.path.basename(core_file), sample_value)
            elif len(intersection) == 1:
                print( 'Intersection lenght is : ', len(intersection))
                # check if there is any paralog in the sample
                multiple_intersection = []
                for contig_key, contig_value in sample_dict.items():
                    for found_seq in re.finditer(sequence_found,contig_value):
                        multiple_intersection.append([contig_key, found_seq.start(), found_seq.end()])
                
                if len(multiple_intersection) > 1:
                    if sample_value in samples_matrix_dict:
                        samples_matrix_dict[sample_value].append('NIPHEM')
                    else:
                        samples_matrix_dict[sample_value] = ['NIPHEM'] 
                    logger.info('Found several exact match . Possible paralog for core gene %s and sample %s',  os.path.basename(core_file), sample_value)

                    print ('found  paralog: length is : ', multiple_intersection)
                else:
                    print ('Exact match  found : ')
                    for key, values in core_dict.items():
                        if values == sequence_found:
                            allele_match = key
                            break
                    if sample_value in samples_matrix_dict:
                        samples_matrix_dict[sample_value].append(allele_match)
                    else:
                        samples_matrix_dict[sample_value] = [allele_match] 
                    logger.info('Found exact match for core gene %s and sample %s',  os.path.basename(core_file), sample_value)

            else:
                # check
                print( 'Intersection lenght is : ', len(intersection))
                logger.info( 'There is no exact match for core gene %s and sample %s',  os.path.basename(core_file), sample_value)
                logger.debug('Adding new allele to new_allele_dict')
                
                
                '''
                # selecting blastdb from gene
                blastdb_gene = os.path.join(os.path.dirname(core_file),'blastdb',os.path.basename(core_file),os.path.basename(core_file))
                # make blast to find out the allele which match wit 90 % of indetity and 100% alignment
                #cline = NcbiblastnCommandline(db=blastdb_gene, evalue=0.001, perc_identity = 90, outfmt=5, max_target_seqs=10, max_hsps=10,num_threads=1, query=sequence_found)
                cline = NcbiblastnCommandline(db=blastdb_gene, evalue=0.001, perc_identity = 90, outfmt=5, max_target_seqs=10, max_hsps=10,num_threads=1)
                out, err = cline(stdin = sequence_found)
                psiblast_xml = StringIO(out)
                blast_records = NCBIXML.parse(psiblast_xml)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for match in alignment.hsps:
                            print ( 'allele is : ', allele_match , ' match  is: ' , match, )
                            
                length_gene = len(core_dict[0])
                length_sample = len(sequence_found)
                
                
                if sample_value in samples_matrix_dict:
                    samples_matrix_dict[sample_value].append(allele_match)
                else:
                    samples_matrix_dict[sample_value] = [allele_match] 
                logger.info('Found exact match for core gene %s and sample %s',  os.path.basename(core_file), sample_value)
            #import pdb; pdb.set_trace()
            '''
        '''   
        s_name = os.path.basename(schema_file)
        reference_allele = os.path.join(schemadir,'first_alleles', s_name)
        # get the list of the contigs for 
        #
        #
        fullAlleleList = []
        fullAlleleNameList = []
        alleleI = 0


        # check if the reference allele is include in the samples contigs
        #
        #
        check_blast ( reference_allele, sample_files, db_name, logger)
        '''
    print ( 'valor de retorno ', samples_matrix_dict)
    result_file = os.path.join ( outputdir, 'result.tsv')
    #saving the reult information to file
    with open (result_file, 'w') as out_fh:
        out_fh.write ('\t'+'\t'.join( full_gene_list) + '\n')
        for key in sorted (samples_matrix_dict):
            out_fh.write (key + '\t' + '\t'.join(samples_matrix_dict[key])+ '\n')
    return True

if __name__ == '__main__' :
    version = ' wgMLST  0.0.1'
    if sys.argv[1] == '-v' or sys.argv[1] == '--version':
        print( version, '\n')
        exit (0)
    arguments = check_arg(sys.argv[1:])
    start_time = datetime.now()
    # open log file
    logger = open_log ('wgMLST.log')
    # check additional programs are installed in your system
    if not check_prerequisites (logger):
        print ('your system does not fulfill the pre-requistes to run the script ')
        exit(0)
    ##############################################
    # Check that directories contatin fasta files
    ##############################################
    valid_core_gene_files = get_fasta_file_list(arguments.coregenedir, logger)
    if not valid_core_gene_files :
        print ('There are not valid  fasta files in ',  arguments.coregenedir , ' directory. Check log file for more information ')
        exit(0)
    
    valid_sample_files = get_fasta_file_list(arguments.inputdir, logger)
    if not valid_sample_files :
        print ('There are not valid  fasta files in ',  arguments.inputdir , ' directory. Check log file for more information ')
        exit(0)
    ###############################
    # Prepare the coreMLST schema .
    ###############################
    tmp_core_gene_dir = os.path.join(arguments.outputdir,'tmp','cgMLST')
    try:
        os.makedirs(tmp_core_gene_dir)
    except:
        logger.info('Deleting the temporary directory for a previous execution without cleaning up')
        shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
        try:
            os.makedirs(tmp_core_gene_dir)
            logger.info ( 'Temporary folder %s  has been created again', tmp_core_gene_dir)
        except:
            logger.info('Unable to create again the temporary directory %s', tmp_core_gene_dir)
            print('Cannot create temporary directory on ', tmp_core_gene_dir)
            exit(0)

    core_gene_dict_files = prepare_core_gene (valid_core_gene_files , tmp_core_gene_dir , logger)
    if not core_gene_dict_files :
        print('There is an error while processing the schema preparation phase. Check the log file to get more information \n')
        logger.info('Deleting the temporary directory to clean up the temporary files created')
        shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
        exit(0)
    
    #######################################################
    # Prepare the samples files
    #######################################################
    tmp_samples_dir = os.path.join(arguments.outputdir,'tmp','samples')
    try:
        os.makedirs(tmp_samples_dir)
    except:
        logger.info('Deleting the temporary directory for a previous execution without properly cleaning up')
        shutil.rmtree(tmp_samples_dir)
        try:
            os.makedirs(tmp_samples_dir)
            logger.info ( 'Temporary folder %s  has been created again', tmp_samples_dir)
        except:
            logger.info('Unable to create again the temporary directory %s', tmp_samples_dir)
            shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
            logger.info('Cleaned up temporary directory ', )
            print('Cannot create temporary directory on ', tmp_samples_dir, 'Check the log file to get more information \n')
            exit(0)
    sample_dict_files = prepare_samples (valid_sample_files, tmp_samples_dir, logger)
    if not sample_dict_files :
        print('There is an error while processing the saving temporary files. Check the log file to get more information \n')
        logger.info('Deleting the temporary directory to clean up the temporary files created')
        shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
        exit(0)
    
    
    
    reference_query_directory = os.path.join(tmp_core_gene_dir,'first_alleles')
    blast_db_directory = os.path.join(tmp_samples_dir,'blastdb')
    if not allele_call_nucleotides( core_gene_dict_files, reference_query_directory, sample_dict_files,  blast_db_directory, arguments.outputdir,  arguments.cpus , logger):
        print('There is an error while processing the allele calling. Check the log file to get more information \n')
        exit(0)



print ('script ends ')
 