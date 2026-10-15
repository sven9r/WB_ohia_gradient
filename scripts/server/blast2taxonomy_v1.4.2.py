#!/usr/bin/env python3
# blast2taxonomy.py
# Author: Yannis Schöneberg <yannis.schoeneberg@gmx.de>
# This script takes in a blast result table and outputs the taxonomy data in a tsv file
# Version 1.4.2
# Upstream documentation and further details:
# https://github.com/yschoeneberg/blast2taxonomy
import getopt
import sys
import os
import logging
import pandas as pd
from multiprocessing import Pool
from ete3 import NCBITaxa
from itertools import repeat


def get_options(argv):
    global version
    global skip_update
    global skip_failed
    global threads
    global ranks
    global tax_column
    global perc_column
    global len_column
    global fail_file
    version = "1.4.2"
    skip_update = False
    skip_failed = False
    threads = 1
    ranks = ['kingdom', 'phylum', 'superclass', 'class', 'subclass', 'order', 'infraorder', 'superfamily', 'family', 'genus', 'species']
    tax_column = 13
    perc_column = 3
    len_column = 4
    fail_file = "failed_taxids.tsv"

    try:
        opts, args = getopt.getopt(argv, "hfsi:o:c:t:p:l:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print(f"Usage: blast2taxonomy_v{str(version)}.py -i <infile> -o <outfile> -c <column taxids> -t <num threads>\n"
              f"Type blast2taxonomy_v{str(version)}.py -h for help")

    for opt, arg in opts:
        if opt == '-h':
            print(f"\nUsage: blast2taxonomy_v{str(version)}.py [options]\n"
                  f"Version: {str(version)}\n"
                  f"\n"
                  f"REQUIRED:\n"
                  f"\t-i\tTabular Blast results input file\n"
                  f"\t-o\tOutput file\n"
                  f"\n"
                  f"OPTIONAL:\n"
                  f"\t-r\tComma seperated list of taxonomic ranks to extract\n"
                  f"\t-c\tColumn number containing the staxids [13]\n"
                  f"\t-p\tColumn number containing the percent identity [3]\n"
                  f"\t-l\tColumn number containing the length of the subject [4]\n"
                  f"\t-t\tNumber of threads [1]\n"
                  f"\t-s\tSkip Taxonomy Database Update\n"
                  f"\t-f\tSkip failed taxIDs and write those to 'failed_taxids.tsv'"
                  f"\t-h\tDisplay this help message\n"
                  f"\n")
            exit()
        elif opt == '-i':
            global blast_infile
            blast_infile = arg
        elif opt == '-o':
            global outfile
            outfile = arg
        elif opt == "-c":
            tax_column = int(arg)
        elif opt == "-f":
            skip_failed = True
        elif opt == "-t":
            threads = int(arg)
        elif opt == "-s":
            skip_update = True
            try:
                os.remove(fail_file)
            except FileNotFoundError:
                pass
            except Exception as e:
                raise e
        elif opt == "-r":
            ranks = arg.split(",")
        elif opt == "-p":
            perc_column = int(arg)
        elif opt == "-l":
            len_column = int(arg)


def update_taxdb():
    logger.info(f"### Updating Taxonomy Database")
    ncbi.update_taxonomy_database()
    logger.info(f"Finished Updating Database")
    logger.info(f"Removing temporary files")
    os.remove("taxdump.tar.gz")


def get_taxonomy (parameters):
    blast_result = parameters[0]
    ranks = parameters[1]
    tax_column = int(parameters[2])
    perc_column = int(parameters[3])
    len_column = int(parameters[4])
    taxids = blast_result[tax_column-1].split(";")
    percid = blast_result[perc_column-1]
    sbjctlen = blast_result[len_column-1]
    tax_annotations = []
    for id in taxids:
        try:
            lineage = ncbi.get_lineage(id)
        except ValueError as e:
            if skip_failed == True:
                logger.warning(f"Taxid {id} not found! Writing to {fail_file}")
                with open(fail_file, "a") as file:
                    file.write("\t".join(blast_result))
                return
            else:
                if skip_update == True:
                    raise ValueError(f"{id} taxid not found.\n"
                                     f"Skipped internal Taxdb update, try rerunning without the option '-s' to update the taxdb")
                else:
                    raise ValueError(f"{id} taxid not found.\n"
                                 f"Did you just update the NCBI-DB? It might not yet be synchronized with the taxonomy db.\n"
                                 f"Consider trying again later or use the '-f' option to write failed taxids to {fail_file}."
                                 f"Warning: When using '-f' you should check the failed taxids manually!") from e
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        desired_taxids = [ranks2lineage.get(rank, 'Nan') for rank in ranks]
        taxonomy = ncbi.get_taxid_translator(lineage)
        tax_inf = [taxonomy.get(taxid, 'Nan') for taxid in desired_taxids]
        tax_annotations.append(tax_inf)
    transp_taxs = list(zip(*tax_annotations))
    lca_taxonomy = []
    for tax_level in transp_taxs:
        if tax_level.count(tax_level[0]) == len(tax_level):
            lca_taxonomy.append(tax_level[0])
        else:
            break
    return [blast_result[0]] + [percid] + [sbjctlen] + lca_taxonomy


if __name__ == '__main__':
    logger = logging.getLogger('my_logger')
    my_handler = logging.StreamHandler()
    my_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %("
                                     "message)s")
    my_handler.setFormatter(my_formatter)
    logger.addHandler(my_handler)
    logger.setLevel(logging.INFO)

    get_options(sys.argv[1:])

    logger.info(f"#### Extracting Taxonomy Information for Blast Results\n"
                f"{'Program Version:':<50} {str(version)}\n"
                f"{'Blast Results File:':<50} {blast_infile}\n"
                f"{'Column With TaxIDs:':<50} {tax_column}\n"
                f"{'Column With Perc Identity:':<50} {perc_column}\n"
                f"{'Column With Subject Length:':<50} {len_column}\n"
                f"{'Output file:':<50} {outfile}\n"
                f"{'Skip Taxonomy DB update:':<50} {skip_update}\n"
                f"{'Skip failed Taxids:':<50} {skip_failed}\n"
                f"{'Number of threads:':<50} {threads}")
    blast_results = pd.read_csv(blast_infile, sep="\t", dtype=str, header=None)
    blast_results = blast_results.values.tolist()
    global ncbi
    ncbi = NCBITaxa()
    if skip_update is False:
        update_taxdb()
    else:
        logger.info("Skipping Taxonomy Database Update")

    logger.info(f"Searching TaxIDs vs Taxonomy DB")
    with Pool(threads) as pool:
        taxlist = pool.map(get_taxonomy, zip(blast_results,repeat(ranks), repeat(tax_column), repeat(perc_column), repeat(len_column)))
    taxlist = [i for i in taxlist if i is not None]

    logger.info(f"Writing Taxonomy Information to: {outfile}")
    headers = ["query", "perc_id", "sbjct_len"] + ranks
    pd.DataFrame(taxlist, columns=headers).to_csv(outfile, sep="\t", index=False)
    logger.info(f"Done")
