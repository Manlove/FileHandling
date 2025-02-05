#!/usr/bin/env python3

# Takes a "mergedregs" file with peak locations and nearest annotated genes,
# a "genes" file with genes and the nearest peaks and outputs a split file with the
# genes from the mergedregs file on seperate lines. Will duplicate peaks regions

# input files
file_in_path = "mergedregs.txt"
file_out_path = "split_file.txt"
gene_file_path = "genes.txt"

def GetGeneData(gene_file_path, gene_name_to_symbol):
    '''Steps through the gene file and saves the gene names and matching gene symbols'''
    with open(gene_file_path, 'r') as gene_file:
        next(gene_file)
        for line in gene_file:
            line = line.strip().split("\t")
            gene_name_to_symbol[line[1].strip().upper()] = line[0]

def SplitMergedRegsFile(file_in_path, file_out_path, gene_name_to_symbol):
    '''Splits the merged regions table with single genes on each line'''
    with open(file_in_path, 'r') as file_in:
        with open(file_out_path, 'w') as file_out:

            header = file_in.readline().strip().split("\t")
            link_index = header.index("UCSC Link")

            # Creates a new header from parts of the original header
            # based on the position of the UCSC link column (link_index).
            file_out.write(
                "\t".join(
                    ["Gene Name", "Gene ID"] + 
                    header[0:link_index-5] + 
                    ["Promoter", header[link_index-2], header[link_index-1]] + 
                    header[link_index+1:]
                )
            )

            for line in file_in:
                line = line.strip().split("\t")

                # Gene names, distances, and relative positions (upstream, promoter, gene body)
                # are included in order within the columns.
                genes = line[link_index - 3].strip('"').split(",")
                distances = line[link_index - 2].strip('"').split(", ")
                positions = line[link_index - 1].strip('"').split(", ")

                if genes[0] != '':
                    for gene, distance, position in zip(genes, distances, positions):
                        gene = gene.strip()
                        if gene.upper() in gene_name_to_symbol:

                            if int(distance.replace(',','')) >= -7500 and int(distance.replace(',','')) <= 2500:
                                promoter = 1
                            else:
                                promoter = 0

                            outstring = "\t".join(
                                [gene, gene_name_to_symbol[gene.upper()]] + 
                                line[0:link_index-5] +
                                [str(promoter), distance, position] + 
                                line[link_index+1:]
                            )

                            file_out.write("\n{}".format(outstring))
                        

if __name__ == "__main__":
    gene_name_to_symbol = {}
    GetGeneData(gene_file_path, gene_name_to_symbol)
    SplitMergedRegsFile(file_in_path, file_out_path, gene_name_to_symbol)

