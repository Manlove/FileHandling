#!/usr/bin/env python3

# Takes a "mergedregs" file with peak locations and nearest annotated genes,
# a "genes" file with genes and the nearest peaks and outputs a split file with the
# genes from the mergedregs file on seperate lines. Will duplicate peaks regions

from re import split

# input files
file_in_path = "/Users/lmanlove/Documents/1_Active_Projects/4930-01SBGWPharma_Integration_45454/Integration/01SBGWPharm_ATAC_mergedregs.txt"
file_out_path = "/Users/lmanlove/Documents/1_Active_Projects/4930-01SBGWPharma_Integration_45454/Integration/01SBGWPharm_ATAC_split_file.txt"
gene_file_path = "/Users/lmanlove/Documents/1_Active_Projects/4930-01SBGWPharma_Integration_45454/Integration/01SBGWPharm_ATAC_genes.txt"

gene_name_to_symbol = {}

# Steps through the gene file and saves the gene names and matching gene symbols
with open(gene_file_path, 'r') as gene_file:
    gene_file.readline()
    for line in gene_file.readlines():

        line = line.strip().split("\t")
        gene_name_to_symbol[line[1].strip(" ").upper()] = line[0]


# Steps thro
with open(file_in_path, 'r') as file_in:
    with open(file_out_path, 'w') as file_out:

        header = file_in.readline().strip().split("\t")
        link_index = header.index("UCSC Link")

        # Creates a new header from parts of the original header
        # based on the position of the UCSC link column (link_index).
        # ** Some of this can be simplified**
        file_out.write("{}\t{}\t".format("Gene Name", "Gene ID"))
        for i in header[0:link_index-5]:
            file_out.write("{}\t".format(i))
        file_out.write("{}".format("Promoter"))
        file_out.write("\t{}".format(header[link_index-2]))
        file_out.write("\t{}".format(header[link_index-1]))
        for i in header[link_index+1:]:
            file_out.write("\t{}".format(i))
        file_out.write("\n") # <- remove this and add to start of gene writing to unnecessary line at end

        
        for line in file_in.readlines():
            line = line.strip().split("\t")

            # Gene names, distances, and relative positions (upstream, promoter, gene body)
            # are included in order within the columns. **this can be simplified**
            genes = line[link_index - 3].strip('"').split(",")
            distances = line[link_index - 2].strip('"').split(", ")
            positions = line[link_index - 1].strip('"').split(", ")

            # positions = split(r',(?=")', positions)
            if genes[0] != '': # **fix this, can be simplified**
                for gene, distance, position in zip(genes, distances, positions):
                    gene = gene.strip(" ")
                    if gene.upper() in gene_name_to_symbol:
                        file_out.write("{}\t{}\t".format(gene, gene_name_to_symbol[gene.upper()]))
                        for i in line[0:link_index-5]:
                            file_out.write("{}\t".format(i))
                        if int(distance.replace(',','')) >= -7500 and int(distance.replace(',','')) <= 2500:
                            promoter = 1
                        else:
                            promoter = 0
                        file_out.write("{}\t{}\t{}".format(promoter, distance, position))
                        for i in line[link_index+1:]:
                            file_out.write("\t{}".format(i))
                        file_out.write("\n")

if __name__ == "__main__":
    pass
